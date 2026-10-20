#!/bin/bash
# Whole-brain LME analysis using AFNI's 3dLMEr.
# Refactored for nested looping and dynamic subbrick extraction.

# =============== CONFIGURATION ===============
RUN_LMER=1           
RUN_SUMMARY=1        
# ==============================================

# Note: The array name is zstat_files, but the items are cope files. 
# Ensure your Python data tables are named lmer_dataTable_cope9.txt, etc.
zstat_files=(cope9 cope10 cope11 cope12)
clinical_vars=("BDI_total")

mask="/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/rois/mni/all_rois_mask.nii.gz"

# Define the root results directory to ensure consistency across model and summary blocks
RESULTS_ROOT="/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/glmm_analysis/Results/univariate_cope"

if [ "$RUN_LMER" -eq 1 ]; then
    for zstat in "${zstat_files[@]}"; do
        for var in "${clinical_vars[@]}"; do
            
            # Define the path to the data table
            DATA_TABLE=/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/glmm_analysis/data_tables/univariate/lmer_dataTable_${zstat}.txt
            
            output="$RESULTS_ROOT/${zstat}/${var}_redo"
            mkdir -p "$output"
            
            echo "=========================================================="
            echo " Running LME Model for: $var ($zstat) "
            echo "=========================================================="
            
            3dLMEr -prefix "$output/lmer_output.nii.gz" \
                   -jobs 8 \
                   -model "${var}+age+sex+(1|scanner)+(1|equipment)" \
                   -qVars "age,${var}" \
                   -mask "${mask}" \
                   -gltCode ${var}_main "${var} :" \
                   -dataTable @$DATA_TABLE
            
            # --- Dynamic Subbrick Extraction ---
            TARGET_CONTRAST="${var}_main"
            
            echo "Locating subbricks for: $TARGET_CONTRAST"
            
            BETA_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST}" "$output/lmer_output.nii.gz")
            STAT_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST} Z" "$output/lmer_output.nii.gz")
            
            # Fallback for t-statistic if Z is not generated
            if [ "$STAT_INDEX" -lt 0 ]; then
                STAT_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST} t" "$output/lmer_output.nii.gz")
            fi
            
            if [ "$STAT_INDEX" -lt 0 ] || [ "$BETA_INDEX" -lt 0 ]; then
                echo "❌ ERROR: Could not find subbricks for $TARGET_CONTRAST."
                continue 
            fi
            
            echo "✅ Found Beta at [$BETA_INDEX] and Stat at [$STAT_INDEX]"
            
            # Extract subbricks dynamically
            3dbucket -prefix "$output/beta.nii.gz" "$output/lmer_output.nii.gz[$BETA_INDEX]"
            3dbucket -prefix "$output/zstat.nii.gz" "$output/lmer_output.nii.gz[$STAT_INDEX]"
            3dcalc -a "$output/lmer_output.nii.gz[$STAT_INDEX]" -expr '-a' -prefix "$output/revzstat.nii.gz"
            
            # Estimate smoothness
            3dFWHMx -mask "${mask}" \
                    -input "$output/lmer_output.nii.gz[$STAT_INDEX]" \
                    -acf "$output/acf.1D" > "$output/acf_params.txt"
            
            ACF_PARAMS=$(sed -n '2p' "$output/acf_params.txt" | awk '{print $1, $2, $3}')
            
            3dClustSim -mask "${mask}" \
                       -acf $ACF_PARAMS \
                       -prefix "$output/ClustSim"
            
            # Apply cluster thresholding
            CLUST_THRESH=$(awk '$1 == 0.001000 {print $3}' "$output/ClustSim.NN2_bisided.1D")
            
            3dClusterize -inset "$output/lmer_output.nii.gz" \
                         -mask "${mask}" \
                         -ithr $STAT_INDEX \
                         -idat $STAT_INDEX \
                         -NN 2 -bisided p=0.001 \
                         -clust_nvox $CLUST_THRESH \
                         -pref_map "$output/zstat_clustered.nii.gz" \
                         -pref_dat "$output/zstat_clustered_masked.nii.gz" \
                         > "$output/cluster_report.txt"
            
        done
    done
    echo "All analyses and extractions complete!"
fi

if [ "$RUN_SUMMARY" -eq 1 ]; then
    echo "------------------------------------------------------"
    echo "Aggregating significant clusters across all models..."

    SUMMARY_FILE="${RESULTS_ROOT}/Global_Cluster_Report.txt"

    printf "Zstat\tModel_Var\tVolume\tCM_RL\tCM_AP\tCM_IS\tminRL\tmaxRL\tminAP\tmaxAP\tminIS\tmaxIS\tMean\tSEM\tMax_Int\tMI_RL\tMI_AP\tMI_IS\n" > "$SUMMARY_FILE"

    for zstat in "${zstat_files[@]}"; do
        for var in "${clinical_vars[@]}"; do
            
            report_path="${RESULTS_ROOT}/${zstat}/${var}/cluster_report.txt"
            
            if [ -f "$report_path" ]; then
                if grep -v '^#' "$report_path" | grep -q '[0-9]'; then
                    grep -v '^#' "$report_path" | \
                    awk -v z="$zstat" -v v="$var" '{print z "\t" v "\t" $0}' >> "$SUMMARY_FILE"
                    
                    echo "Found clusters for: $zstat | $var"
                fi
            else
                echo "Notice: Report file missing for $zstat / $var (likely no significant clusters)"
            fi
            
        done
    done

    echo "Global summary saved to: $SUMMARY_FILE"
    echo "Script complete."
fi

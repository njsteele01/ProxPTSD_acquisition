#!/bin/bash
# Refactored whole-brain LME analysis for testing interactions across conditions.
# Incorporates dynamic subbrick extraction via 3dinfo.

# =============== CONFIGURATION ===============
RUN_LMER=1           
RUN_SUMMARY=1        
# ==============================================

DATA_TABLE=/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/glmm_analysis/data_tables/univariate_int/lmer_dataTable_late.txt
RESULTS_ROOT=/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/glmm_analysis/Results/univariate_cope_int/3way_late

mask="/Volumes/Morey/ProxPTSD.01/Analysis/derivatives/rois/mni/all_rois_mask.nii.gz"
clinical_vars=("CAPS_current" "BDI_total")

if [ "$RUN_LMER" -eq 1 ]; then
    for var in "${clinical_vars[@]}"; do
        
        output="$RESULTS_ROOT/$var"
        mkdir -p "$output"
        
        echo "=========================================================="
        echo " Running Interaction Model for Symptom: $var "
        echo "=========================================================="
        
        # 1. Run the LME Model
        3dLMEr -prefix "$output/lmer_output.nii.gz" \
               -jobs 8 \
               -model "${var}*Threat*Proximity+age+sex+(1|Subj)+(1|scanner)+(1|equipment)" \
               -qVars "age,${var}" \
               -mask "${mask}" \
               -gltCode ${var}_x_Threat_x_Prox "Threat : 1*CSplus -1*CSmin Proximity : 1*Near -1*Far ${var} :" \
               -gltCode ${var}_x_Threat "Threat : 1*CSplus -1*CSmin ${var} :" \
               -gltCode ${var}_x_Prox "Proximity : 1*Near -1*Far ${var} :" \
               -gltCode ${var}_CSplus "Threat : 1*CSplus ${var} :" \
               -gltCode ${var}_CSmin "Threat : 1*CSmin ${var} :" \
               -gltCode ${var}_Near "Proximity : 1*Near ${var} :" \
               -gltCode ${var}_Far "Proximity : 1*Far ${var} :" \
               -dataTable @$DATA_TABLE
        
        # Define an array of the exact strings used in your -gltCode flags
        contrasts=("${var}_x_Threat_x_Prox" "${var}_x_Threat" "${var}_x_Prox")
        
        for TARGET_CONTRAST in "${contrasts[@]}"; do
            
            echo "----------------------------------------------------------"
            echo " Extracting and Thresholding: $TARGET_CONTRAST "
            echo "----------------------------------------------------------"
            
            # 2. Dynamic Subbrick Extraction
            BETA_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST}" "$output/lmer_output.nii.gz")
            STAT_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST} Z" "$output/lmer_output.nii.gz")
            
            # Fallback check
            if [ "$STAT_INDEX" -lt 0 ]; then
                STAT_INDEX=$(3dinfo -label2index "${TARGET_CONTRAST} t" "$output/lmer_output.nii.gz")
            fi
            
            if [ "$STAT_INDEX" -lt 0 ] || [ "$BETA_INDEX" -lt 0 ]; then
                echo "❌ ERROR: Could not find subbricks for $TARGET_CONTRAST."
                continue 
            fi
            
            echo "✅ Found Beta at [$BETA_INDEX] and Stat at [$STAT_INDEX]"
            
            # 3. Extract and Process Subbricks
            3dbucket -prefix "$output/beta_${TARGET_CONTRAST}.nii.gz" "$output/lmer_output.nii.gz[$BETA_INDEX]"
            3dbucket -prefix "$output/zstat_${TARGET_CONTRAST}.nii.gz" "$output/lmer_output.nii.gz[$STAT_INDEX]"
            3dcalc -a "$output/lmer_output.nii.gz[$STAT_INDEX]" -expr '-a' -prefix "$output/revzstat_${TARGET_CONTRAST}.nii.gz"
            
            3dFWHMx -mask "${mask}" \
                    -input "$output/lmer_output.nii.gz[$STAT_INDEX]" \
                    -acf "$output/acf_${TARGET_CONTRAST}.1D" > "$output/acf_params_${TARGET_CONTRAST}.txt"
            
            ACF_PARAMS=$(sed -n '2p' "$output/acf_params_${TARGET_CONTRAST}.txt" | awk '{print $1, $2, $3}')
            
            3dClustSim -mask "${mask}" \
                       -acf $ACF_PARAMS \
                       -prefix "$output/ClustSim_${TARGET_CONTRAST}"
            
            CLUST_THRESH=$(awk '$1 == 0.001000 {print $3}' "$output/ClustSim_${TARGET_CONTRAST}.NN2_bisided.1D")
            
            # 4. Clusterize
            3dClusterize -inset "$output/lmer_output.nii.gz" \
                         -mask "${mask}" \
                         -ithr $STAT_INDEX \
                         -idat $STAT_INDEX \
                         -NN 2 -bisided p=0.001 \
                         -clust_nvox $CLUST_THRESH \
                         -pref_map "$output/zstat_clustered_${TARGET_CONTRAST}.nii.gz" \
                         -pref_dat "$output/zstat_clustered_masked_${TARGET_CONTRAST}.nii.gz" \
                         > "$output/cluster_report_${TARGET_CONTRAST}.txt"
                         
        done
                     
    done
    echo "All 3dLMEr analyses and extractions complete!"
fi

if [ "$RUN_SUMMARY" -eq 1 ]; then
    echo "------------------------------------------------------"
    echo "Aggregating significant clusters across all models..."

    SUMMARY_FILE="${RESULTS_ROOT}/Global_Cluster_Report.txt"
    
    # Added a 'Contrast' column so you know which row belongs to which interaction
    printf "Model_Var\tContrast\tVolume\tCM_RL\tCM_AP\tCM_IS\tminRL\tmaxRL\tminAP\tmaxAP\tminIS\tmaxIS\tMean\tSEM\tMax_Int\tMI_RL\tMI_AP\tMI_IS\n" > "$SUMMARY_FILE"

    for var in "${clinical_vars[@]}"; do
        
        # Define the exact contrast names again to loop through the generated files
        contrasts=("${var}_x_Threat_x_Prox" "${var}_x_Threat" "${var}_x_Prox")
        
        for TARGET_CONTRAST in "${contrasts[@]}"; do
            report_path="${RESULTS_ROOT}/${var}/cluster_report_${TARGET_CONTRAST}.txt"
            
            if [ -f "$report_path" ]; then
                if grep -v '^#' "$report_path" | grep -q '[0-9]'; then
                    # Pass the contrast name into awk to append it as a column
                    grep -v '^#' "$report_path" | \
                    awk -v v="$var" -v c="$TARGET_CONTRAST" '{print v "\t" c "\t" $0}' >> "$SUMMARY_FILE"
                    
                    echo "Found clusters for: $var | $TARGET_CONTRAST"
                fi
            else
                # Optional: stay silent if no clusters found, or keep this echo to verify file paths
                echo "Notice: No report file found/created for $TARGET_CONTRAST (likely no significant clusters)."
            fi
        done
    done

    echo "Global summary saved to: $SUMMARY_FILE"
    echo "Script complete."
fi
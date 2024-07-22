#!/bin/bash

# Define arrays for harmonizations and scales
harmonizations=("combat" "no_harmo") #"min_combat" "no_harmo"
scales=("scaled" "no_scaled") 
time_points=("t1" "t2" "t3")
measures=("area", "sulc", "thk", "vol")

# Paths
path_gcta=/home/cmlr2/rds/hpc-work/ABCD/GCTA
path_grm=${path_gcta}/ABCD_GRM_unrelated
# Path to the ROI names file
# filename=${path_gcta}/roi_names.txt
# Read the ROI names into an array
# -t roi_names < "$filename"

# Initialize an array to store ROI names
roi_names=()

# Read the file into an array using a while loop and read
filename="$path_gcta/roi_names.txt"
while IFS= read -r line; do
    roi_names+=("$line")
done < "$filename"



# Iterate over each harmonization method
for harmo in "${harmonizations[@]}"; do
    
    # Iterate over each scale method
    for scal in "${scales[@]}"; do
        
        directory=${path_gcta}/heritability/${harmo}/${scal}
        
        path_output="$path_gcta/gcta_all_output3/$harmo/$scal/co_heritability"
        mkdir -p "$path_output/files_hsq"
        mkdir -p "$path_output/files_log"
        
        # Loop through the combinations of time points without duplicates.
        for ((i=0; i<${#time_points[@]}; i++)); do
            for ((j=i+1; j<${#time_points[@]}; j++)); do
                tp1=${time_points[i]}
                tp2=${time_points[j]}
                
                # Store the files corresponding to a time point. 
                time_dir1=("$directory"/meas_*_${tp1}.txt)  
                # Find the corresponding "cat" file with the same suffix
                cat_covar=("$directory"/cat_${tp1}.txt) 
                quant_covar=("$directory"/quant_${tp1}.txt)

                # Store the files corresponding to a time point. 
                time_dir2=("$directory"/meas_*_${tp2}.txt)  
                
                # mkdir -p "$path_output/all_roi_results/${tp1}"
                
                # Iterate through the measurements for each time point. 
                for ((l=0; l<${#time_dir1[@]}; l++)); do
                    
                    file1=${time_dir1[l]}
                    #echo $file1
                    echo "~~~~~"

                    file2=${time_dir2[l]}
                    #echo $file2

                    measure=$(echo $file2 | awk -F'/' '{split($NF, a, "_"); print a[2]}')
                    echo $measure

                    time_pt1=$tp1
                    time_pt2=$tp2  # Extract the time point of the file
                    echo " ~~~~~~~~~~  ${time_pt1}_${time_pt2} ~~~~~~~~~~"

                    # Output path
                    
                    output_co_herit="$path_output/all_roi_results/co_herit_${measure}_${time_pt1}${time_pt2}.txt"
                    co_herit_hsq="$path_output/files_hsq/co_herit_${measure}_${time_pt1}${time_pt2}"
                    co_herit_log="$path_output/files_log/co_herit_${measure}_${time_pt}_${time_pt2}"
                    # Create the output directory if it doesn't exist
                    mkdir -p "$path_output/all_roi_results"
                    # Create or overwrite the output file and add a header
                    echo -e "Regions\trG\tSE" > "$output_co_herit"

                    # Create temporary files for the phenotypes
                    output1="${path_output}/output1.txt"
                    output2="${path_output}/output2.txt"
                    awk 'NR==FNR{a[$2]; next} $2 in a' "$file2" "$file1" > "${output1}"
                    awk 'NR==FNR{a[$2]; next} $2 in a' "$file1" "$file2" > "${output2}"

                    # Extract and merge the necessary columns for co-heritability
                    pheno1="${path_output}/co_pheno1.txt"
                    pheno2="${path_output}/co_pheno2.txt"


                    # Loop through the ROIs
                    for ((k=0; k<${#roi_names[@]}; k++)); do
                        echo "ROI Index: $k"
                        
                        # Create or clear the co_pheno.txt file
                        path_co_pheno="$path_output/co_pheno.txt"
                        : > "${path_co_pheno}"

                        awk -v k=$((k+3)) -F" " '{print $1, $2, $k}' "${output1}" > "${pheno1}"
                        awk -v k=$((k+3)) -F" " '{print $k}' "${output2}" > "${pheno2}"
                        paste "${pheno1}" "${pheno2}" > "${path_co_pheno}"

                        # Calculate co-heritability
                        
                            /home/cmlr2/rds/hpc-work/Software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64\
                            --reml-bivar 1 2 \
                            --grm "$path_grm" \
                            --pheno "$path_co_pheno" \
                            --covar "$cat_covar"\
                            --qcovar "$quant_covar"\
                            --out "${co_herit_hsq}_${k}" \
                            --grm-cutoff 0.05 >> $co_herit_log

                        # Extract co-heritability values
                        valueV=$(grep 'rG' "${co_herit_hsq}_${k}.hsq" | awk '{print $2}')
                        valueSE=$(grep 'rG' "${co_herit_hsq}_${k}.hsq" | awk '{print $3}')
                        
                        # Append the ROI name and heritability value to the output file
                        echo -e "${roi_names[$k]}\t${valueV}\t${valueSE}" >> "$output_co_herit"
                    done


                        
                done
                echo "_________________________"
            done    
        done
    done
done


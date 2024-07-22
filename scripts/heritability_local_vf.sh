
# Define arrays for harmonizations and scales
harmonizations=("combat" "no_harmo") #"min_combat" "combat"
scales=("scaled" "no_scaled") 

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
        
        # Directory containing the phenotype files
        directory="$path_gcta/heritability/$harmo/$scal"
        
        # Get the list of files in the directory
        files=($(ls $directory/meas*))
        cat_covs=($(ls $directory/cat*))
        quant_covs=($(ls $directory/quant*))
        
        
        # Output path
        path_output="$path_gcta/gcta_all_output2/$harmo/$scal/heritability"
        # Create the output directory if it doesn't exist
        mkdir -p "$path_output"
        mkdir -p $path_output/files_hsq
        mkdir -p $path_output/files_log
        
        
        
        
        echo "about to start iterating through files"

        # Outer loop: iterate through each file
        for ((i=0; i<${#files[@]}; i++)); do
        
    
            file=${files[i]}  # Get the current file
            
            # Get the base name of the current file without extension
            time_pt=$(basename "$file" | sed 's/\(.*\)\..*/\1/')
            echo $time_pt
            # Define the output file for heritability results
            mkdir -p "$path_output/all_roi_results"
            output_herit="$path_output/all_roi_results/herit_${time_pt}.txt"
            herit_log="$path_output/files_log/herit_${time_pt}"
            
            herit_hsq="$path_output/files_hsq/herit_${time_pt}"
            herit_log="$path_output/files_log/herit_${time_pt}"

            # Create or overwrite the output file and add a header
            echo -e "Regions\tV(G)/Vp\tSE\tPval" > "$output_herit"
            
        
            echo $file
            # Extract the part of the filename after the last underscore
            suffix=$(echo "$file" | sed 's/.*_\(.*\)/\1/')
            
            # Find the corresponding "cat" file with the same suffix
            cat_covar=$(ls $directory/cat_* | grep "_$suffix")
            quant_covar=$(ls $directory/quant_* | grep "_$suffix")
            
            # Loop through the ROIs
            for ((k=0; k<${#roi_names[@]}; k++)); do
                echo "ROI Index: $k"
                
                # Calculate heritability
                herit_hsq="$path_output/files_hsq/herit_${time_pt}"
                /home/cmlr2/rds/hpc-work/Software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64\
                    --grm "$path_grm"\
                    --pheno "$file"\
                    --mpheno $((k+1))\
                    --covar "$cat_covar"\
                    --qcovar "$quant_covar"\
                    --reml --out "$herit_hsq"\
                    --grm-cutoff 0.05 > $herit_log
                
                valueV=$(grep 'V(G)/Vp' "${herit_hsq}.hsq" | awk '{print $2}')
                valueSE=$(grep 'V(G)/Vp' "${herit_hsq}.hsq" | awk '{print $3}')
                Pvalue=$(grep 'Pval' "${herit_hsq}.hsq" | awk '{print $2}')
                
                # Append the ROI name and heritability value to the output file
                echo -e "${roi_names[$k]}\t${valueV}\t${valueSE}\t${Pvalue}" >> "$output_herit"
            
            
            done
        done
    done
done

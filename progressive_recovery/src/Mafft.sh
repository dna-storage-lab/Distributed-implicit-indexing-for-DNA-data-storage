#!/bin/bash

input_folder="$1"             
output_folder="$2"            
threads=48                   

start_time=$(date +%s)
mkdir -p "$output_folder"

fasta_files=("$input_folder"/*.fasta)
if [ ${#fasta_files[@]} -eq 0 ]; then
    echo "ERROR: No fasta files found in input directory!"
    exit 1
fi

echo "Starting parallel processing (threads=$threads)..."
for i in $(seq 1 $threads); do
    mkdir -p "$input_folder/subfolder_$i"
done

for i in "${!fasta_files[@]}"; do
    subfolder_idx=$((i % threads + 1))
    mv "${fasta_files[$i]}" "$input_folder/subfolder_$subfolder_idx/" &
done
wait

process_subfolder() {
    local subfolder="$1"
    local output_folder="$2"
    
    mkdir -p "$output_folder/$(basename "$subfolder")"
    for fasta in "$subfolder"/*.fasta; do
        out_file="$output_folder/$(basename "$subfolder")/$(basename "$fasta")"
        mafft --auto --thread 1 --quiet "$fasta" > "$out_file" 2>/dev/null
    done
}
export -f process_subfolder  

parallel -j "$threads" process_subfolder {} "$output_folder" ::: "$input_folder"/subfolder_*

rm -rf "$input_folder"/subfolder_*

for i in $(seq 1 $threads); do
    mv "$output_folder/subfolder_$i"/*.fasta "$output_folder"/
done


rm -rf "$output_folder"/subfolder_*

end_time=$(date +%s)
runtime=$((end_time - start_time))
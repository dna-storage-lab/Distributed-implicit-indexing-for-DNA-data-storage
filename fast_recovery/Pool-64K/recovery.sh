#!/bin/bash

set -e  
set -u 

export LC_TIME=en_US.UTF-8
echo "=================== Begin Processing ==================="
echo "Start Time: $(date)"
echo "========================================================"

# --- Configuration Parameters ---
THREADS=48
ULIMIT_VAL=655350
TOTAL_CODEWORDS=147
INPUT_FASTQ="data/Pool-64K-NGS.fastq"

COVERAGE_LIST=(1.7) 

# --- Tool Paths ---
BIN_DIR="./bin"
EXTRACT_TOOL="${BIN_DIR}/SubsampleFastq"
PRIMER_TOOL="${BIN_DIR}/PrimerIdentification"
CORR_TOOL="${BIN_DIR}/SlidingCorrelation"
VOTING_TOOL="${BIN_DIR}/MajorityVoting"
DECODER_TOOL="${BIN_DIR}/LDPC_r2_3_decoder_soft"
HAMMING_TOOL="${BIN_DIR}/PostDecHammingDistance"

# --- Initialize Environment ---
ulimit -s $ULIMIT_VAL
mkdir -p ./result

start_time=$(date +%s.%3N)
for COVERAGE in "${COVERAGE_LIST[@]}"; do
    coverage_int=$(echo "$COVERAGE * 10" | bc -l | awk '{print int($1)}')
    L_READS=$((coverage_int * 64827 / 10))

    OUTPUT_DIR="./result/${COVERAGE}x"
    
    # Create output subdirectories
    mkdir -p "${OUTPUT_DIR}/"{PrimerIdentification,Correlation,MajorityVoting,Decoding/RecoveryBitstream}

    echo -e "\n---------------- Coverage ${COVERAGE}x -----------------"

    # ---------------------------------------------------------
    # Step 1: Sub-sampling Reads
    # ---------------------------------------------------------
    echo "[Step 1] Sampling $L_READS reads"

    "$EXTRACT_TOOL" \
        "$INPUT_FASTQ" \
        "$L_READS" \
        "$OUTPUT_DIR"

    # ---------------------------------------------------------
    # Step 2: Primer Identification
    # ---------------------------------------------------------
    echo -e "\n[Step 2] Primer Identification"

    PRIMER_LIST="configure/primer_list.txt"     
    PRIMER_NUM=2                      
    PRIMER_LEN=10                     
    MIN_EDIT=2                        
    MIN_LENGTH=105                   
    MAX_LENGTH=115                   
    OLIGO_LEN=130                     
    MAX_RECURSION_DEPTH=5                                  
    RESULT_PREFIX="${OUTPUT_DIR}/PrimerIdentification"

    "$PRIMER_TOOL" \
        "${OUTPUT_DIR}/sub_reads.fastq" \
        "$PRIMER_LIST" \
        "$PRIMER_NUM" \
        "$PRIMER_LEN" \
        "$MIN_EDIT" \
        "$MIN_LENGTH" \
        "$MAX_LENGTH" \
        "$OLIGO_LEN" \
        "$MAX_RECURSION_DEPTH" \
        "$RESULT_PREFIX" \
        "$THREADS"

    # ---------------------------------------------------------
    # Step 3: Sliding Correlation
    # ---------------------------------------------------------
    echo -e "\n[Step 3] Sliding Correlation"

    SHIFT_RANGE=1
    NUM_GPUS=1

    "$CORR_TOOL" \
        "$COVERAGE" \
        "$SHIFT_RANGE" \
        "$NUM_GPUS"

    # ---------------------------------------------------------
    # Step 4: Majority Voting
    # ---------------------------------------------------------
    echo -e "\n[Step 4] Bit-wise Majority Voting"
    CORR_THRESHOLD=0.40 

    "$VOTING_TOOL" \
        "$COVERAGE" \
        "${OUTPUT_DIR}/MajorityVoting" \
        "$CORR_THRESHOLD" 

    # ---------------------------------------------------------
    # Step 5: LDPC Decoding
    # ---------------------------------------------------------
    echo -e "\n[Step 5] LDPC Decoding"

    DecodedBitstreamFile="${OUTPUT_DIR}/Decoding/recovery_bitstream.txt"
    OutputImageFile="${OUTPUT_DIR}/Decoding/recovery_image.jpg"
    OriginalBitstreamFile="configure/source_bit.txt"
    ErrorRateFile="${OUTPUT_DIR}/Decoding/post_decoding_error.txt"

    for k in $(seq 1 147); do
        Codeword="${OUTPUT_DIR}/MajorityVoting/$k.txt"
        
        # Decoding
        ./bin/LDPC_r2_3_decoder_soft \
            "$Codeword" \
            "$DecodedBitstreamFile" \
            "$OutputImageFile" \
            "S" 
        
        # Hamming distance calculation
        ./bin/PostDecHammingDistance \
            "$DecodedBitstreamFile" \
            "$OriginalBitstreamFile" \
            "$ErrorRateFile" \
            "$k"
        

        percent=$((k * 100 / 147))
        filled=$((percent / 2))      
        empty=$((50 - filled))
        bar=$(printf '%*s' "$filled" | tr ' ' '#')
        space=$(printf '%*s' "$empty" | tr ' ' ' ')
        printf "\r[%-50s] %d%% (%d/147)" "$bar$space" "$percent" "$k"
    done
    echo  

    # --- Summary Statistics ---
    RECOVERED=$(awk '$1==0 {count++} END{print count+0}' "$ErrorRateFile")
    RATE=$(echo "scale=4; $RECOVERED / $TOTAL_CODEWORDS" | bc)

    echo "------------------------------------------------"
    echo "Final Statistics for ${COVERAGE}x:"
    echo "Recovered codeword: ${RECOVERED}/${TOTAL_CODEWORDS}"
    printf "Recovery Success Rate: %.4f\n" "$RATE"
    echo "------------------------------------------------"

    # ---------------------------------------------------------
    # Final Image Reconstruction
    # ---------------------------------------------------------
    if [ "$RECOVERED" -eq "$TOTAL_CODEWORDS" ]; then
        echo "Recovery Success!"
        xdg-open "${OUTPUT_DIR}/Decoding/recovery_image.jpg" >/dev/null 2>&1 &
    fi
done

echo -e "\nEnd Time: $(date)"
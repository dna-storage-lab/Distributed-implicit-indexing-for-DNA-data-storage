#!/bin/bash
set -e
set -u

export LC_TIME=en_US.UTF-8
echo "=================== Begin Processing ==================="
echo "Start Time: $(date)"
echo "========================================================"

# --- Configuration Parameters ---
THREADS=60
NUM_GPUS=1
TOTAL_CODEWORDS=474
ULIMIT_VAL=655350
ulimit -s $ULIMIT_VAL
COVERAGE_LIST=(3.7)

INPUT_FASTQ="data/Pool-120K-ONT-SUP.fastq"

ORIGINAL_BITSTREAM="configure/source_bit.txt"

# --- ALGORITHM HYPERPARAMETERS ---
# Primer Identification 
PRIMER_LIST="configure/primer_list.txt"     
PRIMER_NUM=2
PRIMER_LEN=10
MIN_EDIT=2
MIN_LENGTH=155
MAX_LENGTH=165
OLIGO_LEN=180
MAX_RECURSION_DEPTH=5

# Correlation & FB Parameters 
CORR_THRESHOLD=0.30
SHIFT_RANGE=1
INS_RATE=0.0078
DEL_RATE=0.0087
SUB_RATE=0.0061
FB_MODE=0

# ==============================================================================
# MAIN RECOVERY LOOP
# ==============================================================================
for COVERAGE in "${COVERAGE_LIST[@]}"; do
    COVERAGE_INT=$(echo "$COVERAGE * 10" | bc -l | awk '{print int($1)}')
    L_READS=$((COVERAGE_INT * 119922 / 10))

    for EXP_ID in {1..1}; do
        SUB_DIR="result/${COVERAGE}x"
        STAGE_I_DIR="${SUB_DIR}/StageI_result"
        
        # Initialize Directory Structure 
        mkdir -p "${STAGE_I_DIR}/PrimerIdentification"
        mkdir -p "${STAGE_I_DIR}/Correlation"
        mkdir -p "${STAGE_I_DIR}/FB"
        mkdir -p "${STAGE_I_DIR}/SoftInfo"
        mkdir -p "${STAGE_I_DIR}/Decoding/RecoveryBitstream"

        echo -e "\n---------------- Coverage ${COVERAGE}x -----------------"

        # ---------------------------------------------------------
        # STAGE I RECOVERY 
        # ---------------------------------------------------------
        # [Step I-1] Sub-sampling Reads
        echo "[Step I-1] Sampling $L_READS reads"
        ./bin/SubsampleFastq \
            "$INPUT_FASTQ" \
            "$L_READS" \
            "$EXP_ID" \
            "$SUB_DIR"

        # [Step I-2] Primer Identification
        echo -e "\n[Step I-2] Primer Identification"
        RESULT_PREFIX="${STAGE_I_DIR}/PrimerIdentification"
        ./bin/PrimerIdentification \
            "${SUB_DIR}/sub_reads.fastq" \
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

        # [Step I-3] Sliding Correlation (CUDA Accelerated)
        echo -e "\n[Step I-3] Sliding Correlation"
        CORR_INPUT="${STAGE_I_DIR}/PrimerIdentification/read_without_primer.fastq"
        CORR_OUTPUT="${STAGE_I_DIR}/Correlation/corr_result.txt"
    
        ./bin/SlidingCorrelation \
            "$COVERAGE" \
            "$SHIFT_RANGE" \
            "$NUM_GPUS" \
            "$CORR_INPUT" \
            "$CORR_OUTPUT" 

        # Split Reads by Threshold
        ./bin/FilterThreshold \
            "${STAGE_I_DIR}/Correlation/corr_result.txt" \
            "$CORR_THRESHOLD" \
            "${STAGE_I_DIR}/Correlation/Type-I_reads.txt" \
            "${STAGE_I_DIR}/Correlation/lowthres_reads.fasta"

        # [Step I-4] Forward-Backward Algorithm 
        echo -e "\n[Step I-4] Forward-Backward Algorithm"
        ./bin/FB \
            "${STAGE_I_DIR}/Correlation/Type-I_reads.txt" \
            "${STAGE_I_DIR}/FB/TypeI_symbol_probability.txt" \
            "$INS_RATE" \
            "$SUB_RATE" \
            "$DEL_RATE" \
            "$FB_MODE" \
            "$THREADS"

        # [Step I-5] Multi-read Merging
        echo -e "\n[Step I-5] Multi-read Merging"
        ./bin/MultiReadMerging_StageI \
            "${STAGE_I_DIR}/FB/TypeI_symbol_probability.txt" \
            "${STAGE_I_DIR}/SoftInfo" \
            "$THREADS"

        # [Step I-6] Primary LDPC Decoding 
        echo -e "\n[Step I-6] LDPC Decoding"
        STAGE_I_ERR="${STAGE_I_DIR}/Decoding/post_decoding_error.txt"

        seq 1 $TOTAL_CODEWORDS | parallel -j "$THREADS" \
            "./bin/LDPC_r2_3_decoder_soft ${STAGE_I_DIR}/SoftInfo/{}.txt ${STAGE_I_DIR}/Decoding/RecoveryBitstream/{}.txt S && \
             ./bin/PostDecHammingDistance ${STAGE_I_DIR}/Decoding/RecoveryBitstream/{}.txt $ORIGINAL_BITSTREAM $STAGE_I_ERR {} "

        result=$(awk -v total=474 '$2==0 {count++} END{if(total>0) printf "%d/%d=%.4f", count, total, count/total; else print "0/0=0.0000"}' "$STAGE_I_ERR")
        echo "StageI recovery rate: $result"

        failed_cw=($(awk '$2 > 0 {print $1}' "$STAGE_I_ERR"))
        failed_cw_str=$(IFS=" "; echo "${failed_cw[*]}")
        failed_count=${#failed_cw[@]}
        echo "Total failed Codeword: $failed_count"

        # Verify Stage I Completion
        STAGE_I_SUM=$(awk '{s+=$2} END {print s+0}' "$STAGE_I_ERR")
        if [ "$STAGE_I_SUM" -eq 0 ]; then

            echo "Stage I Recovery Success!"
            ./bin/ImageParsing \
                "${STAGE_I_DIR}/Decoding/RecoveryBitstream" \
                "${STAGE_I_DIR}/Decoding/recovered_image.jpg"

            xdg-open "${STAGE_I_DIR}/Decoding/recovered_image.jpg" >/dev/null 2>&1 &

        else
            echo -e "\nStage I Failed! ➜  Starting Stage II Recovery"

            # ---------------------------------------------------------
            # STAGE II RECOVERY 
            # ---------------------------------------------------------
            STAGE_II_DIR="${SUB_DIR}/StageII_result"
            mkdir -p "${STAGE_II_DIR}/Cluster/cu_more"
            mkdir -p "${STAGE_II_DIR}/Mafft"
            mkdir -p "${STAGE_II_DIR}/Correlation"
            mkdir -p "${STAGE_II_DIR}/FB"
            mkdir -p "${STAGE_II_DIR}/MajorityVoting"
            mkdir -p "${STAGE_II_DIR}/SoftInfo"
            mkdir -p "${STAGE_II_DIR}/Decoding/RecoveryBitstream"

            # Parse Failed Codeword Indices
            FAILED_CW=($(awk '$2 > 0 {print $1}' "$STAGE_I_ERR"))
            FAILED_CW_STR=$(IFS=" "; echo "${FAILED_CW[*]}")

            # [Step II-1] Clustering
            echo -e "\n[Step II-1] Clustering"

            ./bin/usearch_linux_x86_12.0-beta \
                -cluster_fast "${STAGE_I_DIR}/Correlation/lowthres_reads.fasta" \
                -id 0.9 \
                -threads "$THREADS" \
                -uc "${STAGE_II_DIR}/Cluster/Cluster_result.uc" \
                -quiet > /dev/null 2>&1


            # [Step II-2] Multiple Sequence Alignment
            echo -e "\n[Step II-2] Multiple Sequence Alignment"
            ./bin/ClusterSplitter \
                "${STAGE_II_DIR}/Cluster/Cluster_result.uc" \
                "${STAGE_I_DIR}/Correlation/lowthres_reads.fasta" \
                "${STAGE_II_DIR}/Cluster/cu_more"

            ./src/Mafft.sh \
                "${STAGE_II_DIR}/Cluster/cu_more" \
                "${STAGE_II_DIR}/Mafft/mafft_output"
            rm -rf "${STAGE_II_DIR}/Cluster/cu_more"

            ./bin/MajorityVoting \
                -i "${STAGE_II_DIR}/Mafft/mafft_output" \
                -o "${STAGE_II_DIR}/MajorityVoting/majority_voting.fastq"

            # [Step II-3] Sliding Correlation
            echo -e "\n[Step II-3] Sliding Correlation"

            CORR_INPUT="${STAGE_II_DIR}/MajorityVoting/majority_voting.fastq"
            CORR_OUTPUT="${STAGE_II_DIR}/Correlation/corr_result.txt"
        
            ./bin/SlidingCorrelation \
                "$COVERAGE" \
                "$SHIFT_RANGE" \
                "$NUM_GPUS" \
                "$CORR_INPUT" \
                "$CORR_OUTPUT" 

            ./bin/FilterThreshold \
                "${STAGE_II_DIR}/Correlation/corr_result.txt" \
                "$CORR_THRESHOLD" \
                "${STAGE_II_DIR}/Correlation/Type-II_reads.txt" \
                "${STAGE_II_DIR}/Correlation/lowthres_reads.fasta"

            # [Step II-4] Forward-Backward Algorithm 
            echo -e "\n[Step II-4] Forward-Backward Algorithm"
            ./bin/FB \
                "${STAGE_II_DIR}/Correlation/Type-II_reads.txt" \
                "${STAGE_II_DIR}/FB/TypeII_symbol_probability.txt" \
                "$INS_RATE" \
                "$SUB_RATE" \
                "$DEL_RATE" \
                "$FB_MODE" \
                "$THREADS"

            # [Step II-5] Multi-read Merging
            echo -e "\n[Step II-5] Multi-read Merging"
            cat "${STAGE_I_DIR}/FB/TypeI_symbol_probability.txt" \
                "${STAGE_II_DIR}/FB/TypeII_symbol_probability.txt" > "${STAGE_II_DIR}/FB/Symbol_Probability.txt"

            ./bin/MultiReadMerging_StageII \
                "${STAGE_II_DIR}/FB/Symbol_Probability.txt" \
                "${STAGE_II_DIR}/SoftInfo" \
                "$THREADS" \
                "$FAILED_CW_STR"

            # [Step II-6] LDPC Decoding for Remaining Codewords 
            echo -e "\n[Step II-6] LDPC Decoding"
            STAGE_II_ERR="${STAGE_II_DIR}/Decoding/post_decoding_error.txt"
            > "$STAGE_II_ERR"
            for CW_IDX in "${FAILED_CW[@]}"; do

                echo "++++++++++++++++++++++++++++++++++++++++++++" 
                echo "   Codeword $CW_IDX decoding" 
                echo "++++++++++++++++++++++++++++++++++++++++++++"

                ./bin/LDPC_r2_3_decoder_soft \
                    "${STAGE_II_DIR}/SoftInfo/${CW_IDX}.txt" \
                    "${STAGE_II_DIR}/Decoding/RecoveryBitstream/${CW_IDX}.txt" \
                    S
                
                ./bin/PostDecHammingDistance \
                    "${STAGE_II_DIR}/Decoding/RecoveryBitstream/${CW_IDX}.txt" \
                    "$ORIGINAL_BITSTREAM" \
                    "${STAGE_II_ERR}" \
                    "$CW_IDX"
            done  
            awk -v total1="$failed_count" "\$2==0 {count1++} END{printf \"Recovered: %d/%d\nRecoveryRate: %.4f\n\", count1, total1, count1/total1}" "${STAGE_II_ERR}"
            
            # Check whether all failed codewords were recovered
            STAGE_II_RECOVERED=$(awk '$2==0 {count++} END {print count+0}' "$STAGE_II_ERR")

            if [ "$STAGE_II_RECOVERED" -eq "$failed_count" ]; then

                echo "Stage II Recovery Success!"

                FINAL_BITSTREAM_DIR="${STAGE_I_DIR}/Decoding/FinalRecoveryBitstream"
                mkdir -p "$FINAL_BITSTREAM_DIR"

                # Copy all Stage I bitstreams
                cp "${STAGE_I_DIR}/Decoding/RecoveryBitstream/"*.txt \
                "$FINAL_BITSTREAM_DIR/"

                # Replace failed codewords using Stage II results
                for CW_IDX in "${FAILED_CW[@]}"; do
                    cp "${STAGE_II_DIR}/Decoding/RecoveryBitstream/${CW_IDX}.txt" \
                    "${FINAL_BITSTREAM_DIR}/${CW_IDX}.txt"
                done

                # Reconstruct the original image
                ./bin/ImageParsing \
                    "$FINAL_BITSTREAM_DIR" \
                    "${STAGE_II_DIR}/Decoding/recovered_image.jpg"

                xdg-open "${STAGE_II_DIR}/Decoding/recovered_image.jpg" >/dev/null 2>&1 &

            else
                echo "Stage II Recovery Failed!"
            fi
        fi  
    done  
done  

echo -e "\nEnd Time: $(date)"
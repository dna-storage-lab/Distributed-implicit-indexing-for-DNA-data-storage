#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <cuda_runtime.h>
#include <omp.h>
#include <time.h>

#define NUM_WATERMARKS 119922
#define WM_LENGTH 160
#define MAX_READ_LEN 168

#define BATCH_SIZE 6400
#define CORR_BIAS 100000u

typedef struct {
    char seq[MAX_READ_LEN + 1];
    int len;
    unsigned long long best_key;
} ReadResult;

int (*watermarks)[2][WM_LENGTH] = NULL;
int SHIFT_RANGE = 1;    
int NUM_GPUS = 2;       

#define CUDA_CHECK(call) \
do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

void inv(char* x) {
    int n = (int)strlen(x);
    for (int i = 0; i <= (n - 1) / 2; i++) {
        int j = n - 1 - i;
        char temp = x[i];
        x[i] = x[j];
        x[j] = temp;
    }
}

void complement(char *input, char *output) {
    int len = (int)strlen(input);
    for (int i = 0; i < len; i++) {
        if (input[i] == 'A') output[i] = 'T';
        else if (input[i] == 'T') output[i] = 'A';
        else if (input[i] == 'G') output[i] = 'C';
        else if (input[i] == 'C') output[i] = 'G';
        else output[i] = 'N';
    }
    output[len] = '\0';
}

void load_watermarks() {
    watermarks = (int(*)[2][WM_LENGTH])malloc((size_t)NUM_WATERMARKS * sizeof(int[2][WM_LENGTH]));
    FILE *fp = fopen("configure/watermark.txt", "r");
    if (!fp) { perror("[ERROR] Watermark file"); exit(1); }

    char line_buffer[512];
    int wm_idx = 0;
    while (fgets(line_buffer, sizeof(line_buffer), fp) && wm_idx < NUM_WATERMARKS) {
        int count = 0;
        for (int i = 0; line_buffer[i] != '\0' && count < 320; i++) {
            if (line_buffer[i] == '0' || line_buffer[i] == '1') {
                if (count < 160) watermarks[wm_idx][0][count] = (line_buffer[i] == '1' ? 1 : -1);
                else watermarks[wm_idx][1][count - 160] = (line_buffer[i] == '1' ? 1 : -1);
                count++;
            }
        }
        wm_idx++;
    }
    fclose(fp);
}

__device__ __forceinline__ unsigned long long pack_key(int corr, unsigned int linear) {
    unsigned int corr_u = (unsigned int)(corr + (int)CORR_BIAS);
    unsigned int tb     = 0xFFFFFFFFu - linear; 
    return (((unsigned long long)corr_u) << 32) | (unsigned long long)tb;
}

__global__ void compute_and_update_best_kernel(
    const int * __restrict__ d_watermarks,
    const int * __restrict__ d_sequences,
    const int * __restrict__ d_read_lengths,
    unsigned long long * __restrict__ d_best_keys,
    int batch_size, int max_read_len,
    int shift_range  
){
    int wm_id    = blockIdx.x * blockDim.x + threadIdx.x;
    int batch_id = blockIdx.y * blockDim.y + threadIdx.y;
    int num_shifts = shift_range * 2 + 1; 
    int shift_idx = (int)(blockIdx.z % num_shifts);
    int dir       = (int)(blockIdx.z / num_shifts);

    if (batch_id >= batch_size || wm_id >= NUM_WATERMARKS) return;

    int s = shift_idx - shift_range;  
    int L = d_read_lengths[batch_id];
    
    size_t seq_offset = (size_t)batch_id * 4u * max_read_len + (size_t)dir * 2u * max_read_len;
    const int *r0 = &d_sequences[seq_offset];
    const int *r1 = &d_sequences[seq_offset + max_read_len];
    const int *wm0 = &d_watermarks[(size_t)wm_id * 2u * WM_LENGTH];
    const int *wm1 = wm0 + WM_LENGTH;

    int corr = 0;
    if (s >= 0) {
        int overlap = min(L - s, WM_LENGTH);
        for (int j = 0; j < overlap; j++) corr += r0[s + j] * wm0[j] + r1[s + j] * wm1[j];
    } else {
        int drop = -s;
        int overlap = min(L, WM_LENGTH - drop);
        for (int j = 0; j < overlap; j++) corr += r0[j] * wm0[drop + j] + r1[j] * wm1[drop + j];
    }

    unsigned int linear = ((unsigned int)wm_id * (unsigned int)num_shifts + (unsigned int)shift_idx) * 2u + (unsigned int)dir;
    atomicMax((unsigned long long*)&d_best_keys[batch_id], pack_key(corr, linear));
}

static inline void unpack_key(unsigned long long key, int *corr, unsigned int *linear) {
    *corr = (int)(key >> 32) - (int)CORR_BIAS;
    *linear = 0xFFFFFFFFu - (unsigned int)(key & 0xFFFFFFFFu);
}

static inline int calc_overlap_by_new_shift(int L, int shift_out, int shift_range) {
    if (shift_out >= 0) {
        int ov = WM_LENGTH - shift_out;
        if (ov < 0) ov = 0;
        return (L < ov) ? L : ov;
    } else {
        int ov = L + shift_out; 
        if (ov < 0) ov = 0;
        return (ov < WM_LENGTH) ? ov : WM_LENGTH;
    }
}

static long long count_fastq_reads(FILE *fp) {
    char line[1024];
    long long total_reads = 0;

    while (1) {
        if (!fgets(line, sizeof(line), fp)) break;

        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr,
                    "[ERROR] Incomplete FASTQ record after read %lld: missing sequence line.\n",
                    total_reads);
            exit(EXIT_FAILURE);
        }

        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr,
                    "[ERROR] Incomplete FASTQ record after read %lld: missing '+' line.\n",
                    total_reads);
            exit(EXIT_FAILURE);
        }

        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr,
                    "[ERROR] Incomplete FASTQ record after read %lld: missing quality line.\n",
                    total_reads);
            exit(EXIT_FAILURE);
        }

        total_reads++;
    }

    rewind(fp);
    return total_reads;
}

void process_reads_gpu_multi(double coverage, int shift_range, int num_gpus,
                             const char *input_path, const char *output_path) {
    char seqname[200], outname[200];
    sprintf(seqname, input_path, coverage);
    sprintf(outname, output_path, coverage);

    FILE *fseq = fopen(seqname, "r");
    FILE *fout = fopen(outname, "w");
    if (!fseq || !fout) { printf("[ERROR] File access denied\n"); exit(1); }

    long long total_reads = count_fastq_reads(fseq);
    if (total_reads <= 0) {
        fprintf(stderr, "[ERROR] No complete FASTQ reads found in %s\n", seqname);
        fclose(fseq);
        fclose(fout);
        exit(EXIT_FAILURE);
    }

    int **d_watermarks = (int**)malloc(num_gpus * sizeof(int*));
    int **d_sequences = (int**)malloc(num_gpus * sizeof(int*));
    int **d_read_lengths = (int**)malloc(num_gpus * sizeof(int*));
    unsigned long long **d_best_keys = (unsigned long long**)malloc(num_gpus * sizeof(unsigned long long*));

    for (int g = 0; g < num_gpus; g++) {
        cudaSetDevice(g);
        size_t wm_size = (size_t)NUM_WATERMARKS * 2 * WM_LENGTH * sizeof(int);
        CUDA_CHECK(cudaMalloc(&d_watermarks[g], wm_size));
        CUDA_CHECK(cudaMalloc(&d_sequences[g], (size_t)BATCH_SIZE * 4 * MAX_READ_LEN * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_read_lengths[g], BATCH_SIZE * sizeof(int)));
        CUDA_CHECK(cudaMalloc(&d_best_keys[g], BATCH_SIZE * sizeof(unsigned long long)));
        CUDA_CHECK(cudaMemcpy(d_watermarks[g], watermarks, wm_size, cudaMemcpyHostToDevice));
    }

    long long total_processed = 0;
    int batch_id_counter = 0;
    long long reads_per_round = (long long)num_gpus * BATCH_SIZE;
    long long total_batches = (total_reads + reads_per_round - 1) / reads_per_round;
    int num_shifts = shift_range * 2 + 1;  
    
    clock_t start_time = clock();

    while (total_processed < total_reads) {
        ReadResult *batch_results = (ReadResult*)calloc(num_gpus * BATCH_SIZE, sizeof(ReadResult));
        int total_in_this_round = 0;

        for (int i = 0; i < num_gpus * BATCH_SIZE && total_processed + i < total_reads; i++) {
            char line[1024];
            if (!fgets(line, 1024, fseq)) break; 
            if (!fgets(line, 1024, fseq)) break; 
            line[strcspn(line, "\n")] = 0;
            strncpy(batch_results[i].seq, line, MAX_READ_LEN);
            batch_results[i].seq[MAX_READ_LEN] = '\0';
            batch_results[i].len = (int)strlen(batch_results[i].seq);
            if (!fgets(line, 1024, fseq)) break; 
            if (!fgets(line, 1024, fseq)) break; 
            total_in_this_round++;
        }

        if (total_in_this_round == 0) { free(batch_results); break; }

        batch_id_counter++;
        printf("\r[PROGRESS] Batch %d/%lld | Total Reads: %lld/%lld",
                batch_id_counter,
                total_batches,
                total_processed + total_in_this_round,
                total_reads);
        fflush(stdout);

        #pragma omp parallel num_threads(num_gpus)
        {
            int tid = omp_get_thread_num();
            cudaSetDevice(tid);

            int start_idx = tid * BATCH_SIZE;
            int my_batch_size = (total_in_this_round > start_idx) ? (total_in_this_round - start_idx) : 0;
            if (my_batch_size > BATCH_SIZE) my_batch_size = BATCH_SIZE;

            if (my_batch_size > 0) {
                int *h_sequences = (int*)calloc((size_t)my_batch_size * 4 * MAX_READ_LEN, sizeof(int));
                int *h_read_lengths = (int*)malloc(my_batch_size * sizeof(int));

                for (int i = 0; i < my_batch_size; i++) {
                    ReadResult *res = &batch_results[start_idx + i];
                    h_read_lengths[i] = res->len;

                    for (int j = 0; j < res->len; j++) {
                        int off = i * 4 * MAX_READ_LEN + j;
                        char c = res->seq[j];
                        h_sequences[off] = (c == 'G' || c == 'C') ? 1 : -1;
                        h_sequences[off + MAX_READ_LEN] = (c == 'T' || c == 'C') ? 1 : -1;
                    }

                    char rseq[MAX_READ_LEN + 1];
                    strcpy(rseq, res->seq);
                    inv(rseq);
                    complement(rseq, rseq);

                    for (int j = 0; j < res->len; j++) {
                        int off = i * 4 * MAX_READ_LEN + 2 * MAX_READ_LEN + j;
                        char c = rseq[j];
                        h_sequences[off] = (c == 'G' || c == 'C') ? 1 : -1;
                        h_sequences[off + MAX_READ_LEN] = (c == 'T' || c == 'C') ? 1 : -1;
                    }
                }

                CUDA_CHECK(cudaMemcpy(d_sequences[tid], h_sequences,
                    (size_t)my_batch_size * 4 * MAX_READ_LEN * sizeof(int),
                    cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemcpy(d_read_lengths[tid], h_read_lengths,
                    my_batch_size * sizeof(int),
                    cudaMemcpyHostToDevice));
                CUDA_CHECK(cudaMemset(d_best_keys[tid], 0, my_batch_size * sizeof(unsigned long long)));

                dim3 blockDim(16, 16);
                dim3 gridDim((NUM_WATERMARKS + blockDim.x - 1) / blockDim.x,
                             (my_batch_size + blockDim.y - 1) / blockDim.y,
                             num_shifts * 2); 

                compute_and_update_best_kernel<<<gridDim, blockDim>>>(
                    d_watermarks[tid],
                    d_sequences[tid],
                    d_read_lengths[tid],
                    d_best_keys[tid],
                    my_batch_size,
                    MAX_READ_LEN,
                    shift_range  
                );
                CUDA_CHECK(cudaGetLastError());
                CUDA_CHECK(cudaDeviceSynchronize());

                unsigned long long *h_best_keys = (unsigned long long*)malloc(my_batch_size * sizeof(unsigned long long));
                CUDA_CHECK(cudaMemcpy(h_best_keys, d_best_keys[tid],
                    my_batch_size * sizeof(unsigned long long),
                    cudaMemcpyDeviceToHost));

                for (int i = 0; i < my_batch_size; i++) batch_results[start_idx + i].best_key = h_best_keys[i];

                free(h_sequences);
                free(h_read_lengths);
                free(h_best_keys);
            }
        }


        for (int i = 0; i < total_in_this_round; i++) {
            int corr;
            unsigned int linear;
            unpack_key(batch_results[i].best_key, &corr, &linear);

            int dir = linear % 2;
            int s_idx = (linear / 2) % num_shifts;
            int wm_id = (linear / 2) / num_shifts;
            int s_internal = s_idx - shift_range; 
            int shift_out  = -s_internal;   
            int L = batch_results[i].len;
            int overlap_out = calc_overlap_by_new_shift(L, shift_out, shift_range);

            fprintf(fout, "%s %d %d %d %d %d %d\n",
                batch_results[i].seq,
                corr,
                shift_out,
                dir,
                overlap_out,
                wm_id + 1,
                L);
        }

        total_processed += total_in_this_round;
        free(batch_results);
    }

    clock_t end_time = clock();
    for (int g = 0; g < num_gpus; g++) {
        cudaSetDevice(g);
        cudaFree(d_watermarks[g]); 
        cudaFree(d_sequences[g]); 
        cudaFree(d_read_lengths[g]); 
        cudaFree(d_best_keys[g]);
    }
  
    free(d_watermarks);
    free(d_sequences);
    free(d_read_lengths);
    free(d_best_keys);
    
    fclose(fseq);
    fclose(fout);
}


int main(int argc, char *argv[]) {
    if (argc != 6) { 
        printf("Usage: %s <coverage> <shift_range> <num_gpus> <input_path> <output_path>\n", argv[0]); 
        printf("Example: %s 1.0 1 2 \"./result/%%.1fx/1/StepI_result/PrimerIdentification/read_without_primer.fastq\" \"./result/%%.1fx/1/StepI_result/Correlation/corr_result.txt\"\n", argv[0]);
        return 1; 
    }

    double coverage = atof(argv[1]);
    SHIFT_RANGE = atoi(argv[2]);
    NUM_GPUS = atoi(argv[3]);
    const char *input_path = argv[4];
    const char *output_path = argv[5];
    
    if (SHIFT_RANGE < 0) {
        printf("[ERROR] shift_range must be >= 0\n");
        return 1;
    }
    
    if (NUM_GPUS <= 0) {
        printf("[ERROR] num_gpus must be > 0\n");
        return 1;
    }
    
    int dev_count;
    cudaGetDeviceCount(&dev_count);
    printf("------------------------------------------\n");
    printf("Detected %d CUDA devices.\n", dev_count);
    for(int i=0; i<dev_count; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
    }
    
    if (NUM_GPUS > dev_count) {
        printf("[WARNING] Requested %d GPUs but only %d available. Using %d GPUs.\n", 
               NUM_GPUS, dev_count, dev_count);
        NUM_GPUS = dev_count;
    }
    printf("------------------------------------------\n");

    load_watermarks();
    process_reads_gpu_multi(coverage, SHIFT_RANGE, NUM_GPUS, input_path, output_path);
    free(watermarks);
    return 0;
}
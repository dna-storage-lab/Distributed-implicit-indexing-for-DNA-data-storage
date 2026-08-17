#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#ifndef MADV_SEQUENTIAL
#include <sys/mman.h>
#endif

#define MAX_GROUPS 8096
#define FLOATS_PER_GROUP 16
#define NUM_CW 474

typedef struct {
    long double *lastLLR1;
    long double *lastLLR2;
} CodewordData;


typedef struct {
    const char *p_prob;
    int num_groups;
    int truePos;
} RecIndex;

typedef struct {
    RecIndex *arr;
    size_t n;
    size_t cap;
} RecList;

static inline void reclist_push(RecList *rl, RecIndex r) {
    if (rl->n == rl->cap) {
        size_t newcap = rl->cap ? rl->cap * 2 : 1024;
        RecIndex *na = (RecIndex*)realloc(rl->arr, newcap * sizeof(RecIndex));
        if (!na) { perror("realloc RecIndex"); exit(1); }
        rl->arr = na;
        rl->cap = newcap;
    }
    rl->arr[rl->n++] = r;
}


static inline void skip_ws(const char **p, const char *end) {
    const char *q = *p;
    while (q < end) {
        unsigned char c = (unsigned char)*q;
        if (c > ' ') break;
        q++;
    }
    *p = q;
}

static inline long int fast_atol_constp(const char **p, const char *end) {
    skip_ws(p, end);
    char *ep = NULL;
    long int v = strtol(*p, &ep, 10);
    *p = ep;
    return v;
}

static inline long double fast_atold_constp(const char **p, const char *end) {
    skip_ws(p, end);
    char *ep = NULL;
    long double v = strtold(*p, &ep);
    *p = ep;
    return v;
}

static inline const char* next_line(const char *p, const char *end) {
    const void *nl = memchr(p, '\n', (size_t)(end - p));
    if (!nl) return end;
    return (const char*)nl + 1;
}


static int* parse_failed_rows(const char *failed_rows_str, int *failed_rows_count) {
    *failed_rows_count = 0;
    const char *tmp = failed_rows_str;
    while (*tmp != '\0') {
        if (*tmp == ' ') (*failed_rows_count)++;
        tmp++;
    }
    (*failed_rows_count)++;

    int *failed_rows = (int*)malloc((size_t)(*failed_rows_count) * sizeof(int));
    if (!failed_rows) {
        perror("Failed to allocate memory for failed_rows");
        return NULL;
    }

    char *str = strdup(failed_rows_str);
    if (!str) { perror("strdup"); free(failed_rows); return NULL; }

    char *saveptr = NULL;
    char *token = strtok_r(str, " ", &saveptr);
    int idx = 0;
    while (token != NULL) {
        int v = atoi(token) - 1; 
        if (v >= 0 && v < NUM_CW) failed_rows[idx++] = v;
        token = strtok_r(NULL, " ", &saveptr);
    }
    free(str);

    *failed_rows_count = idx;
    return failed_rows;
}

int main(int argc, char const *argv[]) {
    if (argc < 4) {
        printf("Usage: %s <LLR_file> <output_dir> [num_threads] <failed_rows>\n", argv[0]);
        return -1;
    }

    const char *inputfile = argv[1];
    const char *outputdir = argv[2];

    int threads = 64;
    const char *failed_rows_str = NULL;


    if (argc == 4) {
        failed_rows_str = argv[3];
    } else {
        threads = atoi(argv[3]);
        if (threads <= 0) threads = 64;
        failed_rows_str = argv[argc - 1];
    }


    int failed_rows_count = 0;
    int *failed_rows = parse_failed_rows(failed_rows_str, &failed_rows_count);
    if (!failed_rows || failed_rows_count == 0) {
        fprintf(stderr, "Failed to parse failed_rows or empty list.\n");
        return 1;
    }

    unsigned char is_failed[NUM_CW];
    memset(is_failed, 0, sizeof(is_failed));
    for (int i = 0; i < failed_rows_count; i++) is_failed[failed_rows[i]] = 1;

    CodewordData *cw_data = (CodewordData*)calloc(NUM_CW, sizeof(CodewordData));
    if (!cw_data) { perror("calloc cw_data"); return 1; }

    for (int cw = 0; cw < NUM_CW; cw++) {
        if (!is_failed[cw]) continue;
        cw_data[cw].lastLLR1 = (long double*)malloc((size_t)MAX_GROUPS * FLOATS_PER_GROUP * sizeof(long double));
        cw_data[cw].lastLLR2 = (long double*)malloc((size_t)MAX_GROUPS * FLOATS_PER_GROUP * sizeof(long double));
        if (!cw_data[cw].lastLLR1 || !cw_data[cw].lastLLR2) {
            fprintf(stderr, "Memory allocation failed at codeword %d\n", cw);
            return 1;
        }
    }

    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(static)
    for (int cw = 0; cw < NUM_CW; cw++) {
        if (!is_failed[cw]) continue;
        long double *a = cw_data[cw].lastLLR1;
        long double *b = cw_data[cw].lastLLR2;
        size_t n = (size_t)MAX_GROUPS * FLOATS_PER_GROUP;
        for (size_t k = 0; k < n; k++) {
            a[k] = 1.0L;
            b[k] = 1.0L;
        }
    }


    int fd = open(inputfile, O_RDONLY);
    if (fd == -1) { perror("open"); return 1; }

    struct stat sb;
    if (fstat(fd, &sb) != 0) { perror("fstat"); return 1; }

    const char *file_ptr = (const char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (file_ptr == MAP_FAILED) { perror("mmap"); return 1; }

    (void)madvise((void*)file_ptr, (size_t)sb.st_size, MADV_SEQUENTIAL);

    const char *cur = file_ptr;
    const char *end = file_ptr + sb.st_size;


    RecList buckets[NUM_CW];
    memset(buckets, 0, sizeof(buckets));

    long long rec_total = 0;
    long long rec_kept  = 0;

    while (cur < end) {
        const char *line1 = cur;
        if (line1 >= end) break;

        const char *line2 = next_line(line1, end);
        if (line2 >= end) break;

        const char *line3 = next_line(line2, end);
        if (line3 >= end) break;

        const char *next = next_line(line3, end);

        const char *p = line1;
        int position   = (int)fast_atol_constp(&p, end);
        int num_groups = (int)fast_atol_constp(&p, end);

        int cw_num  = (position - 1) / 253;
        int truePos = ((position - 1) % 253 * 160) / 5;

        if ((unsigned)cw_num < (unsigned)NUM_CW && is_failed[cw_num]) {
            RecIndex ri;
            ri.p_prob = line2;      
            ri.num_groups = num_groups;
            ri.truePos = truePos;
            reclist_push(&buckets[cw_num], ri);
            rec_kept++;
        }

        rec_total++;
        cur = next;

        if ((rec_total % 200000) == 0) {
            double pct = 100.0 * (double)(cur - file_ptr) / (double)sb.st_size;
        }

    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (int cw = 0; cw < NUM_CW; cw++) {
        if (!is_failed[cw]) continue;

        RecList *rl = &buckets[cw];
        long double *L1 = cw_data[cw].lastLLR1;
        long double *L2 = cw_data[cw].lastLLR2;

        for (size_t r = 0; r < rl->n; r++) {
            const RecIndex *ri = &rl->arr[r];
            int ng = ri->num_groups;
            int base = ri->truePos;
            if (ng <= 0) continue; 

            const char *pp = ri->p_prob;

            for (int i = 0; i < ng; i++) {
                int row = base + i;
                if ((unsigned)row < (unsigned)MAX_GROUPS) {
                    size_t off = (size_t)row * FLOATS_PER_GROUP;
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) {
                        L1[off + j] *= fast_atold_constp(&pp, end);
                    }
                } else {
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) (void)fast_atold_constp(&pp, end);
                }
            }

        
            for (int i = 0; i < ng; i++) {
                int row = base + i;
                if ((unsigned)row < (unsigned)MAX_GROUPS) {
                    size_t off = (size_t)row * FLOATS_PER_GROUP;
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) {
                        L2[off + j] *= fast_atold_constp(&pp, end);
                    }
                } else {
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) (void)fast_atold_constp(&pp, end);
                }
            }
        }
    }

    #pragma omp parallel
    {
        long double (*bit1)[2] = (long double (*)[2])malloc((size_t)MAX_GROUPS * 4 * sizeof(long double[2]));
        long double (*bit2)[2] = (long double (*)[2])malloc((size_t)MAX_GROUPS * 4 * sizeof(long double[2]));
        if (!bit1 || !bit2) {
            fprintf(stderr, "bit buffer malloc failed in thread %d\n", omp_get_thread_num());
            exit(1);
        }

        #pragma omp for schedule(dynamic, 1)
        for (int cw = 0; cw < NUM_CW; cw++) {
            if (!is_failed[cw]) continue;

            long double *L1 = cw_data[cw].lastLLR1;
            long double *L2 = cw_data[cw].lastLLR2;

            for (int i = 0; i < MAX_GROUPS; i++) {
                long double total1 = 0.0L, total2 = 0.0L;
                size_t off = (size_t)i * FLOATS_PER_GROUP;
                for (int j = 0; j < FLOATS_PER_GROUP; j++) {
                    total1 += L1[off + j];
                    total2 += L2[off + j];
                }
                if (total1 > 0.0L) {
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) L1[off + j] /= total1;
                }
                if (total2 > 0.0L) {
                    for (int j = 0; j < FLOATS_PER_GROUP; j++) L2[off + j] /= total2;
                }
            }


            for (int i = 0; i < MAX_GROUPS; i++) {
                size_t off = (size_t)i * 16;
                for (int ii = 0; ii < 4; ii++) {
                    long double zeroP1 = 0.0L, oneP1 = 0.0L;
                    long double zeroP2 = 0.0L, oneP2 = 0.0L;
                    for (int j = 0; j < 16; j++) {
                        if ((j >> ii) & 1) {
                            oneP1 += L1[off + j];
                            oneP2 += L2[off + j];
                        } else {
                            zeroP1 += L1[off + j];
                            zeroP2 += L2[off + j];
                        }
                    }

                    int out_idx = 4 * i + (3 - ii);
                    long double sum1 = zeroP1 + oneP1;
                    bit1[out_idx][0] = (sum1 == 0) ? 0.0L : zeroP1 / sum1;
                    bit1[out_idx][1] = (sum1 == 0) ? 0.0L : oneP1 / sum1;

                    long double sum2 = zeroP2 + oneP2;
                    bit2[out_idx][0] = (sum2 == 0) ? 0.0L : zeroP2 / sum2;
                    bit2[out_idx][1] = (sum2 == 0) ? 0.0L : oneP2 / sum2;
                }
            }


            char filename[300];
            snprintf(filename, sizeof(filename), "%s/%d.txt", outputdir, cw + 1);

            FILE *fw = fopen(filename, "w");
            if (!fw) {
                fprintf(stderr, "Failed to open output file %s: %s\n", filename, strerror(errno));
                continue;
            }

   
            static const size_t BUF_SZ = 8 * 1024 * 1024;
            char *buf = (char*)malloc(BUF_SZ);
            if (buf) setvbuf(fw, buf, _IOFBF, BUF_SZ);

            for (int block = 0; block < MAX_GROUPS * 4; block += 128) {
                for (int j = block; j < block + 128 && j < MAX_GROUPS * 4; j++) {
                    for (int b = 0; b < 2; b++) {
                        if (bit1[j][b] >= 0.99999) bit1[j][b] = 0.9999;
                        else if (bit1[j][b] <= 0.00001) bit1[j][b] = 0.0001;
                    }
                    fprintf(fw, "%Lf %Lf\n", bit1[j][0], bit1[j][1]);
                }
                for (int j = block; j < block + 128 && j < MAX_GROUPS * 4; j++) {
                    for (int b = 0; b < 2; b++) {
                        if (bit2[j][b] >= 0.99999) bit2[j][b] = 0.9999;
                        else if (bit2[j][b] <= 0.00001) bit2[j][b] = 0.0001;
                    }
                    fprintf(fw, "%Lf %Lf\n", bit2[j][0], bit2[j][1]);
                }
            }

            fclose(fw);
        }

        free(bit1);
        free(bit2);
    }


    for (int cw = 0; cw < NUM_CW; cw++) free(buckets[cw].arr);

    munmap((void*)file_ptr, sb.st_size);
    close(fd);

    for (int cw = 0; cw < NUM_CW; cw++) {
        if (!is_failed[cw]) continue;
        free(cw_data[cw].lastLLR1);
        free(cw_data[cw].lastLLR2);
    }
    free(cw_data);
    free(failed_rows);

    return 0;
}

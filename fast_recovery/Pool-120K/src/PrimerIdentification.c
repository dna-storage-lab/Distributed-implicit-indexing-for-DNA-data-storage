#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include "../lib/edlib/edlib.h"

#define MAXLEN 100000
#define HIST_BINS 100 

static int g_min_edit = 3;
static int g_primer_len = 20;
static int g_oligo_len = 150;
static int g_max_recursion_depth = 5;
static int g_min_length = 155;
static int g_max_length = 165;

static void inv(char *x, int n)
{
    int i, j, m = (n - 1) / 2;
    char temp;
    for (i = 0; i <= m; i++)
    {
        j = n - 1 - i;
        temp = x[i];
        x[i] = x[j];
        x[j] = temp;
    }
}
static void complement(char *input, char *output, int len)
{
    for (int i = 0; i < len; i++)
    {
        char c = input[i];
        output[i] = (c == 'A') ? 'T' : (c == 'T') ? 'A'
                                   : (c == 'G')   ? 'C'
                                   : (c == 'C')   ? 'G'
                                                  : 'N';
    }
}
static void strip_newline(char *s)
{
    size_t n = strlen(s);
    while (n && (s[n - 1] == '\n' || s[n - 1] == '\r'))
        s[--n] = '\0';
}
static int ensure_dir(const char *path)
{
    if (!path || !*path)
        return -1;
    if (mkdir(path, 0755) == 0)
        return 0;
    if (errno == EEXIST)
        return 0;
    return -1;
}
static void concat_file_to(FILE *dest, const char *src_path)
{
    FILE *src = fopen(src_path, "rb");
    if (!src)
        return;
    char buf[1 << 16];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), src)) > 0)
        fwrite(buf, 1, n, dest);
    fclose(src);
}
static double elapsed_sec(struct timeval a, struct timeval b)
{
    return (b.tv_sec - a.tv_sec) + (b.tv_usec - a.tv_usec) / 1e6;
}


typedef struct
{
    char *seq;
    int seqlen;
    int read_num;
} ReadRecord;

static int g_next_read_index = 0;
static pthread_mutex_t g_queue_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct
{
    int thread_id;
    int pair_count;
    int primer_num;
    char **primers;
    ReadRecord *reads;
    int num_reads;
    char *out_prefix;

    FILE **f_output_seq; 
    FILE **f_output_len; 

    int **local_pair_match_hist; 
    int *primer_edit_records;  

    int *local_allpairs_hist;               
    long long local_allpairs_total_matches; 
} WorkerArg;


static void recursive_primer_matching_forward(
    char **primers, char *line, int seqlen, int read_num, int thread_id,
    int *match_count, FILE *f_output, int *primer_edit_record, FILE *f_length_record, int depth)
{
    if (depth > g_max_recursion_depth)
    {
        
        return;
    }
    if (seqlen <= g_min_length || seqlen <= 0)
        return;

    char primerL[g_primer_len + 1];
    char primerR[g_primer_len + 1];
    memcpy(primerL, primers[0], g_primer_len);
    memcpy(primerR, primers[1], g_primer_len);
    primerL[g_primer_len] = '\0';
    primerR[g_primer_len] = '\0';

    EdlibAlignResult result1 = edlibAlign(primerL, g_primer_len, line, seqlen,
                                          edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    EdlibAlignResult result2 = edlibAlign(primerR, g_primer_len, line, seqlen,
                                          edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result1.numLocations <= 0 || result2.numLocations <= 0 ||
        result1.editDistance == -1 || result2.editDistance == -1)
    {
        edlibFreeAlignResult(result1);
        edlibFreeAlignResult(result2);
        return;
    }

    int start = *result1.startLocations;
    int end_forword = *result1.endLocations;
    int end = *result2.endLocations;
    int start_backword = *result2.startLocations;

    if (result1.editDistance <= g_min_edit && result2.editDistance <= g_min_edit)
    {
        if (start_backword > end_forword && start_backword - end_forword - 1 < g_oligo_len)
        {
            int sequence_length = start_backword - end_forword - 1;
            if (sequence_length > 0 && sequence_length <= g_max_length && sequence_length >= g_min_length)
            {
                (*match_count)++;

                fprintf(f_output, "@T%d_Read_%d_Match_%d_forward\n", thread_id, read_num, *match_count);
                for (int j = end_forword + 1; j < start_backword && j < seqlen; j++)
                    fputc(line[j], f_output);
                fputc('\n', f_output);

                fprintf(f_output, "+\n");

                for (int j = end_forword + 1; j < start_backword && j < seqlen; j++)
                    fputc('I', f_output);
                fputc('\n', f_output);

                fflush(f_output);
                (*primer_edit_record)++;
                fprintf(f_length_record, "%d\n", sequence_length);
                fflush(f_length_record);
            }
        }

        if ((start_backword - end_forword - 1) >= g_oligo_len)
        {
            int middle_len = start_backword - end_forword - 1;
            if (middle_len > 0 && middle_len < MAXLEN && middle_len >= g_min_length + 20)
            {
                char *middle_subseq = (char *)malloc(middle_len + 1);
                if (middle_subseq)
                {
                    strncpy(middle_subseq, line + end_forword + 1, middle_len);
                    middle_subseq[middle_len] = '\0';
                    recursive_primer_matching_forward(primers, middle_subseq, middle_len, read_num, thread_id,
                                                      match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                    free(middle_subseq);
                }
            }
        }

        if (start_backword > end_forword)
        {
            if (start >= g_min_length + 20 && start < seqlen)
            {
                char *left_subseq = (char *)malloc(start + 1);
                if (left_subseq)
                {
                    strncpy(left_subseq, line, start);
                    left_subseq[start] = '\0';
                    recursive_primer_matching_forward(primers, left_subseq, start, read_num, thread_id,
                                                      match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                    free(left_subseq);
                }
            }
            if (end >= 0 && end < seqlen - 1)
            {
                int right_len = seqlen - end - 1;
                if (right_len >= g_min_length + 20)
                {
                    char *right_subseq = (char *)malloc(right_len + 1);
                    if (right_subseq)
                    {
                        strncpy(right_subseq, line + end + 1, right_len);
                        right_subseq[right_len] = '\0';
                        recursive_primer_matching_forward(primers, right_subseq, right_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(right_subseq);
                    }
                }
            }
        }
        else if (end < start || (start_backword < start && start < end && end < end_forword))
        {
            if (start_backword > 0 && end >= 0 && end < seqlen)
            {
                int left_len = end + 1;
                if (left_len >= g_min_length + 20)
                {
                    char *left_subseq = (char *)malloc(left_len + 1);
                    if (left_subseq)
                    {
                        strncpy(left_subseq, line, left_len);
                        left_subseq[left_len] = '\0';
                        recursive_primer_matching_forward(primers, left_subseq, left_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(left_subseq);
                    }
                }
            }
            if (end_forword >= 0 && end_forword < seqlen - 1 && start < seqlen)
            {
                int right_len = seqlen - start;
                if (right_len >= g_min_length + 20 && start >= 0)
                {
                    char *right_subseq = (char *)malloc(right_len + 1);
                    if (right_subseq)
                    {
                        strncpy(right_subseq, line + start, right_len);
                        right_subseq[right_len] = '\0';
                        recursive_primer_matching_forward(primers, right_subseq, right_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(right_subseq);
                    }
                }
            }
        }
    }

    edlibFreeAlignResult(result1);
    edlibFreeAlignResult(result2);
}

static void recursive_primer_matching_reverse(
    char **primers, char *line, int seqlen, int read_num, int thread_id,
    int *match_count, FILE *f_output, int *primer_edit_record, FILE *f_length_record, int depth)
{
    if (depth > g_max_recursion_depth)
    {
       
        return;
    }
    if (seqlen < g_min_length || seqlen <= 0)
        return;

    char primerL[g_primer_len + 1];
    char primerR[g_primer_len + 1];
    char primerL_rv[g_primer_len + 1];
    char primerR_rv[g_primer_len + 1];
    char line_temp[MAXLEN] = {0};

    memcpy(primerL, primers[0], g_primer_len);
    memcpy(primerR, primers[1], g_primer_len);
    primerL[g_primer_len] = '\0';
    primerR[g_primer_len] = '\0';

    inv(primerL, g_primer_len);
    complement(primerL, primerL_rv, g_primer_len);
    inv(primerR, g_primer_len);
    complement(primerR, primerR_rv, g_primer_len);
    primerL_rv[g_primer_len] = '\0';
    primerR_rv[g_primer_len] = '\0';

    EdlibAlignResult result4 = edlibAlign(primerL_rv, g_primer_len, line, seqlen,
                                          edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    EdlibAlignResult result3 = edlibAlign(primerR_rv, g_primer_len, line, seqlen,
                                          edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    if (result3.numLocations <= 0 || result4.numLocations <= 0 ||
        result3.editDistance == -1 || result4.editDistance == -1)
    {
        edlibFreeAlignResult(result3);
        edlibFreeAlignResult(result4);
        return;
    }

    int start = *result3.startLocations;
    int end_forword = *result3.endLocations;
    int start_backword = *result4.startLocations;
    int end = *result4.endLocations;

    if (result3.editDistance <= g_min_edit && result4.editDistance <= g_min_edit)
    {
        if (start_backword > end_forword && start_backword - end_forword - 1 < g_oligo_len)
        {
            int sequence_length = start_backword - end_forword - 1;
            int middle_len = sequence_length;
            if (middle_len > 0 && middle_len <= g_max_length && middle_len >= g_min_length)
            {
                for (int j = end_forword + 1; j < start_backword; j++)
                {
                    line_temp[j - end_forword - 1] = line[j];
                }
                inv(line_temp, middle_len);
                complement(line_temp, line_temp, middle_len);

                (*match_count)++;

                fprintf(f_output, "@T%d_Read_%d_Match_%d_reverse\n", thread_id, read_num, *match_count);
                for (int j = 0; j < middle_len; j++)
                    fputc(line_temp[j], f_output);
                fputc('\n', f_output);

                fprintf(f_output, "+\n");
                for (int j = 0; j < middle_len; j++)
                    fputc('I', f_output);
                fputc('\n', f_output);

                fflush(f_output);
                (*primer_edit_record)++;
                fprintf(f_length_record, "%d\n", sequence_length);
                fflush(f_length_record);
            }
        }

        if ((start_backword - end_forword - 1) >= g_oligo_len)
        {
            int middle_len = start_backword - end_forword - 1;
            if (middle_len > 0 && middle_len < MAXLEN && middle_len >= g_min_length + 20)
            {
                char *middle_subseq = (char *)malloc(middle_len + 1);
                if (middle_subseq)
                {
                    strncpy(middle_subseq, line + end_forword + 1, middle_len);
                    middle_subseq[middle_len] = '\0';
                    recursive_primer_matching_reverse(primers, middle_subseq, middle_len, read_num, thread_id,
                                                      match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                    free(middle_subseq);
                }
            }
        }

        if (start_backword > end_forword)
        {
            if (start >= g_min_length + 20 && start < seqlen)
            {
                char *left_subseq = (char *)malloc(start + 1);
                if (left_subseq)
                {
                    strncpy(left_subseq, line, start);
                    left_subseq[start] = '\0';
                    recursive_primer_matching_reverse(primers, left_subseq, start, read_num, thread_id,
                                                      match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                    free(left_subseq);
                }
            }
            if (end >= 0 && end < seqlen - 1)
            {
                int right_len = seqlen - end - 1;
                if (right_len >= g_min_length + 20)
                {
                    char *right_subseq = (char *)malloc(right_len + 1);
                    if (right_subseq)
                    {
                        strncpy(right_subseq, line + end + 1, right_len);
                        right_subseq[right_len] = '\0';
                        recursive_primer_matching_reverse(primers, right_subseq, right_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(right_subseq);
                    }
                }
            }
        }
        else if (end < start || (start_backword < start && start < end && end < end_forword))
        {
            if (start_backword > 0 && end >= 0 && end < seqlen)
            {
                int left_len = end + 1;
                if (left_len >= g_min_length + 20)
                {
                    char *left_subseq = (char *)malloc(left_len + 1);
                    if (left_subseq)
                    {
                        strncpy(left_subseq, line, left_len);
                        left_subseq[left_len] = '\0';
                        recursive_primer_matching_reverse(primers, left_subseq, left_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(left_subseq);
                    }
                }
            }
            if (end_forword >= 0 && end_forword < seqlen - 1 && start < seqlen)
            {
                int right_len = seqlen - start;
                if (right_len >= g_min_length + 20 && start >= 0)
                {
                    char *right_subseq = (char *)malloc(right_len + 1);
                    if (right_subseq)
                    {
                        strncpy(right_subseq, line + start, right_len);
                        right_subseq[right_len] = '\0';
                        recursive_primer_matching_reverse(primers, right_subseq, right_len, read_num, thread_id,
                                                          match_count, f_output, primer_edit_record, f_length_record, depth + 1);
                        free(right_subseq);
                    }
                }
            }
        }
    }

    edlibFreeAlignResult(result3);
    edlibFreeAlignResult(result4);
}

static void *worker_run(void *arg)
{
    WorkerArg *w = (WorkerArg *)arg;

    for (;;)
    {
        pthread_mutex_lock(&g_queue_mutex);
        int idx = g_next_read_index++;
        pthread_mutex_unlock(&g_queue_mutex);
        if (idx >= w->num_reads)
            break;

        ReadRecord *r = &w->reads[idx];

        int matches_this_read_total = 0; 

        for (int p = 0; p < w->pair_count; ++p)
        {
            char *primer_pair[2];
            primer_pair[0] = w->primers[p * 2];
            primer_pair[1] = w->primers[p * 2 + 1];

            int match_count_forward = 0;
            int match_count_reverse = 0;

            recursive_primer_matching_forward(
                primer_pair, r->seq, r->seqlen, r->read_num, w->thread_id,
                &match_count_forward, w->f_output_seq[p], &w->primer_edit_records[p], w->f_output_len[p], 0);

            recursive_primer_matching_reverse(
                primer_pair, r->seq, r->seqlen, r->read_num, w->thread_id,
                &match_count_reverse, w->f_output_seq[p], &w->primer_edit_records[p], w->f_output_len[p], 0);

            int total_match = match_count_forward + match_count_reverse;
            int bin = total_match < HIST_BINS ? total_match : HIST_BINS - 1;

            w->local_pair_match_hist[p][bin]++;

            matches_this_read_total += total_match;
        }


        int all_bin = matches_this_read_total < HIST_BINS ? matches_this_read_total : HIST_BINS - 1;
        w->local_allpairs_hist[all_bin]++;
        w->local_allpairs_total_matches += matches_this_read_total;
    }
    return NULL;
}


int main(int argc, char *argv[])
{
    if (argc != 12)
    {
        fprintf(stderr,
                "Usage: %s <fastq_file> <primer_list> <primer_num> "
                "<primer_len> <min_edit> <min_length> <max_length> "
                "<oligo_len> <max_recursion_depth> <result_output_prefix> <threads>\n",
                argv[0]);
        return 1;
    }

    struct timeval tv_start, tv_after_read, tv_proc_start, tv_after_proc, tv_end;
    gettimeofday(&tv_start, NULL);

    char *fastq_file = argv[1];
    char *primer_list = argv[2];
    int primer_num = atoi(argv[3]);
    g_primer_len = atoi(argv[4]);
    g_min_edit = atoi(argv[5]);
    g_min_length = atoi(argv[6]);
    g_max_length = atoi(argv[7]);
    g_oligo_len = atoi(argv[8]);
    g_max_recursion_depth = atoi(argv[9]);
    char *result_output_prefix = argv[10];
    int num_threads = atoi(argv[11]);
    if (num_threads <= 0)
        num_threads = 1;

    if (g_primer_len <= 0)
    {
        fprintf(stderr, "primer_len must be > 0\n");
        return 1;
    }


    char **primers = (char **)malloc(primer_num * sizeof(char *));
    for (int i = 0; i < primer_num; i++)
    {
        primers[i] = (char *)malloc((g_primer_len + 1) * sizeof(char));
        if (!primers[i])
        {
            fprintf(stderr, "Out of memory.\n");
            return 1;
        }
    }
    FILE *f_primer_list = fopen(primer_list, "r");
    if (!f_primer_list)
    {
        fprintf(stderr, "Error opening primer list file: %s\n", primer_list);
        return 1;
    }
    for (int i = 0; i < primer_num; i++)
    {
        char buf[4096];
        if (fscanf(f_primer_list, "%4095s", buf) != 1)
        {
            fprintf(stderr, "Error reading primer %d\n", i);
            fclose(f_primer_list);
            return 1;
        }
        int L = (int)strlen(buf);
        if (L != g_primer_len)
        {
            fprintf(stderr, "Primer %d length %d != primer_len %d\n", i, L, g_primer_len);
            fclose(f_primer_list);
            return 1;
        }
        memcpy(primers[i], buf, g_primer_len + 1);
    }
    fclose(f_primer_list);


    FILE *f_fastq = fopen(fastq_file, "r");
    if (!f_fastq)
    {
        fprintf(stderr, "Error opening FASTQ file: %s\n", fastq_file);
        exit(1);
    }

    char line1[MAXLEN] = {0}, line2[MAXLEN] = {0}, line3[MAXLEN] = {0}, line4[MAXLEN] = {0};
    size_t cap = 10000, num_reads = 0;
    ReadRecord *reads = (ReadRecord *)malloc(cap * sizeof(ReadRecord));
    if (!reads)
    {
        fprintf(stderr, "Out of memory.\n");
        return 1;
    }

    while (1)
    {
        if (!fgets(line1, MAXLEN, f_fastq))
            break;
        if (!fgets(line2, MAXLEN, f_fastq))
            break;
        if (!fgets(line3, MAXLEN, f_fastq))
            break;
        if (!fgets(line4, MAXLEN, f_fastq))
            break;

        strip_newline(line2);
        int seqlen = (int)strlen(line2);
        if (seqlen <= 0 || seqlen > MAXLEN - 1)
            continue;

        if (num_reads == cap)
        {
            cap *= 2;
            ReadRecord *tmp = (ReadRecord *)realloc(reads, cap * sizeof(ReadRecord));
            if (!tmp)
            {
                fprintf(stderr, "Out of memory during resize.\n");
                return 1;
            }
            reads = tmp;
        }
        reads[num_reads].seq = (char *)malloc(seqlen + 1);
        if (!reads[num_reads].seq)
        {
            fprintf(stderr, "Out of memory for sequence.\n");
            return 1;
        }
        memcpy(reads[num_reads].seq, line2, seqlen + 1);
        reads[num_reads].seqlen = seqlen;
        reads[num_reads].read_num = (int)num_reads + 1;
        num_reads++;

        if (num_reads % 100000 == 0)
        {
            fflush(stdout);
        }
    }
    fclose(f_fastq);
    gettimeofday(&tv_after_read, NULL);
    int pair_count = primer_num / 2;

    if (ensure_dir(result_output_prefix) != 0)
    {
        fprintf(stderr, "Failed to create/access output directory: %s\n", result_output_prefix);
        return 1;
    }


    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    WorkerArg **args = (WorkerArg **)malloc(num_threads * sizeof(WorkerArg *));
    if (!threads || !args)
    {
        fprintf(stderr, "Out of memory.\n");
        return 1;
    }

    for (int t = 0; t < num_threads; ++t)
    {
        args[t] = (WorkerArg *)calloc(1, sizeof(WorkerArg));
        WorkerArg *w = args[t];
        w->thread_id = t;
        w->pair_count = pair_count;
        w->primer_num = primer_num;
        w->primers = primers;
        w->reads = reads;
        w->num_reads = (int)num_reads;
        w->out_prefix = result_output_prefix;

        w->f_output_seq = (FILE **)malloc(pair_count * sizeof(FILE *));
        w->f_output_len = (FILE **)malloc(pair_count * sizeof(FILE *));
        w->primer_edit_records = (int *)calloc(pair_count, sizeof(int));
        w->local_pair_match_hist = (int **)malloc(pair_count * sizeof(int *));
        w->local_allpairs_hist = (int *)calloc(HIST_BINS, sizeof(int));
        w->local_allpairs_total_matches = 0;

        if (!w->f_output_seq || !w->f_output_len || !w->primer_edit_records ||
            !w->local_pair_match_hist || !w->local_allpairs_hist)
        {
            fprintf(stderr, "Out of memory (thread %d).\n", t);
            return 1;
        }

        for (int p = 0; p < pair_count; ++p)
        {
            char path_seq[1024], path_len[1024];
            snprintf(path_seq, sizeof(path_seq), "%s/pair%d_seq_t%d.fastq", result_output_prefix, p + 1, t);
            snprintf(path_len, sizeof(path_len), "%s/pair%d_length_t%d.txt", result_output_prefix, p + 1, t);

            w->f_output_seq[p] = fopen(path_seq, "wb");
            w->f_output_len[p] = fopen(path_len, "wb");
            if (!w->f_output_seq[p] || !w->f_output_len[p])
            {
                fprintf(stderr, "Error opening thread files for pair %d (thread %d)\n", p + 1, t);
                return 1;
            }

            w->local_pair_match_hist[p] = (int *)calloc(HIST_BINS, sizeof(int));
            if (!w->local_pair_match_hist[p])
            {
                fprintf(stderr, "Out of memory for local_pair_match_hist.\n");
                return 1;
            }
        }
    }

    g_next_read_index = 0;
    gettimeofday(&tv_proc_start, NULL);
    for (int t = 0; t < num_threads; ++t)
    {
        if (pthread_create(&threads[t], NULL, worker_run, args[t]) != 0)
        {
            fprintf(stderr, "Failed to create thread %d\n", t);
            return 1;
        }
    }
    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], NULL);
    gettimeofday(&tv_after_proc, NULL);

    for (int t = 0; t < num_threads; ++t)
    {
        WorkerArg *w = args[t];
        for (int p = 0; p < pair_count; ++p)
        {
            if (w->f_output_seq[p])
                fclose(w->f_output_seq[p]);
            if (w->f_output_len[p])
                fclose(w->f_output_len[p]);
        }
    }

    for (int p = 0; p < pair_count; ++p)
    {
        char final_seq[1024];
        snprintf(final_seq, sizeof(final_seq), "%s/read_without_primer.fastq", result_output_prefix);

        FILE *fo_seq = fopen(final_seq, "wb");
        if (!fo_seq)
        {
            fprintf(stderr, "Error opening final file\n");
            return 1;
        }

        for (int t = 0; t < num_threads; ++t)
        {
            char tmp_seq[1024];
            snprintf(tmp_seq, sizeof(tmp_seq), "%s/pair%d_seq_t%d.fastq", result_output_prefix, p + 1, t);
            concat_file_to(fo_seq, tmp_seq);
            remove(tmp_seq);

            char tmp_len[1024];
            snprintf(tmp_len, sizeof(tmp_len), "%s/pair%d_length_t%d.txt", result_output_prefix, p + 1, t);
            remove(tmp_len); 
        }
        fclose(fo_seq);
    }

    int **global_pair_hist = (int **)malloc(pair_count * sizeof(int *));
    for (int p = 0; p < pair_count; ++p)
        global_pair_hist[p] = (int *)calloc(HIST_BINS, sizeof(int));
    int *global_allpairs_hist = (int *)calloc(HIST_BINS, sizeof(int));
    int *primer_edit_records_sum = (int *)calloc(pair_count, sizeof(int));
    if (!global_pair_hist || !global_allpairs_hist || !primer_edit_records_sum)
    {
        fprintf(stderr, "Out of memory in merge.\n");
        return 1;
    }

    long long allpairs_total_matches = 0;
    for (int t = 0; t < num_threads; ++t)
    {
        WorkerArg *w = args[t];
        for (int p = 0; p < pair_count; ++p)
        {
            for (int i = 0; i < HIST_BINS; ++i)
                global_pair_hist[p][i] += w->local_pair_match_hist[p][i];
            primer_edit_records_sum[p] += w->primer_edit_records[p];
        }
        for (int i = 0; i < HIST_BINS; ++i)
            global_allpairs_hist[i] += w->local_allpairs_hist[i];

        allpairs_total_matches += w->local_allpairs_total_matches;
    }

    gettimeofday(&tv_end, NULL);
    double load_time_sec = elapsed_sec(tv_start, tv_after_read);
    double processing_time_sec = elapsed_sec(tv_proc_start, tv_after_proc);
    double total_time_sec = elapsed_sec(tv_start, tv_end);

    {
        char path_sum[1024];
        snprintf(path_sum, sizeof(path_sum), "%s/primer_match_summary.txt", result_output_prefix);
        FILE *fs = fopen(path_sum, "wb");
        if (!fs)
        {
            fprintf(stderr, "Error opening %s\n", path_sum);
            return 1;
        }

        fprintf(fs, "ReadsLoaded\t%zu\n", num_reads);
        fprintf(fs, "Threads\t%d\n", num_threads);
        fprintf(fs, "Pairs\t%d\n", pair_count);
        fprintf(fs, "AllPairsTotalMatches\t%lld\n", allpairs_total_matches);
        fprintf(fs, "LoadTime_sec\t%.3f\n", load_time_sec);
        fprintf(fs, "ProcessingTime_sec\t%.3f\n", processing_time_sec);
        fprintf(fs, "TotalTime_sec\t%.3f\n", total_time_sec);

        for (int p = 0; p < pair_count; ++p)
            fprintf(fs, "Pair_%d\t%d\n", p + 1, primer_edit_records_sum[p]);
        fclose(fs);
    }

    for (size_t i = 0; i < num_reads; ++i)
        free(reads[i].seq);
    free(reads);

    for (int i = 0; i < primer_num; i++)
        free(primers[i]);
    free(primers);

    for (int t = 0; t < num_threads; ++t)
    {
        WorkerArg *w = args[t];
        if (!w)
            continue;
        if (w->f_output_seq)
            free(w->f_output_seq);
        if (w->f_output_len)
            free(w->f_output_len);
        if (w->primer_edit_records)
            free(w->primer_edit_records);
        if (w->local_allpairs_hist)
            free(w->local_allpairs_hist);
        if (w->local_pair_match_hist)
        {
            for (int p = 0; p < pair_count; ++p)
                free(w->local_pair_match_hist[p]);
            free(w->local_pair_match_hist);
        }
        free(w);
    }
    free(args);
    free(threads);

    for (int p = 0; p < pair_count; ++p)
        free(global_pair_hist[p]);
    free(global_pair_hist);
    free(global_allpairs_hist);
    free(primer_edit_records_sum);


    printf("The valid reads ratio: %.2f%%\n",
           ((double)allpairs_total_matches / num_reads) * 100);

    return 0;
}
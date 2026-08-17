#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#define MAX_LINE_LEN 100000
#define LINES_PER_RECORD 4  

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr,
            "Usage: %s <input.fastq> <reads_per_group> <num_groups> <outdir>\n",
            argv[0]);
        return 1;
    }

    const char *in_fastq = argv[1];
    long reads_per_group = atol(argv[2]);
    int num_groups      = atoi(argv[3]);
    const char *outdir  = argv[4];

    long max_records = reads_per_group * num_groups;
    long max_lines   = max_records * LINES_PER_RECORD;

    if (mkdir(outdir, 0755) && errno != EEXIST) {
        perror("mkdir");
        return 1;
    }

    FILE *fin = fopen(in_fastq, "r");
    if (!fin) {
        perror("open input");
        return 1;
    }

    FILE **fouts = calloc(num_groups, sizeof(FILE *));
    if (!fouts) {
        perror("calloc");
        fclose(fin);
        return 1;
    }

    char line[MAX_LINE_LEN];
    long line_no   = 0;
    long seq_count = 0;
    int  record_buffer_size = LINES_PER_RECORD;
    char **record_buffer = malloc(record_buffer_size * sizeof(char *));
    if (!record_buffer) {
        perror("malloc");
        free(fouts);
        fclose(fin);
        return 1;
    }
    

    for (int i = 0; i < record_buffer_size; i++) {
        record_buffer[i] = malloc(MAX_LINE_LEN * sizeof(char));
        if (!record_buffer[i]) {
            perror("malloc for buffer line");
            for (int j = 0; j < i; j++) {
                free(record_buffer[j]);
            }
            free(record_buffer);
            free(fouts);
            fclose(fin);
            return 1;
        }
    }

    while (line_no < max_lines && fgets(line, sizeof(line), fin)) {
        int pos_in_record = line_no % LINES_PER_RECORD;
        
        strncpy(record_buffer[pos_in_record], line, MAX_LINE_LEN - 1);
        record_buffer[pos_in_record][MAX_LINE_LEN - 1] = '\0';  
        
        if (pos_in_record == LINES_PER_RECORD - 1) {
            if (seq_count >= max_records) break;
            int grp = seq_count / reads_per_group; 
            
            if (!fouts[grp]) {
                char path[MAX_LINE_LEN];
                snprintf(path, sizeof(path),
                         "%s/sub_reads.fastq", outdir);
                fouts[grp] = fopen(path, "w");
                if (!fouts[grp]) {
                    perror("open output");
                    break;
                }
            }
            
            for (int i = 0; i < LINES_PER_RECORD; i++) {
                fputs(record_buffer[i], fouts[grp]);
            }
            
            seq_count++;
        }
        
        line_no++;
    }

    for (int i = 0; i < record_buffer_size; i++) {
        free(record_buffer[i]);
    }
    free(record_buffer);

    for (int i = 0; i < num_groups; i++) {
        if (fouts[i]) fclose(fouts[i]);
    }
    free(fouts);
    fclose(fin);

    if (seq_count < max_records) {
        fprintf(stderr,
            "Warning: only %ld sequences written, expected %ld\n",
            seq_count, max_records);
        if (line_no < max_lines) {
            fprintf(stderr,
            "Note: input file may have incomplete FASTQ records at the end\n");
        }
    } 

    return 0;
}
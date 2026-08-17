#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#define MAX_LINE_LEN 100000
#define LINES_PER_RECORD 4  

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr,
            "Usage: %s <input.fastq> <reads_to_extract> <outdir>\n",
            argv[0]);
        return 1;
    }

    const char *in_fastq = argv[1];
    long reads_to_extract = atol(argv[2]);
    const char *outdir  = argv[3];

    long max_lines = reads_to_extract * LINES_PER_RECORD;

    if (mkdir(outdir, 0755) && errno != EEXIST) {
        perror("mkdir");
        return 1;
    }

    FILE *fin = fopen(in_fastq, "r");
    if (!fin) {
        perror("open input");
        return 1;
    }

    char outpath[MAX_LINE_LEN];
    snprintf(outpath, sizeof(outpath), "%s/sub_reads.fastq", outdir);
    FILE *fout = fopen(outpath, "w");
    if (!fout) {
        perror("open output");
        fclose(fin);
        return 1;
    }

    char line[MAX_LINE_LEN];
    long line_no = 0;
    long seq_count = 0;

    char *record_buffer[LINES_PER_RECORD];
    for (int i = 0; i < LINES_PER_RECORD; i++) {
        record_buffer[i] = malloc(MAX_LINE_LEN * sizeof(char));
        if (!record_buffer[i]) {
            perror("malloc for buffer line");
            for (int j = 0; j < i; j++) {
                free(record_buffer[j]);
            }
            fclose(fout);
            fclose(fin);
            return 1;
        }
    }

    while (line_no < max_lines && fgets(line, sizeof(line), fin)) {
        int pos_in_record = line_no % LINES_PER_RECORD;
        
        strncpy(record_buffer[pos_in_record], line, MAX_LINE_LEN - 1);
        record_buffer[pos_in_record][MAX_LINE_LEN - 1] = '\0';
        
        if (pos_in_record == LINES_PER_RECORD - 1) {
            for (int i = 0; i < LINES_PER_RECORD; i++) {
                fputs(record_buffer[i], fout);
            }
            seq_count++;
        }
        
        line_no++;
    }

    for (int i = 0; i < LINES_PER_RECORD; i++) {
        free(record_buffer[i]);
    }

    fclose(fout);
    fclose(fin);

    if (seq_count < reads_to_extract) {
        fprintf(stderr,
            "Warning: only %ld sequences written, expected %ld\n",
            seq_count, reads_to_extract);
        if (line_no < max_lines) {
            fprintf(stderr,
            "Note: input file may have incomplete FASTQ records at the end\n");
        }
    } 

    return 0;
}
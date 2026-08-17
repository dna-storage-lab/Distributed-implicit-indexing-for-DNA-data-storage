#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void deletion_random(int cover, const char *in_fastq, const char *outdir) {
    char buff[600] = {0};
    char subname[200] = {0};
    FILE *fseq = NULL;
    FILE *fsub = NULL;
    
    fseq = fopen(in_fastq, "rt");  
    if (fseq == NULL) {
        printf("Cannot open the file %s!\n", in_fastq);
        return;
    }
    
    int maxLines = 1000000;
    char ***fastq_records = malloc(maxLines * sizeof(char **));
    if (fastq_records == NULL) {
        printf("Memory allocation failed!\n");
        fclose(fseq);
        return;
    }
    
    int recordCount = 0;
    while (recordCount < maxLines) {
        fastq_records[recordCount] = malloc(4 * sizeof(char *));
        if (fastq_records[recordCount] == NULL) {
            printf("Memory allocation for record failed at record %d\n", recordCount);
            break;
        }
        
        int lines_stored = 0;
        for (int line_num = 0; line_num < 4; line_num++) {
            if (fgets(buff, sizeof(buff), fseq) == NULL) {
                break;
            }

            fastq_records[recordCount][line_num] = malloc(strlen(buff) + 1);
            if (fastq_records[recordCount][line_num] == NULL) {
                printf("Memory allocation for line failed\n");
                break;
            }
            strcpy(fastq_records[recordCount][line_num], buff);
            lines_stored++;
        }
        
        if (lines_stored < 4) {
            for (int i = 0; i < lines_stored; i++) {
                free(fastq_records[recordCount][i]);
            }
            free(fastq_records[recordCount]);
            break;
        }
        
        recordCount++;
    }
    
    fclose(fseq);
    srand(123);  
    sprintf(subname, "%s/sub_reads.fastq", outdir);  
    fsub = fopen(subname, "wt");
    if (fsub == NULL) {
        printf("Cannot open file %s for writing!\n", subname);
    }


    for (int j = 0; j < cover; j++) {
        int randomIndex = rand() % recordCount;
        for (int line_num = 0; line_num < 4; line_num++) {
            fprintf(fsub, "%s", fastq_records[randomIndex][line_num]);
        }
    }
    
    fclose(fsub);


    for (int i = 0; i < recordCount; i++) {
        for (int j = 0; j < 4; j++) {
            free(fastq_records[i][j]);
        }
        free(fastq_records[i]);
    }
    free(fastq_records);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <InputFastq> <NumReads> <SubDir>\n", argv[0]);
        printf("  <InputFastq>: Input FASTQ file path\n");
        printf("  <NumReads>: Number of reads to sample per subset\n");
        printf("  <SubDir>: Output directory for FASTQ files\n");
        return 1;
    }
    
    const char *in_fastq = argv[1]; 
    int numReads = atoi(argv[2]);  
    const char *subDir = argv[3];    

    deletion_random(numReads, in_fastq, subDir);
    
    return 0;
}
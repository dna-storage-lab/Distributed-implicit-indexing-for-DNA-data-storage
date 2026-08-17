#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

void split_by_threshold(const char *input_file, double threshold_ratio,
                        const char *output_file_up, const char *output_file_low) {

    int peak = 0, shift = 0, zf1 = 0;
    int overlap = 0, wid = 0, readlen = 0;
    char buff[500] = {0};
    int seq_counter = 1; 

    FILE *fpeak = fopen(input_file, "rt");
    if (fpeak == NULL) {
        printf("Cannot open the input file: %s\n", input_file);
        return;
    }

    FILE *farti_up = fopen(output_file_up, "wt");
    if (farti_up == NULL) {
        printf("Cannot open the output file: %s\n", output_file_up);
        fclose(fpeak);
        return;
    }

    FILE *farti_low = fopen(output_file_low, "wt");
    if (farti_low == NULL) {
        printf("Cannot open the output file: %s\n", output_file_low);
        fclose(fpeak);
        fclose(farti_up);
        return;
    }

    while (1) {
        memset(buff, 0, sizeof(buff));
        
        int ret = fscanf(fpeak, "%s %d %d %d %d %d %d", buff, &peak, &shift, &zf1, &overlap, &wid, &readlen);
        if (ret == EOF) break;
        if (ret != 7) {
            int c;
            while ((c = fgetc(fpeak)) != '\n' && c != EOF) {}
            continue;
        }


        double dynamic_threshold = threshold_ratio * readlen * 2.0;

        if ((double)peak >= dynamic_threshold) {

            fprintf(farti_up, "%s %d %d %d %d %d %d\n", buff, peak, shift, zf1, overlap, wid, readlen);
        } else {
            fprintf(farti_low, ">%d\n", seq_counter);

            fprintf(farti_low, "%s\n", buff);
            
            seq_counter++;  
        }
    }

    fclose(fpeak);
    fclose(farti_up);
    fclose(farti_low);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: ./FilterThreshold <sliding_result> <threshold_ratio> <output_file_up> <output_file_low>\n");
        printf("Example: ./FilterThreshold result.txt 0.32 up.txt low.fasta\n");
        return 1;
    }

    const char *sliding_result = argv[1];
    double threshold_ratio = atof(argv[2]);
    const char *output_file_up = argv[3];
    const char *output_file_low = argv[4];

    if (threshold_ratio < 0) {
        printf("threshold_ratio must be non-negative.\n");
        return 1;
    }

    split_by_threshold(sliding_result, threshold_ratio, output_file_up, output_file_low);
    return 0;
}
#include <stdio.h> 
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define N_SEG   64827
#define SEG_LEN 110
#define CODEWORD_LEN 48510
#define N_CODewORD 147


typedef struct {
    uint32_t cnt0;
    uint32_t cnt1;
} BitCount;


static int base_to_bits(char b, int *up, int *low) {
    switch (toupper((unsigned char)b)) {
        case 'A': *up = 0; *low = 0; return 0;
        case 'T': *up = 0; *low = 1; return 0;
        case 'G': *up = 1; *low = 0; return 0;
        case 'C': *up = 1; *low = 1; return 0;
        case 'N':
            return -1;
        default:
            return -1;
    }
}


static void bit_majority_decide(const BitCount *bc, char *bit, float *prob) {
    uint32_t total = bc->cnt0 + bc->cnt1;
    if (total == 0 || bc->cnt0 == bc->cnt1) {
        *bit  = 'X';
        *prob = 0.5f;  
    } else {
        if (bc->cnt1 > bc->cnt0) {
            *bit  = '1';
            *prob = (float)bc->cnt1 / (float)total; 
        } else {
            *bit  = '0';
            *prob = (float)bc->cnt0 / (float)total; 
        }
    }
}


static char bits_to_base(char ub, char lb) {
    if (ub == 'X' || lb == 'X') return 'X';
    if (ub == '0' && lb == '0') return 'A';
    if (ub == '0' && lb == '1') return 'T';
    if (ub == '1' && lb == '0') return 'G';
    if (ub == '1' && lb == '1') return 'C';
    return 'X';
}

static int ensure_dir(const char *path) { 
    struct stat st; 
    if (stat(path, &st) == 0) { 
        if (S_ISDIR(st.st_mode)) return 0; 
        fprintf(stderr, "Error: %s exists but is not a directory\n", path); 
        return -1; 
    } 
    if (mkdir(path, 0755) != 0) { 
        perror("mkdir"); 
        return -1; 
    } 
    return 0; 
}


static const char *six_patterns[16] = {
    "000000","000001","000010","010000",
    "000100","100000","001010","010100",
    "001000","101000","100100","100010",
    "010001","010010","100001","001001"
};


void getthre(double coverage, double thres){
    int peak = 0,shift = 0,zf1 = 0;
    int overlap = 0,wid = 0;
    int readlen = 0;
    char buff[500] = {0};
    char peakname[200] = {0},readname[200] = {0};
    
    FILE *fpeak = NULL;
    FILE *fread = NULL;

    sprintf(peakname, "./result/%.1fx/Correlation/corr_result.txt",coverage);
    sprintf(readname,"./result/%.1fx/Correlation/corr_upthres.txt",coverage);

    fpeak = fopen(peakname,"rt");
    if(fpeak == NULL){
        printf("ERROR : cannot open the file %s!\n",peakname);
        return;
    }
    fread = fopen(readname,"wt");
    if (!fread) {
        perror("open corr_upthres.txt for write");
        fclose(fpeak);
        return;
    }
    
    while(1){
        memset(buff,'\0',sizeof(buff));
        int nread = fscanf(fpeak,"%s %d %d %d %d %d %d",
                           buff,&peak,&shift,&zf1,&overlap,&wid,&readlen);
                           
        if (nread == EOF) break;
        if (nread != 7) {
            int c;
            while ((c = fgetc(fpeak)) != '\n' && c != EOF) {}
            continue;
        }

        int len = (int)strlen(buff);
        double corr_value = (double)peak / (2.0 * len);

        if(corr_value >= thres){
            fprintf(fread,"%s %d %d %d %d %d %d\n",
                    buff,peak,shift,zf1,overlap,wid,readlen);
        }
    }

    fclose(fread);
    fclose(fpeak);
}


void majority_vote_stats(double coverage, const char *codeword_path)
{

    char inputname[512] = {0};
    sprintf(inputname,"./result/%.1fx/Correlation/corr_upthres.txt",coverage);

    FILE *fin = fopen(inputname, "rt");
    if (!fin) {
        perror("open reads (corr_upthres.txt)");
        return;
    }

    BitCount (*upper_counts)[CODEWORD_LEN] =
        calloc(N_CODewORD, sizeof(*upper_counts));
    BitCount (*lower_counts)[CODEWORD_LEN] =
        calloc(N_CODewORD, sizeof(*lower_counts));
    if (!upper_counts || !lower_counts) {
        perror("alloc bit counts");
        fclose(fin);
        free(upper_counts);
        free(lower_counts);
        return;
    }


    char seq_buf[65536];
    int peak=0, shift=0, dir_ignored=0, overlap=0, wid=0, readlen=0;

    while (1) {
        int nread = fscanf(fin, "%65535s %d %d %d %d %d %d",
                           seq_buf, &peak, &shift,
                           &dir_ignored, &overlap, &wid, &readlen);
        if (nread == EOF) break;
        if (nread != 7) {
            int c;
            while ((c = fgetc(fin)) != '\n' && c != EOF) {}
            continue;
        }

        if (wid < 1 || wid > N_SEG) continue;
        int seg = wid - 1;

        size_t seqlen = strlen(seq_buf);
        for (size_t i = 0; i < seqlen; ++i) {
            char b = (char)toupper((unsigned char)seq_buf[i]);
            int up_bit, low_bit;

            if (base_to_bits(b, &up_bit, &low_bit) != 0) {
                continue;
            }

            int p = shift + (int)i;   
            if (p < 0 || p >= SEG_LEN) continue;

            size_t abspos = (size_t)seg * SEG_LEN + (size_t)p;
            if (abspos >= (size_t)N_CODewORD * CODEWORD_LEN) continue;

            int codeword_idx    = (int)(abspos / CODEWORD_LEN);
            int pos_in_codeword = (int)(abspos % CODEWORD_LEN);
            if (codeword_idx < 0 || codeword_idx >= N_CODewORD) continue;

            BitCount *uc = &upper_counts[codeword_idx][pos_in_codeword];
            BitCount *lc = &lower_counts[codeword_idx][pos_in_codeword];

            if (up_bit == 0) uc->cnt0++; else uc->cnt1++;
            if (low_bit == 0) lc->cnt0++; else lc->cnt1++;
        }
    }
    fclose(fin);


    char  (*upper_bit)[CODEWORD_LEN]  =
        malloc(N_CODewORD * CODEWORD_LEN * sizeof(char));
    char  (*lower_bit)[CODEWORD_LEN]  =
        malloc(N_CODewORD * CODEWORD_LEN * sizeof(char));
    float (*upper_prob)[CODEWORD_LEN] =
        malloc(N_CODewORD * CODEWORD_LEN * sizeof(float));
    float (*lower_prob)[CODEWORD_LEN] =
        malloc(N_CODewORD * CODEWORD_LEN * sizeof(float));

    if (!upper_bit || !lower_bit || !upper_prob || !lower_prob) {
        perror("alloc bit decision arrays");
        free(upper_counts);
        free(lower_counts);
        free(upper_bit);  free(lower_bit);
        free(upper_prob); free(lower_prob);
        return;
    }


    for (int cw = 0; cw < N_CODewORD; ++cw) {
        for (int pos = 0; pos < CODEWORD_LEN; ++pos) {
            upper_bit[cw][pos]  = 'X';
            lower_bit[cw][pos]  = 'X';
            upper_prob[cw][pos] = 0.5f;
            lower_prob[cw][pos] = 0.5f;
        }
    }


    for (int cw = 0; cw < N_CODewORD; ++cw) {
        for (int pos = 0; pos < CODEWORD_LEN; ++pos) {
            bit_majority_decide(&upper_counts[cw][pos],
                                &upper_bit[cw][pos],
                                &upper_prob[cw][pos]);
            bit_majority_decide(&lower_counts[cw][pos],
                                &lower_bit[cw][pos],
                                &lower_prob[cw][pos]);
        }
    }

    free(upper_counts);
    free(lower_counts);


    if (ensure_dir(codeword_path)) {
        free(upper_bit);  free(lower_bit);
        free(upper_prob); free(lower_prob);
        return;
    }

    char   *bits_merged = (char*)malloc(97020 * sizeof(char));
    double *prob_merged = (double*)malloc(97020 * sizeof(double));
    if (!bits_merged || !prob_merged) {
        perror("alloc merged buffers");
        free(bits_merged); free(prob_merged);
        free(upper_bit);  free(lower_bit);
        free(upper_prob); free(lower_prob);
        return;
    }

    for (int cw = 0; cw < N_CODewORD; ++cw) {
        for (int g = 0; g < 441; ++g) {
            int base_src = g * 110;    
            int base_out = g * 220;   

            for (int j = 0; j < 110; ++j) {

                int src_pos = base_src + j; 

                int idx_up = base_out + j; 
                bits_merged[idx_up] = upper_bit[cw][src_pos];
                prob_merged[idx_up] = (double)upper_prob[cw][src_pos];

                int idx_low = base_out + 110 + j; 
                bits_merged[idx_low] = lower_bit[cw][src_pos];
                prob_merged[idx_low] = (double)lower_prob[cw][src_pos];
            }
        }

        char wm_line[256];  
        char wm_buffer[98000]; 

        char wm_path[] = "configure/watermark.txt";
        FILE *fwm = fopen(wm_path, "rt");
        if (!fwm) {
            perror("open watermark");
            fprintf(stderr, "Warning: cannot open watermark file, skip watermark.\n");
        } else {
            int start_line = cw * 441;
            int end_line = (cw + 1) * 441;
            int total_bits = 0;
            int line_count = 0;

            wm_buffer[0] = '\0';

            for (int i = 0; i < start_line; i++) {
                if (!fgets(wm_line, sizeof(wm_line), fwm)) {
                    fprintf(stderr, "Error: watermark file has less than %d lines\n", start_line);
                    fclose(fwm);
                    break;
                }
            }

            if (!feof(fwm)) {
                for (int i = start_line; i < end_line && !feof(fwm); i++) {
                    if (!fgets(wm_line, sizeof(wm_line), fwm)) {
                        fprintf(stderr, "Error: cannot read line %d for cw=%d\n", i+1, cw+1);
                        break;
                    }
                    

                    size_t line_len = strcspn(wm_line, "\r\n");
                    wm_line[line_len] = '\0';

                    if (line_len < 220) {
                        fprintf(stderr, "Warning: line %d length < 220 for cw=%d\n", i+1, cw+1);
                    }

                    int copy_len = (line_len < 220) ? line_len : 220;
                    strncat(wm_buffer, wm_line, copy_len);
                    total_bits += copy_len;
                    line_count++;
                }
            }
            
            fclose(fwm);

            if (total_bits < 97020) {
                fprintf(stderr, "Warning: watermark total length %d < 97020 for cw=%d\n", total_bits, cw+1);
            }

            for (int k = 0; k < 97020 && k < total_bits; ++k) {
                char w = wm_buffer[k];
                if (w == '1' && bits_merged[k] != 'X') {
                    if (bits_merged[k] == '0') bits_merged[k] = '1';
                    else if (bits_merged[k] == '1') bits_merged[k] = '0';
                }
            }
        }


        char out_path[512];
        sprintf(out_path, "%s/%d.txt", codeword_path, cw + 1);
        FILE *fout = fopen(out_path, "wt");
        if (!fout) {
            perror("open output prob file");
            continue;
        }


        for (int g = 0; g < 16170; ++g) {
            int base6 = g * 6;
            double q[6];
            for (int i = 0; i < 6; ++i) {
                char   bit_char = bits_merged[base6 + i]; 
                double p_dec    = prob_merged[base6 + i];  

                double p1;
                if (bit_char == '1') {
                    p1 = p_dec;        
                } else if (bit_char == '0') {
                    p1 = 1.0 - p_dec;    
                } else {
                    p1 = 0.5;        
                }

                if (p1 == 0.0) p1 = 0.0001;
                if (p1 == 1.0) p1 = 0.9999;
                q[i] = p1;
            }

            double p_raw[16];
            double p[16];
            for (int idx = 0; idx < 16; ++idx) {
                const char *pat6 = six_patterns[idx];
                double prob_pat = 1.0;
                for (int i = 0; i < 6; ++i) {
                    if (pat6[i] == '1')
                        prob_pat *= q[i];
                    else
                        prob_pat *= (1.0 - q[i]);
                }
                p_raw[idx] = prob_pat;
                p[idx]     = prob_pat;
            }

 
            double sum = 0.0;
            for (int idx = 0; idx < 16; ++idx) sum += p[idx];
            if (sum <= 0.0) {
                for (int idx = 0; idx < 16; ++idx) p[idx] = 1.0 / 16.0;
            } else {
                for (int idx = 0; idx < 16; ++idx) p[idx] /= sum;
            }

            double p1_0 = p[0] + p[1] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7];
            double p1_1 = 1.0 - p1_0;

            double p2_0 = p[0] + p[1] + p[2] + p[3] + p[8] + p[9] + p[10] + p[11];
            double p2_1 = 1.0 - p2_0;

            double p3_0 = p[0] + p[1] + p[4] + p[5] + p[8] + p[9] + p[12] + p[13];
            double p3_1 = 1.0 - p3_0;

            double p4_0 = p[0] + p[2] + p[4] + p[6] + p[8] + p[10] + p[12] + p[14];
            double p4_1 = 1.0 - p4_0;

            double out_p0[4] = {p1_0, p2_0, p3_0, p4_0};
            double out_p1[4] = {p1_1, p2_1, p3_1, p4_1};

            for (int b = 0; b < 4; ++b) {
                double prob0 = out_p0[b];
                double prob1 = out_p1[b];
                if (prob0 == 0.000000) {
                    prob0 = 0.000100;
                }
                if (prob1 == 0.000000) {
                    prob1 = 0.000100;
                }
                if (prob0 == 1.000000) {
                    prob0 = 0.999900;
                }
                if (prob1 == 1.000000) {
                    prob1 = 0.999900;
                }

                fprintf(fout, "%.6f %.6f\n", prob0, prob1);
            }
        }
        fclose(fout);
    }

    free(bits_merged);
    free(prob_merged);
    free(upper_bit);  free(lower_bit);
    free(upper_prob); free(lower_prob);
}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <coverage> <codeword_path> <thres> \n", argv[0]);
        return 1;
    }
    
    double coverage = atof(argv[1]);
    const char *codeword_path = argv[2];
    double thres = atof(argv[3]);

    getthre(coverage, thres);
    majority_vote_stats(coverage, codeword_path);
    return 0;
}

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static void trim_line_end(char *s, size_t *len) {
    while (*len > 0 && (s[*len - 1] == '\n' || s[*len - 1] == '\r')) {
        s[--(*len)] = '\0';
    }
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <decoded_file_k> <original_file_all> <output_file> <k>\n", argv[0]);
        return 1;
    }

    const char *decoded_path  = argv[1];
    const char *original_path = argv[2];
    const char *out_path      = argv[3];
    
    long long k = atoll(argv[4]);
    FILE *f_dec = fopen(decoded_path, "r");
    if (!f_dec) { perror(decoded_path); return 1; }

    FILE *f_org = fopen(original_path, "r");
    if (!f_org) { perror(original_path); fclose(f_dec); return 1; }

    FILE *f_out = fopen(out_path, "a"); 
    if (!f_out) { perror(out_path); fclose(f_dec); fclose(f_org); return 1; }

    char *line_dec = NULL, *line_org = NULL;
    size_t cap_dec = 0, cap_org = 0;
    ssize_t r1, r2;

    r1 = getline(&line_dec, &cap_dec, f_dec);

    for (long long i = 0; i < k; ++i) {
        r2 = getline(&line_org, &cap_org, f_org);
        if (r2 < 0) {
            fprintf(stderr, "Error: Original file doesn't have line %lld\n", k);
            goto cleanup;
        }
    }

    if (r1 < 0) {
        fprintf(stderr, "Error: Decoded file %s is empty\n", decoded_path);
        goto cleanup;
    }

    size_t len1 = (size_t)r1;
    size_t len2 = (size_t)r2;
    trim_line_end(line_dec, &len1);
    trim_line_end(line_org, &len2);

    if (len1 != len2) {
        fprintf(stderr, "Error: Length mismatch at line %lld (%zu vs %zu)\n", k, len1, len2);
        goto cleanup;
    }

    int dist = 0;
    for (size_t j = 0; j < len1; ++j) {
        if (line_dec[j] != line_org[j]) ++dist;
    }

    fprintf(f_out, "%lld %d\n", k, dist);

cleanup:
    free(line_dec); free(line_org);
    fclose(f_dec); fclose(f_org);
    fclose(f_out);
    return 0;
}
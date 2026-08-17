#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <ctype.h>
#include <sys/stat.h>

#define MAX_FILENAME_LEN 1024
#define BUFFER_SIZE 8192
#define INITIAL_SEQ_CAPACITY 1000
#define INITIAL_FILE_CAPACITY 100

typedef struct {
    char* id;
    char* sequence;
    int length;
} Sequence;

int ends_with(const char* str, const char* suffix) {
    size_t str_len = strlen(str);
    size_t suffix_len = strlen(suffix);
    
    if (str_len < suffix_len) return 0;
    return strcmp(str + str_len - suffix_len, suffix) == 0;
}


void* safe_malloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    return ptr;
}

void* safe_realloc(void* ptr, size_t size) {
    void* new_ptr = realloc(ptr, size);
    if (!new_ptr && size > 0) {
        fprintf(stderr, "Memory reallocation failed\n");
        exit(1);
    }
    return new_ptr;
}

Sequence* read_fasta_file(const char* filename, int* seq_count) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Warning: Cannot open file %s\n", filename);
        *seq_count = 0;
        return NULL;
    }
    
    Sequence* sequences = NULL;
    int capacity = INITIAL_SEQ_CAPACITY;
    *seq_count = 0;
    
    sequences = (Sequence*)safe_malloc(capacity * sizeof(Sequence));
    
    char line[BUFFER_SIZE];
    char current_id[BUFFER_SIZE] = "";
    char* current_seq = NULL;
    int current_seq_len = 0;
    int current_seq_capacity = INITIAL_SEQ_CAPACITY;
    
    while (fgets(line, sizeof(line), file)) {

        size_t len = strlen(line);
        if (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[len-1] = '\0';
            if (len > 1 && line[len-2] == '\r') {
                line[len-2] = '\0';
            }
        }
        
        if (line[0] == '>') {
            if (current_seq) {
                if (*seq_count >= capacity) {
                    capacity *= 2;
                    sequences = (Sequence*)safe_realloc(sequences, capacity * sizeof(Sequence));
                }
                
                sequences[*seq_count].id = strdup(current_id);
                sequences[*seq_count].sequence = current_seq;
                sequences[*seq_count].length = current_seq_len;
                (*seq_count)++;
                
                current_seq = NULL;
                current_seq_len = 0;
                current_seq_capacity = INITIAL_SEQ_CAPACITY;
            }
            
            strcpy(current_id, line + 1);
        } else if (line[0] != '\0') {

            int line_len = strlen(line);
            if (!current_seq) {
                current_seq = (char*)safe_malloc(current_seq_capacity + 1);
            }
            
            if (current_seq_len + line_len + 1 > current_seq_capacity) {
                while (current_seq_len + line_len + 1 > current_seq_capacity) {
                    current_seq_capacity *= 2;
                }
                current_seq = (char*)safe_realloc(current_seq, current_seq_capacity + 1);
            }
            
            strcpy(current_seq + current_seq_len, line);
            current_seq_len += line_len;
        }
    }

    if (current_seq) {
        if (*seq_count >= capacity) {
            sequences = (Sequence*)safe_realloc(sequences, (capacity + 1) * sizeof(Sequence));
        }
        
        sequences[*seq_count].id = strdup(current_id);
        sequences[*seq_count].sequence = current_seq;
        sequences[*seq_count].length = current_seq_len;
        (*seq_count)++;
    }
    
    fclose(file);
    

    if (*seq_count == 0) {
        free(sequences);
        return NULL;
    }
    
    return sequences;
}


char* majority_merge(char** sequences, int seq_count, int seq_length) {
    char* merged = (char*)safe_malloc(seq_length + 1);
    
    for (int col = 0; col < seq_length; col++) {
        int dash_count = 0;
        

        typedef struct {
            int count;
            int first_index;  
        } CharInfo;
        
        CharInfo char_info[256];
        for (int i = 0; i < 256; i++) {
            char_info[i].count = 0;
            char_info[i].first_index = -1;
        }
  
        for (int seq_idx = 0; seq_idx < seq_count; seq_idx++) {
            char c = sequences[seq_idx][col];
            unsigned char uc = (unsigned char)c;
            
            if (c == '-') {
                dash_count++;
            } else {
                if (char_info[uc].count == 0) {
                    char_info[uc].first_index = seq_idx;
                }
                char_info[uc].count++;
            }
        }

        char max_base = '-';
        int max_count = 0;
        int max_first_index = seq_count;  
        
        for (int i = 0; i < 256; i++) {
            if (i != '-' && char_info[i].count > 0) {
                if (char_info[i].count > max_count) {
                    max_count = char_info[i].count;
                    max_base = (char)i;
                    max_first_index = char_info[i].first_index;
                } else if (char_info[i].count == max_count) {

                    if (char_info[i].first_index < max_first_index) {
                        max_base = (char)i;
                        max_first_index = char_info[i].first_index;
                    }
                }
            }
        }
        

        if (max_count * 2 >= dash_count) {
            merged[col] = max_base;
        } else {
            merged[col] = '-';
        }
    }
    
    merged[seq_length] = '\0';
    return merged;
}


char* remove_gaps(const char* sequence) {
    int len = strlen(sequence);
    char* result = (char*)safe_malloc(len + 1);
    int result_len = 0;
    
    for (int i = 0; i < len; i++) {
        if (sequence[i] != '-') {
            result[result_len++] = toupper(sequence[i]);
        }
    }
    result[result_len] = '\0';
    
    return result;
}


void process_fasta_file(const char* filename, FILE* output_handle, int* processed_count, int total_files) {
    int seq_count;
    Sequence* sequences = read_fasta_file(filename, &seq_count);
    
    if (!sequences || seq_count == 0) {
        if (sequences) free(sequences);
        (*processed_count)++;
        return;
    }

    int seq_length = sequences[0].length;
    for (int i = 1; i < seq_count; i++) {
        if (sequences[i].length < seq_length) {
            seq_length = sequences[i].length;
        }
    }
    

    char** seq_array = (char**)safe_malloc(seq_count * sizeof(char*));
    int* need_free = (int*)calloc(seq_count, sizeof(int));  
    
    for (int i = 0; i < seq_count; i++) {
        if (sequences[i].length > seq_length) {
            char* truncated = (char*)safe_malloc(seq_length + 1);
            strncpy(truncated, sequences[i].sequence, seq_length);
            truncated[seq_length] = '\0';
            seq_array[i] = truncated;
            need_free[i] = 1;
        } else {
            seq_array[i] = sequences[i].sequence;
            need_free[i] = 0;
        }
    }
    

    char* merged_seq = majority_merge(seq_array, seq_count, seq_length);
    
    char* final_seq = remove_gaps(merged_seq);

    const char* base_name = strrchr(filename, '/');
    if (base_name) {
        base_name++;
    } else {
        base_name = filename;
    }
    
    char file_id[MAX_FILENAME_LEN];
    strcpy(file_id, base_name);
    char* ext_pos = strstr(file_id, ".fasta");
    if (!ext_pos) ext_pos = strstr(file_id, ".fa");
    if (ext_pos) *ext_pos = '\0';
    

    fprintf(output_handle, "@%s\n", file_id);
    fprintf(output_handle, "%s\n", final_seq);
    fprintf(output_handle, "+\n");
    for (int i = 0; i < (int)strlen(final_seq); i++) {
        fprintf(output_handle, "F");
    }
    fprintf(output_handle, "\n");
    

    (*processed_count)++;
    if (*processed_count % 5000 == 0 || *processed_count == total_files) {
        printf("Processed: %d/%d\n", *processed_count, total_files);
    }
    

    for (int i = 0; i < seq_count; i++) {
        free(sequences[i].id);
        if (need_free[i]) {
            free(seq_array[i]);  
        }
        free(sequences[i].sequence);
    }
    free(sequences);
    free(seq_array);
    free(need_free);
    free(merged_seq);
    free(final_seq);
}


char** get_fasta_files(const char* folder, int* count) {
    DIR* dir;
    struct dirent* entry;
    char** files = NULL;
    int capacity = INITIAL_FILE_CAPACITY;
    *count = 0;
    
    files = (char**)safe_malloc(capacity * sizeof(char*));
    
    dir = opendir(folder);
    if (!dir) {
        fprintf(stderr, "Error: Cannot open directory %s\n", folder);
        return NULL;
    }
    
    while ((entry = readdir(dir)) != NULL) {
        if (ends_with(entry->d_name, ".fasta") || ends_with(entry->d_name, ".fa")) {
            if (*count >= capacity) {
                capacity *= 2;
                files = (char**)safe_realloc(files, capacity * sizeof(char*));
            }
            

            char fullpath[MAX_FILENAME_LEN];
            snprintf(fullpath, sizeof(fullpath), "%s/%s", folder, entry->d_name);
            files[*count] = strdup(fullpath);
            (*count)++;
        }
    }
    
    closedir(dir);
    return files;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s -i <input_folder> -o <output_fastq>\n", argv[0]);
        fprintf(stderr, "Example: %s -i ./aligned_fastas -o output.fastq\n", argv[0]);
        return 1;
    }
    
    char* input_folder = NULL;
    char* output_file = NULL;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_folder = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        }
    }
    
    if (!input_folder || !output_file) {
        fprintf(stderr, "Error: Missing required arguments\n");
        return 1;
    }
    

    int total_files;
    char** files = get_fasta_files(input_folder, &total_files);
    
    if (!files || total_files == 0) {
        fprintf(stderr, "Error: No FASTA files found in %s\n", input_folder);
        if (files) free(files);
        return 1;
    }


    FILE* output = fopen(output_file, "w");
    if (!output) {
        fprintf(stderr, "Error: Cannot create output file %s\n", output_file);
        for (int i = 0; i < total_files; i++) free(files[i]);
        free(files);
        return 1;
    }
    

    int processed_count = 0;
    for (int i = 0; i < total_files; i++) {
        process_fasta_file(files[i], output, &processed_count, total_files);
        free(files[i]);
    }

    fclose(output);
    free(files);
    

    return 0;
}
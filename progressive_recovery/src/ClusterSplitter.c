#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif

#define INITIAL_CAPACITY 16
#define HASH_SIZE 1000003

typedef struct {
    char *key;          
    char *prefix;       
    char **seq_names;   
    int seq_count;
    int seq_cap;
} Cluster;

typedef struct {
    char *name;         
    char *sequence;     
    int group_number;   
} FastaEntry;

typedef struct {
    char *rep_key;    
    char *member_name; 
} HRow;


static char *str_dup(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s);
    char *p = (char *)malloc(len + 1);
    if (!p) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    memcpy(p, s, len + 1);
    return p;
}

static void *xmalloc(size_t size) {
    void *p = malloc(size);
    if (!p) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    return p;
}

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p) {
        perror("realloc");
        exit(EXIT_FAILURE);
    }
    return p;
}

static char *read_line(FILE *fp) {
    if (!fp) return NULL;

    size_t cap = 1024;
    size_t len = 0;
    char *buf = (char *)xmalloc(cap);

    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (len + 1 >= cap) {
            cap *= 2;
            buf = (char *)xrealloc(buf, cap);
        }
        buf[len++] = (char)c;
        if (c == '\n') break;
    }

    if (len == 0 && c == EOF) {
        free(buf);
        return NULL;
    }

    buf[len] = '\0';
    return buf;
}

static void rstrip_newline(char *s) {
    if (!s) return;
    size_t len = strlen(s);
    while (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r')) {
        s[len - 1] = '\0';
        len--;
    }
}

static int split_tab(char *line, char ***fields_out) {
    int cap = 16;
    int count = 0;
    char **fields = (char **)xmalloc(sizeof(char *) * cap);

    char *start = line;
    char *p = line;

    while (1) {
        if (*p == '\t' || *p == '\0') {
            if (count >= cap) {
                cap *= 2;
                fields = (char **)xrealloc(fields, sizeof(char *) * cap);
            }
            if (*p == '\t') {
                *p = '\0';
                fields[count++] = start;
                start = p + 1;
            } else {
                fields[count++] = start;
                break;
            }
        }
        p++;
    }

    *fields_out = fields;
    return count;
}

static void free_split_fields(char **fields) {
    free(fields);
}

static void ensure_dir_exists_single(const char *path) {
    struct stat st;
    if (stat(path, &st) == 0) return;

    if (mkdir(path, 0777) != 0 && errno != EEXIST) {
        perror(path);
        exit(EXIT_FAILURE);
    }
}

static void make_dirs_recursive(const char *path) {
    if (!path || !*path) return;

    char *tmp = str_dup(path);
    size_t len = strlen(tmp);

    if (len > 0 && (tmp[len - 1] == '/' || tmp[len - 1] == '\\')) {
        tmp[len - 1] = '\0';
    }

    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/' || *p == '\\') {
            char old = *p;
            *p = '\0';
            if (*tmp) ensure_dir_exists_single(tmp);
            *p = old;
        }
    }
    if (*tmp) ensure_dir_exists_single(tmp);
    free(tmp);
}

static char *get_dirname(const char *path) {
    if (!path) return str_dup("");

    const char *slash1 = strrchr(path, '/');
    const char *slash2 = strrchr(path, '\\');
    const char *slash = slash1;
    if (slash2 && (!slash1 || slash2 > slash1)) slash = slash2;

    if (!slash) return str_dup("");

    size_t len = (size_t)(slash - path);
    char *dir = (char *)xmalloc(len + 1);
    memcpy(dir, path, len);
    dir[len] = '\0';
    return dir;
}

static void ensure_cluster_seq_capacity(Cluster *c) {
    if (c->seq_count >= c->seq_cap) {
        c->seq_cap = (c->seq_cap == 0) ? 4 : c->seq_cap * 2;
        c->seq_names = (char **)xrealloc(c->seq_names, sizeof(char *) * c->seq_cap);
    }
}

static void append_seq_name(Cluster *c, const char *seq_name) {
    ensure_cluster_seq_capacity(c);
    c->seq_names[c->seq_count++] = str_dup(seq_name);
}

typedef struct HashNode {
    char *key;
    int value;
    struct HashNode *next;
} HashNode;

static unsigned long hash_str(const char *str) {
    unsigned long hash = 5381;
    int c;
    while ((c = (unsigned char)*str++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash;
}

static HashNode **hash_create(void) {
    HashNode **table = (HashNode **)calloc(HASH_SIZE, sizeof(HashNode *));
    if (!table) {
        perror("calloc");
        exit(EXIT_FAILURE);
    }
    return table;
}

static void hash_put(HashNode **table, const char *key, int value) {
    unsigned long h = hash_str(key) % HASH_SIZE;
    HashNode *node = table[h];
    while (node) {
        if (strcmp(node->key, key) == 0) {
            node->value = value;
            return;
        }
        node = node->next;
    }

    HashNode *new_node = (HashNode *)xmalloc(sizeof(HashNode));
    new_node->key = str_dup(key);
    new_node->value = value;
    new_node->next = table[h];
    table[h] = new_node;
}

static int hash_get(HashNode **table, const char *key, int *found) {
    unsigned long h = hash_str(key) % HASH_SIZE;
    HashNode *node = table[h];
    while (node) {
        if (strcmp(node->key, key) == 0) {
            *found = 1;
            return node->value;
        }
        node = node->next;
    }
    *found = 0;
    return -1;
}

static void hash_free(HashNode **table) {
    if (!table) return;
    for (int i = 0; i < HASH_SIZE; i++) {
        HashNode *node = table[i];
        while (node) {
            HashNode *next = node->next;
            free(node->key);
            free(node);
            node = next;
        }
    }
    free(table);
}


static void free_clusters(Cluster *clusters, int cluster_count) {
    for (int i = 0; i < cluster_count; i++) {
        free(clusters[i].key);
        free(clusters[i].prefix);
        for (int j = 0; j < clusters[i].seq_count; j++) {
            free(clusters[i].seq_names[j]);
        }
        free(clusters[i].seq_names);
    }
    free(clusters);
}

static void free_fasta_entries(FastaEntry *entries, int count) {
    for (int i = 0; i < count; i++) {
        free(entries[i].name);
        free(entries[i].sequence);
    }
    free(entries);
}

static void free_hrows(HRow *hrows, int count) {
    for (int i = 0; i < count; i++) {
        free(hrows[i].rep_key);
        free(hrows[i].member_name);
    }
    free(hrows);
}


static void uc_trans_dict(const char *uc_path,
                          Cluster **clusters_out,
                          int *cluster_count_out) {
    FILE *fp = fopen(uc_path, "r");
    if (!fp) {
        perror(uc_path);
        exit(EXIT_FAILURE);
    }

    Cluster *clusters = NULL;
    int cluster_count = 0;
    int cluster_cap = 0;

    HRow *hrows = NULL;
    int h_count = 0;
    int h_cap = 0;

    HashNode **cluster_map = hash_create();

    char *line;
    while ((line = read_line(fp)) != NULL) {
        rstrip_newline(line);

        if (line[0] == '\0') {
            free(line);
            continue;
        }

        char **fields = NULL;
        int nfields = split_tab(line, &fields);
        if (nfields <= 0) {
            free_split_fields(fields);
            free(line);
            continue;
        }

        if (strcmp(fields[0], "S") == 0) {
            if (nfields > 8) {
                char *key = fields[8];
                int found = 0;
                (void)hash_get(cluster_map, key, &found);

                if (!found) {
                    if (cluster_count >= cluster_cap) {
                        cluster_cap = (cluster_cap == 0) ? INITIAL_CAPACITY : cluster_cap * 2;
                        clusters = (Cluster *)xrealloc(clusters, sizeof(Cluster) * cluster_cap);
                    }

                    clusters[cluster_count].key = str_dup(key);
                    clusters[cluster_count].prefix = str_dup(fields[1]);
                    clusters[cluster_count].seq_names = NULL;
                    clusters[cluster_count].seq_count = 0;
                    clusters[cluster_count].seq_cap = 0;

                    append_seq_name(&clusters[cluster_count], key);

                    hash_put(cluster_map, key, cluster_count);
                    cluster_count++;
                }
            }
        } else if (strcmp(fields[0], "H") == 0) {
            if (nfields > 9) {
                if (h_count >= h_cap) {
                    h_cap = (h_cap == 0) ? INITIAL_CAPACITY : h_cap * 2;
                    hrows = (HRow *)xrealloc(hrows, sizeof(HRow) * h_cap);
                }
                hrows[h_count].rep_key = str_dup(fields[9]);
                hrows[h_count].member_name = str_dup(fields[8]);
                h_count++;
            }
        }

        free_split_fields(fields);
        free(line);
    }
    fclose(fp);

    for (int i = 0; i < h_count; i++) {
        int found = 0;
        int idx = hash_get(cluster_map, hrows[i].rep_key, &found);
        if (found) {
            append_seq_name(&clusters[idx], hrows[i].member_name);
        }
    }

    free_hrows(hrows, h_count);
    hash_free(cluster_map);

    *clusters_out = clusters;
    *cluster_count_out = cluster_count;
}

static void fasta_process(const char *fasta_file_name,
                          FastaEntry **entries_out,
                          int *count_out,
                          HashNode ***fasta_map_out) {
    FILE *fp = fopen(fasta_file_name, "r");
    if (!fp) {
        perror(fasta_file_name);
        exit(EXIT_FAILURE);
    }

    FastaEntry *entries = NULL;
    int count = 0;
    int cap = 0;
    int group_number = 1;

    HashNode **fasta_map = hash_create();

    while (1) {
        char *header = read_line(fp);
        if (!header) break;

        char *sequence = read_line(fp);
        if (!sequence) {
            free(header);
            break;
        }

        rstrip_newline(header);
        rstrip_newline(sequence);

        if (count >= cap) {
            cap = (cap == 0) ? INITIAL_CAPACITY : cap * 2;
            entries = (FastaEntry *)xrealloc(entries, sizeof(FastaEntry) * cap);
        }

        const char *name = header[0] == '>' ? header + 1 : header;

        entries[count].name = str_dup(name);
        entries[count].sequence = str_dup(sequence);
        entries[count].group_number = group_number;

        hash_put(fasta_map, entries[count].name, count);

        count++;
        group_number++;

        free(header);
        free(sequence);
    }

    fclose(fp);

    *entries_out = entries;
    *count_out = count;
    *fasta_map_out = fasta_map;
}

static void cluster_cu(Cluster *clusters, int cluster_count,
                       FastaEntry *fasta_entries, int fasta_count,
                       HashNode **fasta_map,
                       const char *path_cu_more) {
    (void)fasta_count;

    make_dirs_recursive(path_cu_more);

    for (int i = 0; i < cluster_count; i++) {
        Cluster *c = &clusters[i];
        int seq_num = c->seq_count;      
        int total_len = seq_num + 1;  

        if (total_len > 2) {
            size_t path_len = strlen(path_cu_more) + strlen(c->prefix) + 64;
            char *out_path = (char *)xmalloc(path_len);

#ifdef _WIN32
            snprintf(out_path, path_len, "%s\\%s_%d.fasta", path_cu_more, c->prefix, seq_num);
#else
            snprintf(out_path, path_len, "%s/%s_%d.fasta", path_cu_more, c->prefix, seq_num);
#endif

            FILE *out = fopen(out_path, "w");
            if (!out) {
                perror(out_path);
                free(out_path);
                exit(EXIT_FAILURE);
            }

            for (int j = 0; j < seq_num; j++) {
                int found = 0;
                int fidx = hash_get(fasta_map, c->seq_names[j], &found);
                if (found) {
                    fprintf(out, ">%s\n", c->seq_names[j]);
                    fprintf(out, "%s\n", fasta_entries[fidx].sequence);
                }
            }

            fclose(out);
            free(out_path);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr,
                "Usage: %s <uclust_output.uc> <low_corr_reads.fasta> <cu_more_dir>\n",
                argv[0]);
        return 1;
    }

    const char *uclust_output = argv[1];
    const char *low_corr_reads = argv[2];
    const char *cu_more_dir = argv[3];

    Cluster *clusters = NULL;
    int cluster_count = 0;

    FastaEntry *fasta_entries = NULL;
    int fasta_count = 0;
    HashNode **fasta_map = NULL;

    uc_trans_dict(uclust_output, &clusters, &cluster_count);
    fasta_process(low_corr_reads, &fasta_entries, &fasta_count, &fasta_map);

    cluster_cu(clusters, cluster_count,
               fasta_entries, fasta_count, fasta_map,
               cu_more_dir);

    free_clusters(clusters, cluster_count);
    free_fasta_entries(fasta_entries, fasta_count);
    hash_free(fasta_map);

    return 0;
}
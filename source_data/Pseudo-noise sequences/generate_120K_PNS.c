/*
 * generate_120K_PNS.c
 *
 * Compile:
 *   gcc -o generate_120K_PNS generate_120K_PNS.c
 *
 * Run:
 *   ./generate_120K_PNS.c
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OUTPUT_FILE          "PNS-Pool-120K.txt"
#define ROW_BITS             320u
#define NUM_ROWS             119922u
#define BLOCK_ROWS           8192u
#define LINE_BYTES           (ROW_BITS + 1u)
#define IO_BUFFER_SIZE       (8u * 1024u * 1024u)
#define EXPECTED_BITS        ((uint64_t)NUM_ROWS * (uint64_t)ROW_BITS)
#define EXPECTED_FILE_BYTES  ((uint64_t)NUM_ROWS * (uint64_t)LINE_BYTES)

static const uint8_t INIT1[3] = {
    1, 0, 1
};

static const uint8_t INIT2[4] = {
    1, 0, 0, 1
};

static const uint8_t INIT3[23] = {
    1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1,
    0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1
};


static uint32_t pack_state(const uint8_t *init, unsigned order)
{
    uint32_t state = 0u;
    unsigned i;

    for (i = 0u; i < order; ++i) {
        state |= ((uint32_t)(init[i] & 1u)) << i;
    }
    return state;
}

static double current_time_seconds(void)
{
    struct timespec ts;

    if (timespec_get(&ts, TIME_UTC) == TIME_UTC) {
        return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
    }

    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static void print_progress(uint32_t completed_rows)
{
    const unsigned width = 50u;
    const double ratio = (double)completed_rows / (double)NUM_ROWS;
    unsigned filled = (unsigned)(ratio * (double)width);
    unsigned i;

    if (filled > width) {
        filled = width;
    }

    fputc('\r', stderr);
    fputc('[', stderr);
    for (i = 0u; i < width; ++i) {
        fputc(i < filled ? '=' : ' ', stderr);
    }
    fprintf(stderr,
            "] %6.2f%% %u/%u",
            ratio * 100.0,
            completed_rows,
            NUM_ROWS);
    fflush(stderr);
}

int main(void)
{
    FILE *fp = NULL;
    char *output_block = NULL;
    char *stdio_buffer = NULL;

    uint32_t state1 = pack_state(INIT1, 3u);
    uint32_t state2 = pack_state(INIT2, 4u);
    uint32_t state3 = pack_state(INIT3, 23u);

    const uint32_t mask1 = (1u << 3u) - 1u;
    const uint32_t mask2 = (1u << 4u) - 1u;
    const uint32_t mask3 = (1u << 23u) - 1u;

    uint32_t row_start;
    double start_time;
    double elapsed;

    start_time = current_time_seconds();

    fp = fopen(OUTPUT_FILE, "wb");
    if (fp == NULL) {
        perror("Failed to create wm.txt");
        return EXIT_FAILURE;
    }

    stdio_buffer = (char *)malloc(IO_BUFFER_SIZE);
    if (stdio_buffer != NULL) {
        if (setvbuf(fp, stdio_buffer, _IOFBF, IO_BUFFER_SIZE) != 0) {
            fprintf(stderr, "Warning: Failed to set large file buffer.\n");
            free(stdio_buffer);
            stdio_buffer = NULL;
        }
    }

    output_block = (char *)malloc((size_t)BLOCK_ROWS * (size_t)LINE_BYTES);
    if (output_block == NULL) {
        fprintf(stderr, "Error: Failed to allocate output buffer.\n");
        fclose(fp);
        free(stdio_buffer);
        return EXIT_FAILURE;
    }

    for (row_start = 0u; row_start < NUM_ROWS; row_start += BLOCK_ROWS) {
        uint32_t rows_in_block = NUM_ROWS - row_start;
        uint32_t row;

        if (rows_in_block > BLOCK_ROWS) {
            rows_in_block = BLOCK_ROWS;
        }

        for (row = 0u; row < rows_in_block; ++row) {
            char *line = output_block + (size_t)row * (size_t)LINE_BYTES;
            uint32_t column;

            for (column = 0u; column < ROW_BITS; ++column) {
                const uint32_t bit1 = (state1 >> 2u) & 1u;
                const uint32_t bit2 = (state2 >> 3u) & 1u;
                const uint32_t bit3 = (state3 >> 22u) & 1u;

                const uint32_t feedback1 =
                    ((state1 >> 2u) ^ state1) & 1u;
                const uint32_t feedback2 =
                    ((state2 >> 3u) ^ state2) & 1u;
                const uint32_t feedback3 =
                    ((state3 >> 22u) ^ (state3 >> 4u)) & 1u;

                line[column] = (char)('0' + (bit1 ^ bit2 ^ bit3));

                state1 = ((state1 << 1u) & mask1) | feedback1;
                state2 = ((state2 << 1u) & mask2) | feedback2;
                state3 = ((state3 << 1u) & mask3) | feedback3;
            }

            line[ROW_BITS] = '\n';
        }

        {
            const size_t bytes_to_write =
                (size_t)rows_in_block * (size_t)LINE_BYTES;
            const size_t bytes_written =
                fwrite(output_block, 1u, bytes_to_write, fp);

            if (bytes_written != bytes_to_write) {
                perror("\nFailed to write to wm.txt");
                free(output_block);
                fclose(fp);
                free(stdio_buffer);
                return EXIT_FAILURE;
            }
        }

        print_progress(row_start + rows_in_block);
    }

    if (fflush(fp) != 0) {
        perror("\nFailed to flush wm.txt");
        free(output_block);
        fclose(fp);
        free(stdio_buffer);
        return EXIT_FAILURE;
    }

    {
        const long file_bytes = ftell(fp);

        if (file_bytes < 0L) {
            perror("\nFailed to check wm.txt size");
            free(output_block);
            fclose(fp);
            free(stdio_buffer);
            return EXIT_FAILURE;
        }

        if ((uint64_t)file_bytes != EXPECTED_FILE_BYTES) {
            fprintf(stderr,
                    "\nError: Output file size anomaly, actual %ld bytes, "
                    "expected %llu bytes.\n",
                    file_bytes,
                    (unsigned long long)EXPECTED_FILE_BYTES);
            free(output_block);
            fclose(fp);
            free(stdio_buffer);
            return EXIT_FAILURE;
        }
    }

    if (fclose(fp) != 0) {
        perror("\nFailed to close wm.txt");
        free(output_block);
        free(stdio_buffer);
        return EXIT_FAILURE;
    }
    fp = NULL;

    free(output_block);
    free(stdio_buffer);

    elapsed = current_time_seconds() - start_time;

    fputc('\n', stderr);
    printf("Generation complete!\n");

    return EXIT_SUCCESS;
}
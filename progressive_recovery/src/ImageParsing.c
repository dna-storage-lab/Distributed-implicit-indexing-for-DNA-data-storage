#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_FILES 474
#define BIT_LEN 43168

int main(int argc, char *argv[])
{
    if (argc != 3) {
        printf("Usage: %s <input_directory> <output_image>\n", argv[0]);
        printf("Example: %s RecoveryBitstream recovered_image.jpg\n", argv[0]);
        return 1;
    }

    const char *input_dir = argv[1];      
    const char *output_file = argv[2];

    FILE *out = fopen(output_file, "wb");

    if (!out)
    {
        perror("Cannot open output file");
        return 1;
    }

    char filename[256];
    char *bit_buffer = (char *)malloc(BIT_LEN);

    if (!bit_buffer)
    {
        printf("Memory allocation failed\n");
        fclose(out);
        return 1;
    }

    for (int file_id = 1; file_id <= NUM_FILES; file_id++)
    {
        snprintf(
            filename,
            sizeof(filename),
            "%s/%d.txt",
            input_dir,
            file_id
        );

        FILE *fp = fopen(filename, "r");

        if (!fp)
        {
            printf("Cannot open %s\n", filename);
            fclose(out);
            free(bit_buffer);
            return 1;
        }

        int idx = 0;

        while (idx < BIT_LEN)
        {
            int c = fgetc(fp);

            if (c == EOF)
                break;

            if (c == '0' || c == '1')
            {
                bit_buffer[idx++] = (char)c;
            }
        }

        fclose(fp);
        
        if (idx != BIT_LEN)
        {
            printf(
                "%s length error: %d bits (expected %d)\n",
                filename,
                idx,
                BIT_LEN
            );

            fclose(out);
            free(bit_buffer);
            return 1;
        }
        unsigned char byte = 0;
        int bit_count = 0;

        for (int i = 0; i < BIT_LEN; i++)
        {
            byte <<= 1;

            if (bit_buffer[i] == '1')
            {
                byte |= 1;
            }

            bit_count++;

            if (bit_count == 8)
            {
                if (fwrite(&byte, 1, 1, out) != 1)
                {
                    printf("Error writing output file\n");
                    fclose(out);
                    free(bit_buffer);
                    return 1;
                }

                byte = 0;
                bit_count = 0;
            }
        }
    }

    fclose(out);
    free(bit_buffer);
    return 0;
}
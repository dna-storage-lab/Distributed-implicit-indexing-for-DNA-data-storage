/*
 * generate_64K_PNS.c
 *
 * Compile:
 *   gcc -o generate_64K_PNS generate_64K_PNS.c
 *
 * Run:
 *   ./generate_64K_PNS.c
 */


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OUTPUT_FILE "PNS-Pool-64K.txt"
#define ROW_BITS 220u
#define NUM_ROWS 64827u
#define BLOCK_ROWS 8192u
#define LINE_BYTES (ROW_BITS+1u)
#define TOTAL_BITS \
((uint64_t)ROW_BITS*(uint64_t)NUM_ROWS)

#define IO_BUFFER_SIZE \
(8u*1024u*1024u)

static const uint8_t INIT1[7]=
{
    1,0,0,
    1,0,0,1
};

static const uint8_t INIT2[8]=
{
    0,0,1,1,
    0,0,1,0
};

static const uint8_t INIT3[9]=
{
    1,0,1,
    1,0,1,
    1,0,1
};

static uint32_t pack_state(
        const uint8_t *init,
        unsigned order
)
{
    uint32_t state=0;


    for(unsigned i=0;i<order;i++)
    {
        state |=
        ((uint32_t)init[i])<<i;
    }


    return state;
}

static void progress(
        uint32_t rows
)
{
    const int width=50;
    double ratio =
    (double)rows/(double)NUM_ROWS;

    int filled =
    (int)(ratio*width);

    fprintf(stderr,"\r[");

    for(int i=0;i<width;i++)
    {
        fprintf(stderr,
        "%c",
        i<filled?'=':' ');
    }
    fprintf(stderr,
    "] %.2f%% %u/%u",
    ratio*100,
    rows,
    NUM_ROWS);

    fflush(stderr);
}


int main()
{
    FILE *fp;
    char *buffer;
    char *stdio_buffer;

    uint32_t state1 =
    pack_state(INIT1,7);

    uint32_t state2 =
    pack_state(INIT2,8);

    uint32_t state3 =
    pack_state(INIT3,9);

    uint32_t mask1 =
    (1u<<7)-1;

    uint32_t mask2 =
    (1u<<8)-1;

    uint32_t mask3 =
    (1u<<9)-1;

    fp=fopen(
        OUTPUT_FILE,
        "wb"
    );

    if(!fp)
    {
        perror("PNS-Pool-64K.txt");
        return -1;
    }


    stdio_buffer =
    malloc(IO_BUFFER_SIZE);

    if(stdio_buffer)
    {
        setvbuf(
            fp,
            stdio_buffer,
            _IOFBF,
            IO_BUFFER_SIZE
        );
    }

    buffer =
    malloc(
        BLOCK_ROWS*LINE_BYTES
    );

    if(!buffer)
    {
        printf("Memory error\n");
        return -1;
    }

    uint32_t row_start=0;
    while(row_start<NUM_ROWS)
    {
        uint32_t rows =
        NUM_ROWS-row_start;
        if(rows>BLOCK_ROWS)
            rows=BLOCK_ROWS;
        for(uint32_t r=0;r<rows;r++)
        {
            char *line =
            buffer+r*LINE_BYTES;

            for(uint32_t c=0;c<ROW_BITS;c++)
            {

                uint32_t b1 =
                (state1>>6)&1;

                uint32_t b2 =
                (state2>>7)&1;


                uint32_t b3 =
                (state3>>8)&1;

                uint32_t fb1 =
                ((state1>>6)
                 ^
                 (state1>>2))
                &1;

                uint32_t fb2 =
                ((state2>>7)
                ^
                (state2>>3)
                ^
                (state2>>2)
                ^
                (state2>>1))
                &1;

                uint32_t fb3 =
                ((state3>>8)
                ^
                (state3>>3))
                &1;

                line[c]=
                '0'+
                (b1^b2^b3);

                state1 =
                ((state1<<1)&mask1)
                |
                fb1;

                state2 =
                ((state2<<1)&mask2)
                |
                fb2;

                state3 =
                ((state3<<1)&mask3)
                |
                fb3;
            }
            line[ROW_BITS]='\n';
        }

        fwrite(
            buffer,
            LINE_BYTES,
            rows,
            fp
        );
        row_start+=rows;
        progress(row_start);

    }

    fflush(fp);
    fclose(fp);
    free(buffer);
    free(stdio_buffer);

    printf("\n\nGeneration complete!\n");

    return 0;
}
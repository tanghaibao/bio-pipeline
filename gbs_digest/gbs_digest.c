#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <argp.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    int l;
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s in.fq.gz CAGT\n", argv[0]);
        return 1;
    }

    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);
    char *p = 0;
    char *r = argv[2];

    int offset = 0;
    while ((l = kseq_read(seq)) >= 0)
    {
        p = strstr(seq->seq.s, r);
        if (p != NULL)
        {
            offset = p - seq->seq.s;
            printf("@%s\n%s\n", seq->name.s, p);
            printf("+\n%s\n", seq->qual.s + offset);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

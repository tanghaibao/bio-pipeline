#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <libgen.h>
#include <limits.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


int make_outfile(char *s, char *sep, char *r, char *outfile)
{
    char *d = dirname(s);
    char *f = basename(s);
    char *p = strstr(f, sep);
    int offset = p - f;
    sprintf(outfile, "%s/%.*s.%s%s", d, offset, f, r, p);

    return 0;
}

int main(int argc, char *argv[])
{
    gzFile fp, fw;
    kseq_t *seq;
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s motif in.fq.gz\n\n\t"
                "Split reads using motif (e.g. CAGT).\n", argv[0]);
        return 1;
    }

    int l;
    char *p = 0;
    char *r = argv[1];
    unsigned int a = 0, b = 0;
    unsigned int nreads = 0, nmotifs = 0, noutreads = 0;

    char *infile = argv[2];
    char outfile[PATH_MAX] = { 0 };
    make_outfile(infile, ".", r, outfile);

    fp = gzopen(infile, "r");
    fw = gzopen(outfile, "wb");
    seq = kseq_init(fp);

    while ((l = kseq_read(seq)) >= 0)
    {
        a = b = 0;
        p = seq->seq.s;
        nreads ++;
        while (p != NULL)
        {
            p = strstr(p, r);
            if (p == NULL)
            {
                b = seq->seq.l;
            }
            else
            {
                b = p - seq->seq.s;
                p ++;
                nmotifs ++;
            }
            int slen = b - a;
            if (slen >= 30)
            {
                gzprintf(fw, "@%s_%d_%d\n%.*s\n", seq->name.s, a, b, slen, seq->seq.s + a);
                gzprintf(fw, "+\n%.*s\n", slen, seq->qual.s + a);
                noutreads ++;
            }
            a = b;
        }
    }
    fprintf(stderr, "Output written to `%s`\n", outfile);
    fprintf(stderr, "Total input reads = %d\n", nreads);
    fprintf(stderr, "Total output reads = %d\n", noutreads);
    fprintf(stderr, "Total motifs = %d (avg %.1f per read)\n", nmotifs, 1. * nmotifs / nreads);

    kseq_destroy(seq);
    gzclose(fp);
    gzclose(fw);

    return 0;
}

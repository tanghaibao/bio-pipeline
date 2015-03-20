#include <string.h>
#include <zlib.h>
#include <libgen.h>
#include <limits.h>
#include <regex.h>
#include "cmdline.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


typedef struct args_info args_info;


static int endswith(const char *str, const char *suffix) {
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

static int make_outfile(char *s, char *sep, char *r, char *outfile) {
    char *d = dirname(s);
    char *f = basename(s);
    char *p = strstr(f, sep);
    int offset = p - f;
    sprintf(outfile, "%s/%.*s.GATC%s", d, offset, f, p);

    if (endswith(outfile, ".gz") == 0) {
        char *t = strdup(outfile);
        sprintf(outfile, "%s.gz", t);
        free(t);
    }

    return 0;
}

int main(int argc, char *argv[]) {
    gzFile fp, fw;
    kseq_t *seq;
    args_info args;
    if (parser(argc, argv, &args) != 0 || args.inputs_num != 2) {
        parser_print_help();
        exit(1);
    }

    int l;
    char *r = args.inputs[0];
    unsigned int a = 0, b = 0, g = 0, i = 0;
    unsigned int nreads = 0, nmotifs = 0, noutreads = 0;

    unsigned int max_matches = 10,  max_groups = 2;
    regex_t regex;
    regmatch_t groups[max_groups];
    int reti = 0;
    reti = regcomp(&regex, r, REG_EXTENDED);
    if (reti) {
        fprintf(stderr, "Could not compile regex %s\n", r);
        exit(1);
    }

    char *infile = args.inputs[1];
    char outfile[PATH_MAX] = { 0 };
    make_outfile(infile, ".", r, outfile);

    if (access(infile, F_OK) == -1) {
        fprintf(stderr, "File `%s` not found!\n", infile);
        return -1;
    }
    fp = gzopen(infile, "r");
    fw = gzopen(outfile, "wb");
    seq = kseq_init(fp);

    unsigned int m[100];  // stores the locations of motifs
    unsigned int offset = 0;
    char *cursor = 0;

    // http://stackoverflow.com/questions/2577193/how-do-you-capture-a-group-with-regex
    while ((l = kseq_read(seq)) >= 0) {
        nreads ++;
        int c = 0;
        m[c ++] = 0;
        cursor = seq->seq.s;
        for (i = 0; i < max_matches; i ++) {
            if (regexec(&regex, cursor, max_groups, groups, 0))
                break;

            offset = groups[0].rm_eo;
            for (g = max_groups - 1; g >= 0; g --) {
                if (groups[g].rm_so == (size_t) -1)
                    continue;  // no more groups
#if 0
                printf("Match %u, Group %u: [%2u-%2u]: %.*s\n",
                       i, g, groups[g].rm_so, groups[g].rm_eo,
                       groups[g].rm_eo - groups[g].rm_so, cursor + groups[g].rm_so);
#endif
                m[c ++] = groups[g].rm_so + cursor - seq->seq.s;
                nmotifs ++;
                break;  // stop when the group is matched
            }
            cursor += offset;
        }

        m[c ++] = seq->seq.l;
        for (i = 0; i < c - 1; i ++) {
            a = m[i];
            b = m[i + 1];
            int slen = b - a;
            if (slen >= args.minlen_arg) {
                if (seq->qual.s == NULL) { // FASTA format
                    gzprintf(fw, ">%s_%d_%d\n%.*s\n", seq->name.s, a, b, slen, seq->seq.s + a);
                }
                else {
                    gzprintf(fw, "@%s_%d_%d\n%.*s\n", seq->name.s, a, b, slen, seq->seq.s + a);
                    gzprintf(fw, "+\n%.*s\n", slen, seq->qual.s + a);
                }
                noutreads ++;
            }
        }
    }
    fprintf(stderr, "Output written to `%s`\n", outfile);
    fprintf(stderr, "Total input reads = %d\n", nreads);
    fprintf(stderr, "Total output reads = %d\n", noutreads);
    fprintf(stderr, "Total motifs = %d (avg %.1f per read)\n", nmotifs, 1. * nmotifs / nreads);

    kseq_destroy(seq);
    regfree(&regex);
    gzclose(fp);
    gzclose(fw);
    parser_free (&args);  // release allocated memory

    return 0;
}

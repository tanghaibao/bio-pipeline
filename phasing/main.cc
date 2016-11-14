#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include "vcf.h"


using namespace seqan;
using namespace std;


struct CPRA
{
    CharString seq;
    int32_t pos;
    CharString ref;
    CharString alt;
};

struct PhasingOptions
{
    CharString bam_file;
    CharString vcf_file;
};


ArgumentParser::ParseResult
parseCommandLine(PhasingOptions &opts, int argc, char const **argv)
{
    ArgumentParser parser("phasing");
    addArgument(parser, ArgParseArgument(
                        ArgParseArgument::STRING, "bam_file"));
    addArgument(parser, ArgParseArgument(
                        ArgParseArgument::STRING, "vcf_file"));

    addUsageLine(parser,
                 "[\\fBoptions\\fP] \\fBbam_file\\fP \\fBvcf_file\\fP");
    addDescription(parser,
                   "Parses BAM and VCF file to prepare phasing");

    setShortDescription(parser, "SeqAn-based phasing");
    setVersion(parser, "0.6.11");
    setDate(parser, __DATE__);

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(opts.bam_file, parser, 0);
    getArgumentValue(opts.vcf_file, parser, 1);

    return ArgumentParser::PARSE_OK;
}


int parse_vcf_file(vector<CPRA> variants, CharString &vcf_file)
{
    htsFile *inf = bcf_open(toCString(vcf_file), "r");
    if (inf == nullptr)
        return EXIT_FAILURE;

    bcf_hdr_t * hdr = bcf_hdr_read(inf);
    int nsamples = bcf_hdr_nsamples(hdr);
    if (nsamples != 1)
    {
        fprintf(stderr, "File %s contains %d samples\n", toCString(vcf_file), nsamples);
        return EXIT_FAILURE;
    }

    const char **seqnames = nullptr;
    int nseq;
    seqnames = bcf_hdr_seqnames(hdr, &nseq);

    bcf1_t *rec = bcf_init();
    if (rec == NULL)
        return EXIT_FAILURE;

    /* Conditions: bi-allelic SNPs */
    while (bcf_read(inf, hdr, rec) == 0)
    {
        bcf_unpack(rec, BCF_UN_STR);  // up to ALT inclusive
        if (rec->n_allele > 2)
            continue;
        if (strlen(rec->d.allele[0]) > 1)
            continue;
        if (strlen(rec->d.allele[1]) > 1)
            continue;
        fprintf(stderr, "seqname:%s pos:%d n_allele:%d ref:%s alt:%s\n",
                seqnames[rec->rid], rec->pos, rec->n_allele,
                rec->d.allele[0], rec->d.allele[1]);
    }

    // Cleanup
    if (seqnames != nullptr)
        free(seqnames);
    bcf_destroy(rec);
    bcf_close(inf);

    return 0;
}


int main(int argc, char const **argv)
{
    PhasingOptions opts;
    ArgumentParser::ParseResult res = parseCommandLine(opts, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    cout << "bam_file=" << opts.bam_file << endl;
    cout << "vcf_file=" << opts.vcf_file << endl;

    vector<CPRA> variants;
    parse_vcf_file(variants, opts.vcf_file);

    return 0;
}

#ifndef MAIN_H
#define MAIN_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include "htslib/vcf.h"
#include "htslib/sam.h"

using namespace std;
using namespace seqan;

struct CPRA
{
    CharString seq;
    int32_t pos;
    CharString ref;
    CharString alt;
    CPRA (const char *_seq, int32_t _pos, char *_ref, char *_alt)
    : seq(_seq), pos(_pos), ref(_ref), alt(_alt) { }
};

struct aux_t
{   // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
};

struct PhasingOptions
{
    CharString bam_file;
    CharString vcf_file;
};

ArgumentParser::ParseResult
parseCommandLine(PhasingOptions &opts, int argc, char const **argv);

#endif

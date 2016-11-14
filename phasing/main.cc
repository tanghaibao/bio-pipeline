#include <iostream>
#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>


using namespace seqan;
using namespace std;


struct CPRA
{
    CharString seq;
    int32_t pos;
    CharString ref;
    CharString alt;

    CPRA(VcfRecord &r)
    {
    }
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
    VcfFileIn vcfIn(toCString(vcf_file));
    // Attach to standard output.
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, cout, Vcf());

    // Copy over header
    VcfHeader header;
    readHeader(header, vcfIn);
    auto contigs = contigNamesCache(context(vcfIn));
    // for (auto &i : contigs) { cout << i << endl; }
    cout << "length(contigs)=" << length(contigs) << endl;

    VcfRecord record;
    CharString rName = "chr12";
    int rID = 0;
    while (!atEnd(vcfIn))
    {
        getIdByName(rID, contigs, rName);
        readRecord(record, vcfIn);
        if ((length(record.ref) > 1) || (length(record.alt) > 1))
            continue;
        cout << rID << "\t" << record.beginPos << "\t"
             << record.ref << "\t" << record.alt << endl;
    }

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

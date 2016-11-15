#include <seqan/arg_parse.h>
#include "main.h"


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

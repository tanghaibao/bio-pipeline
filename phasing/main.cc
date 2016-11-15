#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "main.h"


/* Collect all bi-allelic SNPs in a vcf_file
 * Conditions: 1) bi-allelic ; 2) SNPs */
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

    while (bcf_read(inf, hdr, rec) == 0)
    {
        bcf_unpack(rec, BCF_UN_STR);  // up to ALT inclusive
        if (rec->n_allele > 2)
            continue;
        if (strlen(rec->d.allele[0]) > 1)
            continue;
        if (strlen(rec->d.allele[1]) > 1)
            continue;
        variants.push_back(CPRA(seqnames[rec->rid], rec->pos,
                                rec->d.allele[0], rec->d.allele[1]));
    }

    for (auto &i : variants)
    {
        fprintf(stderr, "seqname:%s pos:%d ref:%s alt:%s\n",
                toCString(i.seq), i.pos, toCString(i.ref), toCString(i.alt));
    }

    // Cleanup
    if (seqnames != nullptr)
        free(seqnames);
    bcf_destroy(rec);
    bcf_close(inf);

    return EXIT_SUCCESS;
}


/* Function copied from `https://github.com/samtools/samtools/blob/develop/bam2depth.c`
 * This function reads a BAM alignment from one BAM file. */
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        break;
    }
    return ret;
}


/* For each variant position, check the read counts supporting each variant in
 * the bam_file */
int parse_bam_file(vector<CPRA> variants, CharString &bam_file)
{
    int tid, pos, n_plp;
    const bam_pileup1_t *ret;
    aux_t *data = (aux_t *) calloc(1, sizeof(aux_t));
    char *filename = toCString(bam_file);
    // Initialize BAM
    htsFormat *informat = (htsFormat *) calloc(1, sizeof(htsFormat));
    hts_parse_format(informat, filename);
    data->fp = sam_open_format(filename, "r", informat);
    if (data->fp == nullptr)
    {
        return EXIT_FAILURE;
    }
    data->hdr = sam_hdr_read(data->fp);
    if (data->hdr == nullptr)
    {
        fprintf(stderr, "Could not read header for `%s`\n", filename);
        return EXIT_FAILURE;
    }
    hts_idx_t *idx = sam_index_load(data->fp, filename);
    if (idx == nullptr)
    {
        fprintf(stderr, "Could not load index header for `%s`\n", filename);
        return EXIT_FAILURE;
    }

    // Pileup
    auto plp = bam_plp_init(read_bam, (void *) data);
    for (auto &i : variants)
    {
        //data->iter = sam_itr_querys(idx, data->hdr, reg);
        while ((ret=bam_plp_auto(plp, &tid, &pos, &n_plp)) != nullptr)
        {
        }
    }

    bam_plp_destroy(plp);

    return EXIT_SUCCESS;
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
    parse_bam_file(variants, opts.bam_file);

    return EXIT_SUCCESS;
}

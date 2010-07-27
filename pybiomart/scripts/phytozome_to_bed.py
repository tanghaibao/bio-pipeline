#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
This script queries the phytozome database to get gene annotations (chr, start, stop, gene name), commonly known as the .bed format.
"""


import urllib as U
import sys

xml = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName = "zome_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "phytozome" interface = "default" >
        <Filter name = "orgid" value = "%s" />
        <Attribute name = "chr_name1" />
        <Attribute name = "gene_chrom_start" />
        <Attribute name = "gene_chrom_end" />
        <Attribute name = "gene_name1" />
    </Dataset>
</Query>
""".replace("\n", "")
url = "http://www.phytozome.net/biomart/martservice"
# check available filters and attributes at following address
filters = "%s?type=filters;dataset=phytozome" % url
attributes = "%s?type=attributes;dataset=phytozome" % url

# species-id mapping
species = "Arabidopsis lyrata,Arabidopsis thaliana,Brachypodium distachyon,Carica papaya,Chlamydomonas reinhardtii,Cucumis sativus,Glycine max,Manihot esculenta,Medicago truncatula,Mimulus guttatus,Oryza sativa,Physcomitrella patens,Populus trichocarpa,Ricinus communis,Selaginella moellendorffii,Sorghum bicolor,Vitis vinifera,Zea mays".split(',')
ids = [107,130,114,113,136,122,109,137,135,118,120,77,129,119,91,79,68,121]
species_id = dict(zip(species, ids))


class PhytozomeMart(object):

    def __init__(self, xml, url, filter=ids):
        xml = xml % ",".join(str(x) for x in filter)
        self.xml = xml
        self.url = url

    def send_query(self):
        data = U.urlencode({"query": self.xml})
        return U.urlopen(self.url + "?" + data).read()


if __name__ == '__main__':

    mart = PhytozomeMart(xml, url, filter=[species_id["Arabidopsis thaliana"]])
    print mart.send_query()


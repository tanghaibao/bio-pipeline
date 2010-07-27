"""
Builds the queries for the BioMart service
Certain portion of the codes are ported from R package biomaRt (thanks)
"""

import sys
import urllib
from xml.etree.ElementTree import ElementTree, Element, SubElement, tostring


class MartXMLParser(ElementTree):

    def __init__(self, xml_data):
        self.parse(xml_data)

    def parse_marts(self):
        for t in self.getiterator("MartURLLocation"):
            if t.attrib["visible"]=="1":
                yield Mart(**t.attrib)

    def parse_configuration(self):
        # the attributes
        for t in self.getiterator("AttributeDescription"):
            yield Attribute(**t.attrib)

        # the filters 
        for t in self.getiterator("FilterDescription"):
            f = Filter(**t.attrib)
            options = [Option(**x.attrib) for x in t.getiterator("Option")]
            f.add_options(options)
            yield f


class Mart(dict):

    def __init__(self, host="www.biomart.org", path="/biomart/martservice",
            port="80", name="ensembl", virtual_schema="default", **attrib):

        self.__dict__ = attrib.copy()
        self.__dict__.update(x for x in locals().items() \
                                if x[0] not in ("self", "attrib"))

        self.registry = {}
        self.url = "http://%s:%s%s" % (self.host, self.port, path)
        self.display_name = self.__dict__.get("displayName", "")
        self.virtual_schema = self.__dict__.get("serverVirtualSchema", 
                self.virtual_schema)

    def __str__(self):
        return "\t".join((self.name, self.display_name, self.virtual_schema))
    
    def get_registry(self, archive=False):
        type = "registry_archive" if archive else "registry"
        params = urllib.urlencode(dict(type=type))
        xml_data = urllib.urlopen(self.url, params)

        parser = MartXMLParser(xml_data)
        for t in parser.parse_marts():
            self.registry[t.name] = t

    def list_registry(self):
        if len(self.registry)==0: self.get_registry()
        for m in sorted(self.registry.values()): 
            print m

    def get_datasets(self):
        params = urllib.urlencode(dict(type="datasets", mart=self.name))
        web_data = urllib.urlopen(self.url, params)

        for row in web_data:
            atoms = row.strip().split("\t")
            if atoms[0]=="TableSet":
                name, description, last_updated = atoms[1], atoms[2], atoms[-1]
                self[name] = Dataset(name, description, last_updated, self)

    def list_datasets(self):
        if len(self)==0: self.get_datasets()
        for m in sorted(self.values(), key=str): 
            print m


class Dataset(object):
    """ connect to a specified dataset in the database
    """
    def __init__(self, name, description, last_updated, mart):
        self.name = name
        self.description = description
        self.last_updated = last_updated
        self.mart = mart

        self.attributes = {} 
        self.filters = {} 

    def __str__(self):
        return "\t".join((self.name, self.description, self.last_updated))

    def get_configuration(self):
        params = urllib.urlencode(dict(type="configuration", dataset=self.name))
        xml_data = urllib.urlopen(self.mart.url, params)
        
        parser = MartXMLParser(xml_data)
        for t in parser.parse_configuration():
            if isinstance(t, Attribute):
                self.attributes[t.internalName] = t
            elif isinstance(t, Filter):
                self.filters[t.internalName] = t
                #print t.options.keys()

    def list_attributes(self):
        if len(self.attributes)==0: self.get_configuration()
        for m in sorted(self.attributes.values()):
            print m

    def list_filters(self):
        if len(self.filters)==0: self.get_configuration()
        for m in sorted(self.filters.values()):
            print m

    def query(self, filters={}, attributes=()):
        q = MartQuery(dataset=self)
        q.add_filters(**filters)
        q.add_attributes(attributes)
        return q.execute()


class MartQuery(object):
    
    def __init__(self, dataset=None, formatter="TSV", header="0",
            unique_rows="0", count="0"):
        self.dataset = dataset
        self.url = dataset.mart.url
        self.virtual_schema = dataset.mart.virtual_schema
        self.formatter = formatter 
        self.header = header 
        self.unique_rows = unique_rows 
        self.count = count 
        self.name = dataset.name
        self.attributes = [] 
        self.filters = {} 

    def add_filters(self, **filters):
        for key, val in filters.items():
            """
            if key not in self.dataset.filters:
                print >>sys.stderr, "[Warning] %s not in filters, ignored.." % key
            else:
                self.filters[key] = str(val)
            """
            self.filters[key] = str(val)

    def add_attributes(self, attributes):
        for key in attributes:
            """
            if key not in self.dataset.attributes:
                print >>sys.stderr, "[Warning] %s not in attributes, ignored.." % key
            else:
                self.attributes.append(key)
            """
            self.attributes.append(key)

    def set_header(self, flag):
        self.header = str(flag) 

    def set_formatter(self, format="TSV"):
        self.formatter = format 

    def build_query(self):
        query_t = Element("Query", dict(virtualSchemaName=self.virtual_schema,
            formatter=self.formatter, header=self.header, uniqueRows=self.unique_rows,
            count=self.count, datasetConfigVersion="0.6"))
        dataset_t = SubElement(query_t, "Dataset", dict(name=self.name,
            interface="default")) 
        for key, val in self.filters.items():
            filter_t = SubElement(dataset_t, "Filter", dict(name=key, value=val))
        for attribute in self.attributes:
            attribute_t = SubElement(dataset_t, "Attribute", dict(name=attribute))

        return tostring(query_t)

    def execute(self):
        xml_data = self.build_query()
        print xml_data
        data = urllib.urlencode(dict(query=xml_data))
        return urllib.urlopen(self.url, data)


class MartArgument(object):

    def __init__(self, **attrib):
        self.__dict__ = attrib.copy()

    def __str__(self):
        return self.__class__.__name__ + str(self.__dict__)
    

class Attribute(MartArgument):
    """ attributes define the values that we are retrieving
    For example, the gene start, stop, or chromosomes it belongs to
    """
    pass


class Filter(MartArgument):
    """ filters define a restriction on the query.
    For example, you can restrict output to all genes located on chr. 1
    then use the filter chromosome_name with value `1`
    """
    def add_options(self, options):
        self.options = dict((x.displayName, x) for x in options) 


class Option(MartArgument):
    """ available options to go with filter
    """
    pass


class Sequence(object):
    def __init__(self, seq):
        self.seq = seq

    def export_fasta(self):
        pass


def test_biomart():
    bm = Mart()
    bm.list_registry()
    bm.list_datasets()
    return bm

def test_ensembl_dataset():
    bm = Mart()
    ensembl = bm.registry["ensembl"]
    ensembl.get_datasets()
    dataset = ensembl["mmusculus_gene_ensembl"]
    return dataset

def test_phytozome_dataset():
    # either of the following method is okay
    #bm = Mart()
    #bm.get_registry()
    #phytozome = bm.registry["phytozome_mart"]

    # or
    phytozome = Mart(host="www.phytozome.net", port="80", 
            name="phytozome_mart", virtual_schema="zome_mart")

    phytozome.get_datasets()
    dataset = phytozome["phytozome"]
    return dataset


if __name__ == '__main__':
    
    # test_biomart()

    dataset = test_phytozome_dataset()
    #dataset.list_filters()
    #dataset.list_attributes()
    dataset.get_configuration()

    filters = dict(gene_name_filter1="AT5G54690")
    attributes = "chr_name1,gene_chrom_start,gene_chrom_end,gene_name1,go_id,ko_id,pfam_id".split(",") 
    print dataset.filters.keys()
    print dataset.attributes.keys()

    #q = MartQuery(dataset=dataset)
    #q.add_attributes(attributes)
    #q.add_filters(orgid=130)
    #data = q.execute()

    data = dataset.query(filters=filters, attributes=attributes)
    for row in data: print row,

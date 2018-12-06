import csv
import requests


with open('/sybig/home/kco/my_interproscan/interproscan-5.32-71.0/knirps_test.fasta_1.tsv') as fasta_output:
    for line in fasta_output:
        cols = line.split("\t")
        col_count = len(line.split("\t"))
        if col_count < 14:
            continue
        DB = 'iBB'
        DB_object_id = cols[0]
        DB_object_symbol = cols[0]
        Qualifier = ''
        GO_id = cols[13].split('|')
        DB_ref = 'GO_REF:0000002'
        Evidence = 'IEA'
        With = 'Interpro:' + cols[11]
        # For aspect information
        Aspect_url = 'http://ibeetle-base.uni-goettingen.de/ibp/go/GO2aspect/'
        Aspect = requests.get(Aspect_url)
        DB_object_name = ''
        DB_object_synonym = ''
        DB_object_type = 'Protein'
        Taxon = 'taxon:7070'
        Date = cols[10].replace('-', '')
        Assign_DB = 'iBB'
        Annotation_extension = ''
        Gene_product_form_id = ''

        for i, line in enumerate(fasta_output):
            with open('/sybig/home/kco/my_interproscan/interproscan-5.32-71.0/Knirps_test/knirps_GO_%i.tsv' %i, 'w') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(['!gaf-version: 2.0'])
                tsv_writer.writerow([DB, DB_object_id, DB_object_symbol, Qualifier, GO_id, DB_ref,
                                     Evidence, With, Aspect, DB_object_name, DB_object_synonym,
                                     DB_object_type, Taxon, Date, Assign_DB, Annotation_extension,
                                     Gene_product_form_id])


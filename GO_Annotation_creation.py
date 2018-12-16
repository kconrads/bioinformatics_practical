import requests
import csv

filename = 'knirps_test.fasta_1.tsv'
go_ID_list = {}
go_ID = {}

def get_Aspect(x):
    GO_id = x
    Aspect_url = 'http://ibeetle-base.uni-goettingen.de/ibp/go/GO2aspect/' + GO_id
    Aspect_raw = requests.get(Aspect_url).text
    Aspect_list = Aspect_raw.split(':')
    Aspect_long = Aspect_list[2].replace('"', '').replace("}", '')
    if Aspect_long == 'cellular_component':
        Aspect = 'C'
    elif Aspect_long == 'biological_process':
        Aspect = 'P'
    else: Aspect = 'F'
    return Aspect


with open(filename) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        # sometimes there is no go_term (i.e. column 13)
        try:
            go_term = cols[13]
        except:
            continue
        ID = cols[11]
        gene_name = cols[0]
        date = cols[10]
        # now loop through each go term and write one per line
        go_terms = go_term.split("|")
        for terms in go_ID:
            go_ID[terms] = ID
        go_ID_list[ID] = go_terms
        for k in go_ID_list:
            for x in go_ID_list[k]:
                go_ID[x] = k

save_name = gene_name + '_%i.tsv'

for i, k in enumerate(go_ID):
    with open(save_name % i, 'w') as out_file:
        DB = 'iBB'
        DB_object_id = gene_name
        DB_object_symbol = gene_name
        Qualifier = ''
        GO_id = k
        DB_ref = 'GO_REF:0000002'
        Evidence = 'IEA'
        With = 'Interpro:' + go_ID[k]
        # For aspect information
        Aspect = get_Aspect(k)
        DB_object_name = ''
        DB_object_synonym = ''
        DB_object_type = 'Protein'
        Taxon = 'taxon:7070'
        Date = ''.join(date.split('-')[::-1])
        Assign_DB = 'iBB'
        Annotation_extension = ''
        Gene_product_form_id = ''
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['!gaf-version: 2.0'])
        tsv_writer.writerow([DB, DB_object_id, DB_object_symbol, Qualifier, GO_id, DB_ref,
                             Evidence, With, Aspect, DB_object_name, DB_object_synonym,
                             DB_object_type, Taxon, Date, Assign_DB, Annotation_extension,
                             Gene_product_form_id])


#for item in go_ID:
#    print(go_ID[item])

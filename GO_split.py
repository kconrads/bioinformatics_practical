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
        GO_id_list = cols[13].split('|')
        GO_id_str = ','.join(GO_id_list)
        Aspect_url = 'http://ibeetle-base.uni-goettingen.de/ibp/go/GO2aspect/' + GO_id_str
        Aspect_raw = requests.get(Aspect_url).text
        Aspect_split = Aspect_raw.split(':')
        Aspect_count = len(Aspect_split)
        if Aspect_count < 3:
            continue
        Aspect_long = Aspect_split[2].replace('"', '').replace("}", '')
        print(Aspect_raw)













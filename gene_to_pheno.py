from flask import Flask
import json
import re
from collections import defaultdict
import requests
from mysql.connector import MySQLConnection, Error
from db_config import read_db_config

app = Flask(__name__)

phenotype = defaultdict(list)
genotype = defaultdict(list)
ontology_terms = []
developemental_stages = ['pupa', 'pupa_male', 'pupa_female', 'adult', 'adult_female', 'adult_male', 'egg', 'mature_egg',
                         'larva', 'L1_larva', 'late_larva_(L5+)', 'L7_larva', 'embryo', 'prepupa']

# Splits text at any special character, takes first output
def extract_morph_structure(x):
    morph_structure = re.split(r"[^-\w\s]", x)
    #print(x +" -> " +morph_structure[0])
    if morph_structure[0].strip() in get_ontology_terms():
        return morph_structure[0].strip()

# Only gets the ontolotgy if it hasn't already been done
def get_ontology_terms():
    if (len(ontology_terms) < 1):
        with open('fly_anatomy.obo') as fly_ontology:
            for line in fly_ontology:
                if (line.startswith('name:')):
                    ontology_terms.append(line.split(':')[1].strip())
    return ontology_terms

def get_TCN(FBgn):
    TCN_request_url = 'http://ibeetle-base.uni-goettingen.de/ibp/idmapper/FBgn2TC/' + FBgn
    TCN = requests.get(TCN_request_url).text
    TCN_list = json.loads(TCN)
    return TCN_list

def get_tron(TrOn):
    TrOn_request_url = 'http://smallsrv-vm:9998/tribolium/cls/' + TrOn
    TrOn_info = requests.get(TrOn_request_url).text
    TrOn_list = json.loads(TrOn_info)
    return TrOn_list

def get_tron_descriptors(tron_list):
    tron_descriptors = []
    for i in range(len(tron_list['annotations'])):
        if tron_list['annotations'][i]['name'] == 'label':
            tron_descriptors.append(tron_list['annotations'][i]['value'])
    return tron_descriptors

def query_with_fetchone(TCN):
    try:
        dbconfig = read_db_config()
        conn = MySQLConnection(**dbconfig)
        cursor = conn.cursor(prepared=True)
        stmt = '''SELECT f.ontologyID  
        FROM Feature f 
          JOIN CoreData cd ON cd.ID = f.coreData 
          JOIN CoreData_TriboliumGene cd2tg on cd.ID = cd2tg.coreData_ID JOIN 
          TriboliumGene tg ON cd2tg.triboliumGenes_ID = tg.ID 
        WHERE tg.tCNumber = %s'''
        cursor.execute(stmt, (TCN,))
        rows = cursor.fetchall()

        #for i in range(len(pheno_table[0])):
        # pheno_list.append(pheno_table[0][i][5].decode())
        output = []
        for row in rows:
            value =row[0]
            if (value is not None) and (value.decode() not in output):
                output.append(value.decode())
        return output

    except Error as error:
        print(error)

    finally:
        cursor.close()
        conn.close()


# Creating phenotype
with open('allele_phenotypic_data_fb_2018_05.tsv') as allele_pheno:
    for line in allele_pheno:
        if (line.startswith('#') or line.startswith('\n')): # skipping # and blankspace
            continue
        cols = line.split("\t")
        if (cols[1] not in phenotype): #if the entry doesn't exist, creates it
            phenotype[cols[1]] = []
        structure = extract_morph_structure(cols[2])
        if (structure is not None):
            phenotype[cols[1]].append(structure) #otherwise the list of pheno-traits appended

# Used defaultdict from collections
with open('fbal_to_fbgn_fb_2018_05.tsv') as allele_gene:
    for line in allele_gene:
        if (line.startswith('#') or line.startswith('\n')):
            continue
        cols = line.split("\t")
        genotype[cols[2]].append(cols[0])

@app.route('/')
def getStart():
    return 'Please enter FBgn!'


@app.route('/<string:FBgn>/')
def getPheno(FBgn):

    TCN_dict = get_TCN(FBgn)
    TCN_list = TCN_dict[FBgn]

    tron_list = []
    for tcn in TCN_list:
        tron_list.append(query_with_fetchone(tcn))

    tron_list2 = []
    for i in range(len(tron_list)):
        for j in tron_list[i]:
            tron_list2.append(get_tron(j))

    tron_descriptors = []
    for tron in tron_list2:
        tron_term = get_tron_descriptors(tron)[0]
        if tron_term is not None and tron_term not in tron_descriptors and tron_term not in developemental_stages:
            tron_descriptors.append(tron_term)

    geno_list = genotype[FBgn]

    pheno_list = []
    pheno_list_match = []
    for gene in geno_list:
        if phenotype.get(gene) is None:
            continue
        for structure in phenotype.get(gene):
            if structure is not None and structure not in pheno_list:
                pheno_list.append(structure)
                pheno_list_match.append(structure.replace("pupal ", "").replace("adult ", '').replace("larva ", "").replace("embryo ",''))

    match_list = []

    for tcas_structure in tron_descriptors:
        if tcas_structure not in developemental_stages:
            tcas_structure = tcas_structure.replace("pupal_", "").replace("adult_", '').replace("larva_", "").replace("embryo_",'')
            if tcas_structure in pheno_list_match:
                match_list.append(tcas_structure)

    data = {'FBgn_morph_pheno': pheno_list,
            'FBgn_to_TCN': TCN_list,
            'TrOn_list': tron_list,
            'TrOn_descriptors': tron_descriptors,
            'Common': match_list}
    response = app.response_class(
        response=json.dumps(data),
        status=200,
        mimetype='application/json'
    )

    return response

if __name__ == '__main__':
    app.run(debug=True)

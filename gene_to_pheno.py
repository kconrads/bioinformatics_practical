from flask import Flask
import json
import re
from collections import defaultdict
import requests
from mysql.connector import MySQLConnection, Error
from db_config import read_db_config
import gzip
import ftplib
from ftplib import FTP
import os
import fnmatch

app = Flask(__name__)

## File Update
with FTP('ftp.flybase.org') as ftp:
    try:
        ftp.login()
        ftp.cwd('/releases/current/')
        files = ftp.nlst()
        current_name = files[0].strip()

        with open('current_release.txt', 'r+') as current_release:
            current_dl_name = current_release.readlines()[-1].strip()

            if current_dl_name != current_name:
                print('Updating necessary files...')

                ftp.cwd('precomputed_files/alleles/')
                alleles_files = ftp.nlst()

                for f in alleles_files:
                    if fnmatch.fnmatch(f, 'allele_phenotypic_data_fb_*'):
                        local_dir = os.path.dirname(os.path.realpath(__file__))
                        local_filename = os.path.join(local_dir, f)

                        with open(local_filename, 'wb') as output:
                            print('Downloading ' + f + '...')
                            ftp.retrbinary('RETR %s' % f, output.write)
                            print('Downloading complete!')

                    elif fnmatch.fnmatch(f, 'fbal_to_fbgn_fb_*'):
                        local_dir = os.path.dirname(os.path.realpath(__file__))
                        local_filename = os.path.join(local_dir, f)

                        with open(local_filename, 'wb') as output:
                            print('Downloading ' + f + '...')
                            ftp.retrbinary('RETR %s' % f, output.write)
                            print('Downloading complete!')

                ftp.cwd('../ontologies/')

                ontology_files = ftp.nlst()

                for f in ontology_files:
                    if fnmatch.fnmatch(f, 'fly_anatomy*'):
                        local_dir = os.path.dirname(os.path.realpath(__file__))
                        local_filename = os.path.join(local_dir, f)

                        with open(local_filename, 'wb') as output:
                            print('Downloading ' + f + '...')
                            ftp.retrbinary('RETR %s' % f, output.write)
                            print('Downloading complete!')

                current_release.write(current_name + '\n')

            else:
                print('Up to date!')

    except ftplib.all_errors as e:
        print('FTP error:', e)

phenotype = defaultdict(list)
genotype = defaultdict(list)
ontology_terms = set()
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
        with gzip.open('fly_anatomy.obo.gz','rt', encoding="utf-8") as fly_ontology:
            for line in fly_ontology:
                if (line.startswith('name:')):
                    ontology_terms.add(line.split(':')[1].strip())
    return ontology_terms

def get_FBgn(TCN):
    FBgn_request_url = 'http://ibeetle-base.uni-goettingen.de/ibp/idmapper/TC2FBgn/' + TCN
    FBgn = requests.get(FBgn_request_url).text
    FBgn_list = json.loads(FBgn)
    return FBgn_list

def get_TCN(FBgn):
    TCN_request_url = 'http://ibeetle-base.uni-goettingen.de/ibp/idmapper/FBgn2TC/' + FBgn
    TCN = requests.get(TCN_request_url).text
    TCN_list = json.loads(TCN)
    return TCN_list

def get_tron(TrOn):
    TrOn_request_url = 'http://smallsrv-vm:9998/tribolium/cls/' + TrOn
    TrOn_info = requests.get(TrOn_request_url, headers = {'Accept': 'application/json'}).text
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


with open('current_release.txt') as current_release:
    current_dl_name = current_release.readlines()[-1].strip()
    current_file_version = current_dl_name.split('FB')[1].strip()


# Creating phenotype
with gzip.open('allele_phenotypic_data_fb_' + current_file_version + '.tsv.gz', 'rt', encoding="utf-8") as allele_pheno:
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
with gzip.open('fbal_to_fbgn_fb_' + current_file_version + '.tsv.gz', 'rt', encoding="utf-8") as allele_gene:
    for line in allele_gene:
        if (line.startswith('#') or line.startswith('\n')):
            continue
        cols = line.split("\t")
        genotype[cols[2]].append(cols[0])

@app.route('/')
def getStart():
    return 'Please enter FBgn!'


@app.route('/FBgn2TCN/<string:FBgn>')
def get_TCN_Pheno(FBgn):

    TCN_dict = get_TCN(FBgn)
    TCN_list = TCN_dict[FBgn]

    tron_list = []
    for tcn in TCN_list:
        tron_list.append(query_with_fetchone(tcn))

    tron_list2 = []
    for i in range(len(tron_list)):
        for j in tron_list[i]:
            tron_list2.append(get_tron(j))

    tron_list_formatted = [j for i in tron_list for j in i]

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
                pheno_list_match.append(structure.replace("pupal ", "").replace("adult ", '').replace("larva ", "").replace("embryo ",'').replace("embryonic", ""))

    match_list = []

    for tcas_structure in tron_descriptors:
        if tcas_structure not in developemental_stages:
            tcas_structure = tcas_structure.replace("pupal_", "").replace("adult_", '').replace("larva_", "").replace("embryo_",'')
            if tcas_structure in pheno_list_match:
                match_list.append(tcas_structure)

    data = {'FBgn_morph_pheno': pheno_list,
            'FBgn_to_TCN': TCN_list,
            'TrOn_list': tron_list_formatted,
            'TrOn_descriptors': tron_descriptors,
            'Common': match_list}
    response = app.response_class(
        response=json.dumps(data),
        status=200,
        mimetype='application/json'
    )

    return response

@app.route('/TCN2FBgn/<string:TCN>')
def get_FBgn_Pheno(TCN):

    FBgn_dict = get_FBgn(TCN)
    FBgn_list = FBgn_dict[TCN]

    tron_list = query_with_fetchone(TCN)

    tron_list2 = []
    for tron in tron_list:
        tron_list2.append(get_tron(tron))

    '''tron_list_formatted = [j for i in tron_list for j in i]'''

    tron_descriptors = []
    for tron in tron_list2:
        tron_term = get_tron_descriptors(tron)[0]
        if tron_term is not None and tron_term not in tron_descriptors and tron_term not in developemental_stages:
            tron_descriptors.append(tron_term)

    geno_list = []
    for gene in FBgn_list:
        geno_list = genotype[gene]

    pheno_list = []
    pheno_list_match = []
    for gene in geno_list:
        if phenotype.get(gene) is None:
            continue
        for structure in phenotype.get(gene):
            if structure is not None and structure not in pheno_list:
                pheno_list.append(structure)
                pheno_list_match.append(structure.replace("pupal ", "").replace("adult ", '').replace("larva ", "").replace("embryo ",'').replace('embryonic', ''))

    match_list = []

    for tcas_structure in tron_descriptors:
        if tcas_structure not in developemental_stages:
            tcas_structure = tcas_structure.replace("pupal_", "").replace("adult_", '').replace("larva_", "").replace("embryo_",'')
            if tcas_structure in pheno_list_match:
                match_list.append(tcas_structure)

    data = {'TrOn_list': tron_list,
            'TrOn_descriptors': tron_descriptors,
            'TCN_to_FBgn': FBgn_list,
            'FBgn_morph_pheno': pheno_list,
            'Common': match_list}
    response = app.response_class(
        response=json.dumps(data),
        status=200,
        mimetype='application/json'
    )

    return response

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')

from flask import Flask, jsonify
import re
from collections import defaultdict
import requests

app = Flask(__name__)

phenotype = defaultdict(list)
genotype = defaultdict(list)
ontology_terms = []

# Splits text at any special character, takes first output
def extract_morph_structure(x):
    morph_structure = re.split(r"[^-\w\s]", x)
    if morph_structure[0] in get_ontology_terms():
        return morph_structure[0]

# Only gets the ontolotgy if it hasn't already been done
def get_ontology_terms():
    if (len(ontology_terms) < 1):
        with open('fly_anatomy.obo') as fly_ontology:
            for line in fly_ontology:
                if (line.startswith('name:')):
                    ontology_terms.append(line.split(':')[1].strip())
    return ontology_terms

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

    phenotype_output = list(filter(None, map(phenotype.get, genotype[FBgn])))

    if genotype.get(FBgn, 0) == 0:
        return jsonify('Gene not found!')
    elif not phenotype_output:
        return jsonify('No phenotype found for this gene')
    else:
        return jsonify(phenotype_output)

if __name__ == '__main__':
    app.run(debug=True)
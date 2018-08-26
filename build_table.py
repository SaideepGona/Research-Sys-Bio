import pandas as pd
import glob
import copy
import csv
import sys
import requests

# Builds TF tables based on given ChIP-Seq data

distance_threshold_upstream = -10000                         # Permissable default peak distances from transcription start site for inclusion
distance_threshold_downstream = 1000

if len(sys.argv) > 4:
    distance_threshold_upstream = int(sys.argv[4])
    distance_threshold_downstream = int(sys.argv[5])

enc_metadata = sys.argv[1]
geo_metadata = sys.argv[2]
all_genes_name = sys.argv[3]                                 # Path to file of all genes

def add_experiment_to_table(table, experiment):

    '''
    Takes in an annotation file defined as:
        File Name: Experiment_ID
        Top Line: TF Name
        Next lines: Genes of interest
    And adds to a given table based on annotation information

    '''
    # print(experiment)
    experiment_file = open(experiment, 'r')
    exp_acc = experiment.rstrip(".tfgenes")
    experiment_contents = experiment_file.readlines()
    experiment_contents = [x.rstrip("\n") for x in experiment_contents]
    # print(experiment_contents)

    tf = experiment_contents[0]
    genes = experiment_contents[1:]

    table[exp_acc+"-"+tf] = copy.copy(empty_row)
    for gene in genes:
        if gene not in all_genes:
            extra_genes.append(gene)
            continue
        # print(gene)
        # print(table[exp_acc+"-"+tf])
        table[exp_acc+"-"+tf][gene] += 1

def dictionary_to_csv(table, thresholds):
    '''
    Converts dictionary table to CSV file format
    '''
    table_list = []
    table_list.append(["Experiment Accession", "Transcription Factor"]+all_genes)

    for tf, gene_bools in table.items():
        new_row = tf.split("-")
        for gene in all_genes:
            new_row.append(gene_bools[gene])
        table_list.append(new_row)
    print(len(table_list), len(table_list[1]))
    with open("tftable_"+thresholds[1]+"_"+thresholds[2]+".tsv", "w") as f:
        writer = csv.writer(f, delimiter = "\t")
        writer.writerows(table_list)

def get_request(accession):
    '''
    performs a get request given an experiment accession. returns experiment metadata as a dict
    '''

    headers = {'accept': 'application/json'}
    url = "http://www.encodeproject.org/biosample/" + accession + "/?frame=object"
    response = requests.get(url, headers = headers)
    response_dict = response.json()

    return response_dict

def process_annotation(accession_id, geo_matadata):

    '''
    Extracts information from annotation files annotation of a given experiment file.
    Output is a .tfgenes file named with accession id which also has the tf on first line and
    all genes with peaks on the second
    '''
    # print(accession_id + "_peaks.xls.anno", "adfs")
    if os.path.exists(accession_id+".tfgenes"):
        return
    annotation_df = pd.read_csv(accession_id + "_peaks.xls.anno", sep='\t')
    # print(annotation_df)
    threshold_file = open(accession_id+"_"+str(distance_threshold_upstream)+"_"+str(distance_threshold_downstream)+".tableprops", 'w') #Make tableprops file
    threshold_file.close()
    # print(annotation_df)
    # Insert filtering criteria below
    distance_df = (annotation_df[(annotation_df["Distance to TSS"] > distance_threshold_upstream) \
                        & (annotation_df["Distance to TSS"] < distance_threshold_downstream) \
                        ])
    # print(distance_df)
    all_genes = []
    for index, row in distance_df.iterrows():
        all_genes.append(str(row["Gene Name"]))

    # Get TF data for experiment

    if accession_id[0] == 'E':              # For encode data

        info = get_request(accession_id)
        target_string = info['target']
        target_string = target_string.lstrip('/targets/')
        target_string = target_string.rstrip('-human/')

    if accession_id[0] == 'G':              # For geo data

        target_string = geo_matadata[accession_id]

    write_list = [target_string] + all_genes
    # print('yo')
    # print(target_string, 't')
    out_file = open(accession_id + ".tfgenes", 'w')
    sys.stdout.flush()
    out_file.write("\n".join(write_list)+"\n")
    out_file.close()



def readin_geo(geo):
    '''
    Reads in a geo metadata file and converts to dictionary {exp_id:tf}
    '''

    exp_tfs = {}

    with open(geo, 'r') as g:
        for line in g:
            p_line = line.rstrip("\n").split("\t")
            gse = p_line[1]
            tf = p_line[4].lstrip("anti-").rstrip("-tag").rstrip("(Thermo").rstrip("(abcam").lstrip("eGFP-")
            exp_tfs[gse] = tf

    return exp_tfs


extra_genes = []
all_tsv_files = glob.glob("*.anno")         #Collect list of all local annotation files
all_annotation_files = []
for tsv in all_tsv_files:
    if tsv[0:3] == "ENC" or tsv[0:3] == "GSE":
        all_annotation_files.append(tsv)

geo_m = readin_geo(geo_metadata)
# print(geo_m)

for an_file in all_annotation_files:
    # print(an_file, 'an')
    try:
        process_annotation(an_file.rstrip("_peaks.xls.anno"), geo_m)
    except:
        # print(an_file)
        continue

all_tfgene_files = glob.glob("*.tfgenes")                            # Collects names of all contributing experiment annotation results
# print(all_tfgene_files)
table_properties = glob.glob("*.tableprops")                         # Finds the annotation thresholds for this current build
# print(table_properties)
thresholds = table_properties[0].rstrip(".tableprops").split("_")

all_genes_file = open(all_genes_name, "r")
all_genes = all_genes_file.readlines()                              # Collects list of all human genes
all_genes = [x.rstrip("\n") for x in all_genes]

empty_row = {}                                                      # Creates an instance of empty row
for gene in all_genes:
    empty_row[gene] = 0

full_table = {}                                                     # Starts a full table dictionary
for single_file in all_tfgene_files:
    print(single_file)
    add_experiment_to_table(full_table, single_file)

#print(full_table)
dictionary_to_csv(full_table, thresholds)
print(list(set(extra_genes)))

# Example Experiment Metadata
'''

{   '@id': '/experiments/ENCSR000EWS/',
    '@type': ['Experiment', 'Dataset', 'Item'],
    'accession': 'ENCSR000EWS',
    'aliases': [],
    'alternate_accessions': [],
    'assay_slims': ['DNA binding'],
    'assay_term_id': 'OBI:0000716',
    'assay_term_name': 'ChIP-seq',
    'assay_title': 'ChIP-seq',
    'assembly': ['hg19'],
    'award': '/awards/U54HG004558/',
    'biosample_summary': 'MCF-7',
    'biosample_synonyms': ['MCF7 cell', 'MCF7', 'MCF-7 cell'],
    'biosample_term_id': 'EFO:0001203',
    'biosample_term_name': 'MCF-7',
    'biosample_type': 'immortalized cell line',
    'category_slims': ['protein and DNA interaction'],
    'contributing_files': [],
    'date_created': '2013-12-12T07:13:41.707994+00:00',
    'date_released': '2012-05-14',
    'dbxrefs': ['UCSC-ENCODE-hg19:wgEncodeEH002812', 'GEO:GSM935445'],
    'description': 'GATA3 ChIP-seq on human MCF-7',
    'developmental_slims': [],
    'documents': [],
    'files': [   '/files/ENCFF000ZKP/',
                 '/files/ENCFF000ZKR/',
                 '/files/ENCFF000ZKS/',
                 '/files/ENCFF000ZKT/',
                 '/files/ENCFF000ZKW/',
                 '/files/ENCFF000ZLE/',
                 '/files/ENCFF001VQG/',
                 '/files/ENCFF002CZK/',
                 '/files/ENCFF000ZKV/',
                 '/files/ENCFF093ETA/'],
    'hub': '/experiments/ENCSR000EWS/@@hub/hub.txt',
    'internal_status': 'unreviewed',
    'internal_tags': ['ENCYCLOPEDIAv4'],
    'lab': '/labs/peggy-farnham/',
    'month_released': 'May, 2012',
    'objective_slims': [   'protein and DNA interaction identification '
                           'objective',
                           'epigenetic modification identification objective'],
    'organ_slims': ['mammary gland'],
    'original_files': [   '/files/ENCFF000ZKP/',
                          '/files/ENCFF000ZKR/',
                          '/files/ENCFF000ZKS/',
                          '/files/ENCFF000ZKT/',
                          '/files/ENCFF000ZKU/',
                          '/files/ENCFF000ZKW/',
                          '/files/ENCFF000ZLE/',
                          '/files/ENCFF001VQG/',
                          '/files/ENCFF002CZK/',
                          '/files/ENCFF000ZKV/',
                          '/files/ENCFF093ETA/'],
    'possible_controls': ['/experiments/ENCSR000EWW/'],
    'references': [],
    'related_files': [],
    'related_series': [],
    'replicates': [   '/replicates/c84bc4d4-7cdf-4de2-a116-1f2ef3a2e9c0/',
                      '/replicates/71c91366-c4f3-469e-8306-6ce6d00bdbf3/'],
    'replication_type': 'isogenic',
    'revoked_files': ['/files/ENCFF000ZKU/'],
    'schema_version': '15',
    'status': 'released',
    'submitted_by': '/users/d48be354-153c-4ca8-acaa-bf067c1fd836/',
    'superseded_by': [],
    'supersedes': [],
    'system_slims': [],
    'target': '/targets/GATA3-human/',
    'type_slims': ['immunoprecipitation assay'],
    'uuid': 'eed3ee6f-7f59-499d-b9f5-e3ebefdd8d43'}

'''
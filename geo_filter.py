import pandas as pd
import xmltodict
import lxml.etree as etree
import xml.etree.ElementTree as ET
import os
from fuzzysearch import find_near_matches
from fuzzywuzzy import process
import sys

# OUTPUT SHOULD BE A TSV ONLY WITH USABLE EXPERIMENTS. CONTROLS AND REPLICATES ARE LABELED AS SUCH WITH A COLUMN ENTRY
#LOOKING FOR TF BINDING
    #<Characteristics tag="chip antibody">

filter_file = sys.argv[1]           # filtered_get.csv
geo_query_out = sys.argv[2]         # series_data.tsv

tf_tags = ['chip antibody', 'antibody']
illegal_tfs = [None,'None','none','RNA','IgG','Input','input']
#LOOKING FOR CONTROLS
    #IgG mentioned
    #Control mentioned
control_checks = ['IgG','Control','Input','input']
mods = [["b", "("],["a","."]]
for mod in mods:
    if mod[0] == "b":
        new = [(mod[1]+x) for x in control_checks]
        control_checks = control_checks + new
    elif mod[0] == "a":
        new = [(x+mod[1]) for x in control_checks]
        control_checks = control_checks + new
print(control_checks)
illegal_controls = ['ING','IPG']
miniML = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'
human_checks = ['Mus musculus']


# ref1 = "hg38"
# ref1_illegal = [None,"hg37","hg19"]
# ref2 = "GRCh38"
# ref2_illegal = [None,"GRCm38", "GRCh37"]

# ref19 = "hg19"
# ref19_illegal = [None,"hg38","hg37"]

# igg = "IgG"
# control = "Control"

def remove_tail(str):
    split = str.split(' ')
    return split[0]

    # char = root.findall(miniML+"Channel")
    # for node in char:
    #     print(node.text)
    #     if is_in(node.attrib['tag'], tf_tags):
    #         tf = node.text
    # return tf

def fuzzy_extract(qs, ls, threshold):
    '''
    fuzzy matches 'qs' in 'ls' and returns list of
    tuples of (word,index)
    '''
    for word, _ in process.extractBests(qs, (ls,), score_cutoff=threshold):
        # print('word {}'.format(word))
        for match in find_near_matches(qs, word, max_l_dist=1):
            match = word[match.start:match.end]
            # print('match {}'.format(match))
            index = ls.find(match)
            yield match

def not_in(string, anti_list):
    if anti_list == None:
        return True
    for x in anti_list:
        if x == string:
            return False
    return True

def is_in(string, c_list):
    for x in c_list:
        if x == string:
            return True
    return False

# ************************************ XML NODE CHECK FUNCTIONS *****************************************************************
def find_tf(root):
    '''
    Checks entire document for nodes which indicate which TF is being studied
    '''
    tf = None
    for node in root.iter():
        # print(node.tag)
        if node.tag == (miniML+'Characteristics') and is_in(node.attrib['tag'], tf_tags):
            tf = node.text
    return tf


def check_word_willegal(word,root,illegal,thresh):
    '''
    Given root node checks if word exists in a fuzzy sense over entire xml node. Can also provide any known 'illegal' matches
    Return True if word exists
    '''
    # print(thresh)
    for child in root.iter():                                                           # Check match in all
        for match in fuzzy_extract(word,child.text,thresh):
            # print(match)
            if not_in(match, illegal):                  # if the match is not illegal return true
                print(match)
                return True
    return False

# ************************************ END SERIES CHECK FUNCTIONS *************************************************************
# ************************************* SAMPLE CHECK FUNCTIONS ****************************************************************

def check_control(root):
    '''
    Given a sample xml node, checks if it corresponds to a control experiment
    '''
    is_control_loc = False
    for check in control_checks:
        is_control_loc = check_word_willegal(check,root,illegal_controls,50)
    # print(is_control)
    return is_control_loc

def sample_is_not_valid(root):
    # return False
    total_bool = False

    non_human = False
    for hcheck in human_checks:
        non_human = check_word_willegal(hcheck,root,illegal_controls,75)
    if non_human:
        total_bool = True

    return total_bool

# ************************************ END SAMPLE CHECK FUNCTIONS *************************************************************


full_table = pd.read_csv(geo_query_out)                                             # Read in final output of geo query
cwd = os.getcwd()
xml_path = cwd + '/metadata/'

current_series = None


# xml_full_path = xml_path + str(current_series) + '/' + str(current_series) + '_family.xml'
# current_xml = ET.parse(xml_full_path)
# current_root = current_xml.getroot()

#print(current_xml['Contributor'])
total_series = 1
correct_series = []
incorrect_series = []
igg_series = []
control_series = []

overcount = 0
valid_count = 0
tf_count = 0

new_table = [['iid','series','is_control','sra_path','tf']]

for index, row in full_table.iterrows():                                                            # Check references and

    try:    # CHECKS METADATA FOR ENTIRE SERIES AFTER ENCOUNTERING FIRST SERIES ELEMENT
        if row['series'] != current_series:
            print('new_exp')
            control_dict = {}
            tf = None                                                      # For each experiment DO something
            valid = True
            # sra = True
            total_series += 1
            current_series = row['series']

            xml_full_path = xml_path + str(current_series) + '/' + str(current_series) + '_family.xml'
            current_xml = ET.parse(xml_full_path)
            current_root = current_xml.getroot()                              # Gets root node of xml for parsing


            all_samples = current_root.findall(miniML+'Sample')
            if len(all_samples) > 50:
                overcount += 1
                continue
            for sample in all_samples:                                       # Loops through all 'Sample' nodes
                # print(sample)
                print(valid)
                if sample_is_not_valid(sample):                              # if a sample is not valid, neither is the series
                    valid = False
                    continue
                is_control = 'False'
                if check_control(sample):
                    is_control = 'True'
                # print(is_control)
                if is_control == 'False' and tf == None:         #If sample is not a control and tf isn't yet identified, find tf in replicate
                    # print('checking tf')
                    out_find = find_tf(sample)
                    if out_find != None:
                        tf = remove_tail(out_find).strip('\n').rstrip(',')

                current_iid = sample.attrib['iid']
                control_dict[current_iid] = is_control

            # print(tf)
            if is_in(tf, illegal_tfs):
                tf_pass = False                                 # If tf cannot be indentified, skip the whole experiment
            else:
                tf_pass = True
    except:
        print('err')
        tf_pass = False
        valid = False
        continue

    # print(valid)
    if not tf_pass:                                                         # Skip if no tf was found or if there was an error
        tf_count += 1
        continue
    if not valid:
        valid_count += 1
        continue
    print(row)
    if row['iid'] in control_dict and row['sra_ids'] != None:
        row_new = []
        row_new.append(str(row['iid']))
        row_new.append(str(row['series']))
        row_new.append(str(control_dict[row['iid']]))
        row_new.append(str(row['sra_ids']))
        row_new.append(str(tf))
        new_table.append(row_new)
    else:
        continue

print(len(new_table))
print(overcount, 'overcount')
if os.path.exists(filter_file):
    os.remove(filter_file)

with open(filter_file, 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in new_table)



# print(new_table)
        # current_study_ref = check_word_willegel(ref1,current_root,ref1_illegal)
        # current_study_ref = check_word_willegel(ref2,current_root,ref2_illegal)
        # current_study_ref_19 = check_word_willegel(ref2,current_root,ref19_illegal)
        # current_study_ref_igg = check_word_willegel(ref2,current_root,'AGG')

        # # for child in current_root.iter():                                                           # Check match in all
        # #     for match in fuzzy_extract(ref1,child.text,50):
        # #         print(match)
        # #         if not_in(match, ref1_illegal) :
        # #             current_study_ref = True
        # # for child in current_root.iter():
        # #     for match in fuzzy_extract(ref2,child.text,50):
        # #         print(match)
        # #         if not_in(match, ref2_illegal):
        # #             current_study_ref = True
        # # for child in current_root.iter():
        # #     for match in fuzzy_extract(ref2,child.text,50):
        # #         print(match)
        # #         if not_in(match, ref19_illegal):
        # #             current_study_ref_19 = True

        # # for child in current_root.iter():
        # #     for match in fuzzy_extract(igg,child.text,100):
        # #         print(match)
        # #         if not_in(match, ref19_illegal):
        # #             current_study_ref_igg = True
        # #print(current_xml)
        # # if current_study_ref:
        # #     correct_series.append(current_series)
        # # if current_study_ref_19:
        # #     incorrect_series.append(current_series)

        # if current_study_ref_igg:
        #     igg_series.append(current_series)






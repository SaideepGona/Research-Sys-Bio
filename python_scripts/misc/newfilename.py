# Given a bed file of genome locations, returns bed file of features and labels if training set

import os
import sys
import numpy as np
import pdb
pythonversion=sys.version_info[0]

if pythonversion<3:
        import urllib2
        from urllib2 import urlopen
        from urllib import urlencode
else:
        from urllib.request import urlopen
        from urllib.parse import urlencode

input_file = sys.argv[1]        # Name of input bed file
training = sys.argv[2]          # Boolean true if it is training data
range_val = sys.argv[3]         # Range around each peak to be considered
output_data = sys.argv[4]       # Name of output file

def convert_ranges(input_file, range_val):
    '''
    Takes as input a chip-seq bed file (must have summit info) but replaces the given interval with a
    new expanded one based on given input ranges.
    '''

    with open(input_file, "r") as bed_in:
        bed_lines = bed_in.read().split("\n")

    bed_lines_split = []
    for line in bed_lines:
        bed_lines_split.append(line.split("\t"))

    new_lines = []                                  # Consider edge cases
    for line in bed_lines_split:
        new_bottom = line[4] - range_val
        if new_bottom < 0:
            new_bottom = 0
        new_top = line[4] + range_val
        new_line = line[:]
        new_line[1] = new_bottom
        new_line[2] = new_top
        new_lines.append(new_line)

    with open(input_file+".ranged", "w") as bed_out:
        out_string = ""
        for line in new_lines:
            out_string += "\t".join(line) + "\n"
        bed_out.write(out_string.rstrip("\n"))

    return input_file + ".ranged"


def get_conservation(pos):
    url="https://genome.ucsc.edu/cgi-bin/hgTables?"
    form_data={'hgsid':'666516159_sXwdmHvVwfpAowRRzGB5IcIEdQNY',
                         'clade':'mammal',
                          'org':'human',
                          'db':'hg38',
                          'hgta_group':'compGeno',
                          'hgta_track':'cons100way',
                          'hgta_table':'phyloP100way',
                          'hgta_regionType':'range',
                          'position':pos,
                          'hgta_outputType':'wigData',
                          'hgta_outFileName':'1',
                          'hgta_doTopSubmit':'get output'}
        params=urlencode(form_data).encode('utf-8')
        response=urlopen(url,params)
        data=response.read()
        data=data.decode('utf-8')
        #pdb.set_trace()
        data=data.split('\n')
        data=data[1:]
        out=[]
        #pdb.set_trace()
        for i in data:
               if (len(i)>0) and (i[0]!="#")and (i[0]!='v'):
                     ii=i.split()
                     ii=[int(ii[0]),float(ii[1])]
                     out.append(ii)
    # scores = [x[1] for x in out]
    mean = sum(scores)/len(scores)
        return out, mean

ranged_file = convert_ranges(input_file, range_val)

with open(ranged_file, 'r') as ranged:
    with open(output_data, 'w') as out:
        counter = 0
        for line in ranged:
            print(counter)
            counter += 1
            if counter == 1000:
                break
            split_l = line.rstrip('\n').split('\t')
            new_l = split_l[0:2]
            position_str = split_l[0]+':'+split_l[1]+'-'+split_l[2]
            cons, mean_cons = get_conservation(position_str)
            new_l += str(mean_cons) + '\n'
            out.write(new_l)




pos="chrY:10000-15000"  # format of genomic range
s=get_gc(pos) # output  (position,score)

print(s)


# def get_gc(pos):
#       url="https://genome.ucsc.edu/cgi-bin/hgTables?"
#       form_data={'hgsid':'666320311_Bf9AWctuHMXX5QeBndjSHJ95M99Z',
#                 'clade':'mammal',
#                 'org':'Human',
#                 'db':'hg38',
#                 'hgta_group':'map',
#                 'hgta_track':'gc5BaseBw',
#                 'hgta_table':'0',
#                 'hgta_regionType':'range',
#                 'position':pos,
#                 'hgta_outputType':'wigData',
#                 'hgta_outFileName':'1',
#                 'hgta_doTopSubmit':'get output'}
#       params=urlencode(form_data).encode('utf-8')
#       response=urlopen(url,params)
#       data=response.read()
#       data=data.decode('utf-8')
#       #pdb.set_trace()
#       data=data.split('\n')
#       data=data[1:]
#     print(data)
#       out=[]
#       #pdb.set_trace()
#       for i in data:
#            if (len(i)>0) and (i[0]!="#")and (i[0]!='v'):
#                ii=i.split()
#                ii=[int(ii[0]),float(ii[1])]
#                out.append(ii)
#       return out


# def get_data_tables(pos, form):
#     url="https://genome.ucsc.edu/cgi-bin/hgTables?"

#     # https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=666320311_Bf9AWctuHMXX5QeBndjSHJ95M99Z&clade=mammal&org=Human&db=hg38&hgta_group=compGeno&hgta_track=cons100way&hgta_table=phastCons100way&hgta_regionType=genome&position=chr1%3A11102837-11267747&hgta_outputType=wigData&hgta_outFileName=
#     # https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=666332037_HPrIAOA8MgMT3adqEDLSUXySSRcC&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=gc5BaseBw&hgta_table=0&hgta_regionType=genome&position=chr1%3A11102837-11267747&hgta_outputType=primaryTable&hgta_outFileName=

#     forms = {

#         "cons_form": {'hgsid':'666320311_Bf9AWctuHMXX5QeBndjSHJ95M99Z',
#                 'clade':'mammal',
#                 'org':'Human',
#                 'db':'hg38',
#                 'hgta_group':'compGeno',
#                 'hgta_track':'cons100way',
#                 'hgta_table':'phastCons100way',
#                 'hgta_regionType':'range',
#                 'position':pos,
#                 'hgta_outputType':'wigData',
#                 'hgta_outFileName':'1',
#                 'hgta_doTopSubmit':'get output'}counts how many coverage counts are 0

#         "gc_form": {'hgsid':'666320311_Bf9AWctuHMXX5QeBndjSHJ95M99Z',
#                 'clade':'mammal',
#                 'org':'Human',
#                 'db':'hg38',
#                 'hgta_group':'map',
#                 'hgta_track':'gc5BaseBw',
#                 'hgta_table':'0',
#                 'hgta_regionType':'range',
#                 'position':pos,
#                 'hgta_outputType':'wigData',
#                 'hgta_outFileName':'1',
#                 'hgta_doTopSubmit':'get output'}
#     }

#     params=urllib.parse.urlencode(forms[form]).encode("utf-8")
#     response=urlopen(url,params)
#     data=response.read()
#     data.decode('UTF-8')
#     split_data=data.split('\n')

#     return data


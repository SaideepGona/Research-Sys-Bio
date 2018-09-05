'''
Author: Saideep Gona

Given a bed file of genome locations, queries UCSC and returns a bed file with features and labels
NEEDS TO BE MODIFIED TO NOT USE WEB QUERYING
'''

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

# IO ************************************************************

if len(sys.argv) > 1:

    input_file = sys.argv[1]        # Name of input bed file
    training = sys.argv[2]          # Boolean true if it is training data
    range_val = sys.argv[3]         # Max range around each peak to be considered
    range_val = int(range_val)
    tf_density_data = sys.argv[4]   # File with tf_density values
    start_line = sys.argv[5]        # Line of input file to start on
    end_line = sys.argv[6]          # Line of input file to end on

else:

    input_file = None
    training = True
    range_val = None
    tf_density_data = None
    start_line = None
    end_line = None

# IO END *********************************************************

output_data = start_line + "-" + end_line + "_" + input_file[-3:]
start_line = int(start_line)
end_line = int(end_line)

def convert_ranges(input_file, range_val):
    '''
    Takes as input a bed file but replaces the given interval with a
    new expanded one based on given input ranges.
    '''

    with open(input_file, "r") as bed_in:
        with open(input_file + ".ranged", 'w') as out:

            for in_line in bed_in:

                line = in_line.rstrip('\n').split('\t')
                new_bottom = int(line[1]) - range_val
                if new_bottom < 0:
                    new_bottom = 0
                new_top = int(line[2]) + range_val
                new_line = [line[0]]
                new_line.append(str(new_bottom))
                new_line.append(str(new_top))

                out.write("\t".join(new_line)+"\n")

    return input_file + ".ranged"

def position_mod(pos, direction, amount):
    '''
    Convert a given position by reducing/increasing(direction-less/more) the range on both sides
    by an (amount)
    pos="chrY:10000-15000"
    '''

    ranges = pos.split(":")
    chrom = ranges[0]
    ranges = ranges[1].split("-")
    low = int(ranges[0])
    high = int(ranges[1])
    cur_dif = abs(high-low)
    if direction == "less":
        if (2*amount) > cur_dif:
            print("range too small")
            return
        else:
            new_low = low + amount
            new_high = high - amount
            new_pos = chrom + ":" + str(new_low) + "-" + str(new_high)
            print(new_pos)
            return new_pos

    if direction == "more":
        new_low = low - amount
        if new_low < 0:
            new_low = 0
        new_high = high + amount
        new_pos = chrom + ":" + str(new_low) + "-" + str(new_high)
        return new_pos


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
    scores = np.array([x[1] for x in out])
    mean = np.mean(scores)
    var = np.var(scores)
    return out, mean, var

def get_gc(pos):
    url="https://genome.ucsc.edu/cgi-bin/hgTables?"
    form_data={'hgsid':'666320311_Bf9AWctuHMXX5QeBndjSHJ95M99Z',
            'clade':'mammal',
            'org':'Human',
            'db':'hg38',
            'hgta_group':'map',
            'hgta_track':'gc5BaseBw',
            'hgta_table':'gc5BaseBw',
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
    data=data[1:-1]
    percents = [int(x.split('\t')[1]) for x in data if len(x.split('\t'))>1]
    #pdb.set_trace()
    mean = (sum(percents)/len(percents))/100

    return mean

def get_histone(pos, histone):
    print('hist')
    url="https://genome.ucsc.edu/cgi-bin/hgTables?"
    form_data={'hgsid': '666806887_se33NbCPAKBp23Zk2OAYrQ0p4ZgT',
            'clade':'mammal',
            'org':'Human',
            'db':'hg38',
            'hgta_group':'regulation',
            'hgta_track':'wgEncodeRegMark'+histone,
            'hgta_table':'wgEncodeBroadHistoneGm12878'+histone+'StdSig',
            'hgta_regionType':'range',
            'position':pos,
            'hgta_outputType':'wigData',
            'hgta_outFileName':'1',
            'hgta_doTopSubmit':'get output'}
    params=urlencode(form_data).encode('utf-8')
    # print(params)
    response=urlopen(url,params)
    data=response.read()
    data=data.decode('utf-8')
    #pdb.set_trace()
    data=data.split('\n')
    if len(data) < 2:
        print('no out,', pos)
        return 0, 0
    data=data[1:]
    print('mid')
    # print(data)
    vals = np.array([float(x.split('\t')[-1]) for x in data if len(x.split('\t'))>1])
    if len(vals) == 0:
        return 0, 0
    # mean = (sum(vals)/len(vals))
    mean = np.mean(vals)
    var = np.var(vals)
    #pdb.set_trace()

    return mean, var

ranged_file = convert_ranges(input_file, range_val)
if not os.path.exists(ranged_file+".rand"):
    randomize = "sort -R " +ranged_file+" > "+ranged_file+".rand"
    os.system(randomize)


with open(ranged_file+".rand", 'r') as ranged:
    with open(output_data, 'w') as out:
        for i, line in enumerate(ranged):
            if not (start_line < i < end_line):
                continue
            print(i)

            split_l = line.rstrip('\n').split('\t')
            new_l = split_l[0:3]
            position_str = split_l[0]+':'+split_l[1]+'-'+split_l[2]

            y = [10000,9900,5000]
            for x in range(len(y)):

                cons, mean_cons, var_cons = get_conservation(position_mod(position_str,'less', 9900))
                gc = get_gc(position_mod(position_str,'less',5000))
                print('me')
                H3k4me1, H3k4me1_var = get_histone(position_mod(position_str,'less',y[x]), 'H3k4me1')
                # H3k4me3 = get_histone(position_mod(position_str,'less',x), 'H3k4me3', '666806887_se33NbCPAKBp23Zk2OAYrQ0p4ZgT')
                print('ac')
                H3k27ac, H3k27ac_var = get_histone(position_mod(position_str,'less',y[x]), 'H3k27ac')

                new_l.append(str(mean_cons))        # Conservation mean score
                new_l.append(str(var_cons))         # Conservation variance
                new_l.append(str(gc))               # GC percentage from surrounding 5000 each way
                new_l.append(str(H3k4me1))          # Histones
                new_l.append(str(H3k4me1_var))
                new_l.append(str(H3k27ac))
                new_l.append(str(H3k27ac_var))

            out.write("\t".join(new_l)+"\n")




























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


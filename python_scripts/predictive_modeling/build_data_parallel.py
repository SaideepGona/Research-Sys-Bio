'''
Author:Saideep Gona

Parallelization scheme for build_data.py
'''


import os
import sys
from multiprocessing import Process

input_file = sys.argv[1]        # Name of input bed file
training = sys.argv[2]          # Boolean true if it is training data
range_val = sys.argv[3]         # Max range around each peak to be considered
tf_density_data = sys.argv[4]   # File with tf_density values
spacing = sys.argv[5]
spacing = int(spacing)
start_line = sys.argv[6]
start_line = int(start_line)
instances = sys.argv[7]
instances = int(instances)


def run_build(start, end):
    process_core = "python3 build_data.py "+input_file+" "+training+" "+range_val+" "+tf_density_data+" "
    process_full = process_core +str(start)+" "+str(end)
    os.system(process_full)

start = start_line

for x in range(instances):

    end = start + spacing
    p = Process(target=run_build, args=(start,end))

    start = start + spacing
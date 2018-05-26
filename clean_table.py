import os
import sys

table = sys.argv[1]
out = sys.argv[2]

with open(table, 'w') as inp:

    for line in inp:
        p_line = line.rstrip('\n').split('\t')
        print(p_line)
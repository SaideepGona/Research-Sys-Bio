import sys
import os

in_file = sys.argv[1]
out_file = sys.argv[2]
col = sys.argv[3]
col = int(col)

zero_count = 0
nonzero_count = 0

with open(in_file, 'r') as inp:
    with open(out_file, 'w') as out:
        for line in inp:
            # print(line)
            split_l = line.rstrip('\n').split("\t")
            # print(split_l)
            print(split_l[col-1])
            if split_l[col-1] == '0':
                print('zero')
                zero_count += 1
                continue
            else:
                out.write("\t".join(split_l)+"\n")
                nonzero_count += 1

print(zero_count, 'zeros')
print(nonzero_count, 'nonzeros')

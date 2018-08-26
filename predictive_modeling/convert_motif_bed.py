import os
import sys

fasta = sys.argv[1]
fimo_count = sys.argv[2]

pwd = os.getcwd()
fimo_dir = pwd + "/fimo_"

with open(fasta, 'r') as fast:
    loc_dict = {}
    for line in fast:
        if line[0] == ">":
            split_l = line.lstrip(">").rstrip("\n").split(" ")
            ident = split_l[0]
            loc = split_l[1]
            loc_dict[ident] = loc.lstrip("range=")

print(loc_dict)

for count in range(fimo_count):
    cur_dir = fimo_dir + str(count) + "/"
    with open(cur_dir+"fimo.txt", 'r') as cur:
        with open(cur_dir+"fimo.bed", 'w') as out:
            for line in cur:

                split_l = line.rstrip("\n").split("\t")
                idnt = split_l[1]
                start = split_l[2]
                end = split_l[3]

                loc = loc_dict[idnt]
                chrom = loc.split(':')[0]
                posits = loc.split(':')[1]
                start_ref = posits('-')[0]
                end_ref = posits('-')[1]

                new_start = str(int(start) + int(start_ref))
                new_end = str(int(start) + int(end_ref))

                new_bed_line = chrom + "\t" + new_start + "\t" +new_end + "\n"

                out.write(new_bed_line)





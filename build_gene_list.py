import os
import sys

# Given a GFF file, constructs a list of all unique gene names from that file

input_gff = sys.argv[1]
output_list = sys.argv[2]

with open(input_gff, "r") as gff:
    all_genes = []
    for line in gff:
        if line[0] == "#":
            continue
        p_line = line.rstrip("\n").split(";")
        # print(p_line)
        for sect in p_line:
            sect_split = sect.split("=")
            # print(sect_split)
            if sect_split[0] == "gene":
                # print(sect)
                gene_name = sect_split[1]
                all_genes.append(gene_name)

all_genes = list(set(all_genes))
print(len(all_genes))

with open(output_list, "w") as out:

    for gene in all_genes:
        out.write(gene+"\n")
import sys
import os

in_file = sys.argv[1]
sub_count = sys.argv[2]
sub_count = int(sub_count)
def chunks(l,n):
    for i in range(0,len(l),n):
        yield l[i:i+n]

with open(in_file, "r") as inp:
    contents = inp.read()

    split = contents.split("\n\n")

    header = split[0:4]
    motifs = split[4:]

    motif_chunks = chunks(motifs, sub_count)

    for sub_part in range(sub_count):

        write = "\n\n".join(header)
        write += "\n\n" + "\n\n".join(next(motif_chunks))

        with open("motif_chunk_"+str(sub_part)+".meme",'w') as out:
            out.write(write)

for chunk in range(sub_count):

    command = "fimo -o fimo_"+str(chunk)+ " motif_chunk_"+str(chunk)+".meme " + "promoter_clean.fasta"
    os.system(command)




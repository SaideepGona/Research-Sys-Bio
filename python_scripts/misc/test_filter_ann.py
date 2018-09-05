import pandas
import sys
distance_threshold_upstream = -1000
distance_threshold_downstream = 100

annotation_df = pandas.read_csv(sys.argv[1], sep='\t')
distance_df = (annotation_df[(annotation_df["Distance to TSS"] > distance_threshold_upstream) \
                            & (annotation_df["Distance to TSS"] < distance_threshold_downstream) \
                            ])

all_genes = []
for index, row in distance_df.iterrows():
    all_genes.append(str(row["Gene Name"]))

write_list = ["A TF"] + all_genes

out_file = open("123" + ".tfgenes", 'w')
out_file.write("\n".join(write_list))

print(annotation_df)
print(distance_df)
print(len(annotation_df), len(distance_df))
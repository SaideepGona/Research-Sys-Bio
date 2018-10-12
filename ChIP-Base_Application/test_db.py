from flask_app import db
from flask_app import Peaks, ChIP_Meta, Query_History

# x = Peaks.query.filter(Peaks.chrom=="chrX").all() 
# x = Query_History.query.all() 
# x = ChIP_Meta.query.filter(ChIP_Meta.tissue_types=="brain").all() 
x = ChIP_Meta.query.all()
for entry in x:
    print(entry.date_posted)
    print(entry.log_p)
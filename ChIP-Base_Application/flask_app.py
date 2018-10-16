'''
Author: Saideep Gona

This is the core script for the ChIP-Base application which hosts large-scale ChIP-Seq data. It allows for parameter specification and generation of binary
tf-gene binding tables. Built using the python Flask framework.
'''
import os
from datetime import datetime
from multiprocessing import Process, Queue
import pickle
import time
import glob

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

from flask import Flask
from flask import render_template, send_file, redirect, url_for
from flask_wtf import FlaskForm
from flask_sqlalchemy import SQLAlchemy
from wtforms import StringField, PasswordField, BooleanField, SubmitField, FloatField, IntegerField
from wtforms.validators import DataRequired, NumberRange, Length
from wtforms.fields.html5 import IntegerRangeField 
from wtforms_html5 import AutoAttrMeta


app = Flask(__name__)
app.config['SECRET_KEY'] = 'sai-key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///ChIP_Base.db'

pwd = os.getcwd()

all_genes_path = pwd + "/static_lists/all_genes.txt"
all_tfs_path = pwd + "/static_lists/all_tfs.txt"

all_genes = []
with open(all_genes_path, "r") as ag:
    for line in ag:
        all_genes.append(line.rstrip("\n"))
print("NUMBER OF GENES: ", len(all_genes))

all_tfs = []
with open(all_tfs_path, "r") as at:
    for line in at:
        all_tfs.append(line.rstrip("\n"))
print("NUMBER OF TFS: ", len(all_tfs))

class DownloadFiles():

    def __init__(self):
        self.make_time = datetime.utcnow

    def collect_peaks(self, peak_dir):
        peak_files = glob.glob(pwd + "/" + peak_dir + "/*")
        unpath_peak_files = ["#".join(x.split("/")) for x in peak_files]
        peak_files_strip = [x.split("/")[-1] for x in peak_files]
        self.peak_files = peak_files
        self.unpath_peak_files = unpath_peak_files
        self.peak_files_strip = peak_files_strip
        self.num_peak_files = len(peak_files)

# Forms

class ParameterForm(FlaskForm):
    # log_p = StringField('-log(p) Value', validators=[DataRequired()])
    # fold_enrichment = StringField('Fold Enrichment', validators=[DataRequired()])
    # promoter_bool = BooleanField('Promoters')
    # enhancer_bool = BooleanField('Enhancers')

    class Meta(AutoAttrMeta):
        pass

    promoter_bool = BooleanField('Promoters')
    enhancer_bool = BooleanField('Enhancers')

    transcription_factors = StringField('Transcription Factors')
    tissue_types = StringField('Tissue Types')

    log_p = FloatField('-log(p) Value', validators=[DataRequired(message="logp not right")])
    fold_enrichment = FloatField('Fold Enrichment', validators=[DataRequired(message="logp not right")])

    # distance_from_TSS = IntegerField('Distance from TSS', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_upstream = IntegerField('Distance from TSS Upstream', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_downstream = IntegerField('Distance from TSS Downstream', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    
    email = StringField("Email", validators=[DataRequired(message="email not right")])

    submit = SubmitField('Submit')

    # New params


# Database Specs

db = SQLAlchemy(app)

meta_column_list = ["experiment_accession",
            "transcription_factors",
            "tissue_types"
            ]

class ChIP_Meta(db.Model):
    '''
    Stores metadata on ChIP-Seq Studies
    '''
    id = db.Column(db.Integer, primary_key=True)
    experiment_accession = db.Column(db.String(50), unique=True, nullable=False)
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return ('<Experiment Accession {}>'.format(self.experiment_accession)
        + '<Transcription Factor {}>'.format(self.transcription_factors)
        + '<Tissue Type {}>'.format(self.tissue_types) 
        )


# The order of this list determines write order 
peaks_column_list = [
            "chrom",
            "start",
            "end",
            "id",
            "length",
            "summit",
            "pileup",
            "log_p",
            "fold_enrichment",
            "log_q",
            "transcription_factors",
            "tissue_types",
            "experiment_accession"
            ]

class Peaks(db.Model):
    '''
    This stores all the peaks that have been collected, and uses the experiment accession # to relate them to each other
    '''

    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String, nullable=False)
    start = db.Column(db.Integer, nullable=False)
    end = db.Column(db.Integer, nullable=False)
    length = db.Column(db.Integer, nullable=False)
    summit = db.Column(db.Integer, nullable=False)
    pileup = db.Column(db.Float, nullable=False) 
    log_p = db.Column(db.Float, nullable=False)
    fold_enrichment = db.Column(db.Float, nullable=False)
    log_q = db.Column(db.Float, nullable=False)
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)
    experiment_accession = db.Column(db.String, nullable=False)

    def __repr__(self):
        return ('<Experiment Accession {}>'.format(self.experiment_accession)
        + '<Transcription Factor {}>'.format(self.transcription_factors)
        + '<Tissue Type {}>'.format(self.tissue_types) 
        + '<Chrom {}>'.format(self.chrom) 
        + '<Start {}>'.format(self.start) 
        + '<end {}>'.format(self.end) 
        )

class Footprints(db.Model):

    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String, nullable=False)
    start = db.Column(db.Integer, nullable=False)
    end = db.Column(db.Integer, nullable=False) 
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)

class Query_History(db.Model):
    '''
    Stores a history of queries and query data
    '''
    id = db.Column(db.Integer, primary_key=True)
    promoter = db.Column(db.Boolean)
    enrichment = db.Column(db.Boolean)

    transcription_factors = db.Column(db.String)
    tissue_types = db.Column(db.String)

    log_p = db.Column(db.Float)
    fold_enrichment = db.Column(db.Float)

    distance_from_TSS_upstream = db.Column(db.Integer)
    distance_from_TSS_downstream = db.Column(db.Integer)

    email = db.Column(db.String)

    date_posted = db.Column(db.DateTime, nullable = False, default = datetime.utcnow)


class Presets(db.Model):
    '''
    Store a set of preset TF matrices for rapid download
    '''
    id = db.Column(db.Integer, primary_key=True)
    log_p = db.Column(db.Float)
    fold_enrichment = db.Column(db.Float)
    promoter = db.Column(db.Boolean)
    enrichment = db.Column(db.Boolean)
    table_file_path = db.Column(db.String) 


# General Functions

def run_pipeline(user_params):
    # 1.) Query database for appropriate peaks
    # 2.) Toss them into a bed file and sort
    # 3.) Run HOMER annotation with bed file input
    # 4.) Build table using homer output
    # 5.) Return completion signal email output  

    # Step 1

    print("BEGIN PIPELINE *************************************************************************************************************")

    time_string = user_params["time"]

    peak_subset = (
                    Peaks.query.filter(Peaks.tissue_types.in_(user_params["tissue_types"]))
                    .filter(Peaks.transcription_factors.in_(user_params["transcription_factors"]))
                    )
            
    study_list = [x.experiment_accession for x in peak_subset]

    # for x in peak_subset:
    #     print(x)

    print("||||||||||||||||||passed peak subsetting")

    # End Step 1
    # Step 2

    temp_peaks_file = pwd + "/intermediates/" + "temppeaks_" + time_string + ".bed"
    peaks_list = convert_query_to_file(peaks_column_list, peak_subset, user_params,  temp_peaks_file)

    # os.system("sort -k 1,1 -k2,2n " + temp_peaks_file)

    print("||||||||||||||||||temp peaks bed file created and sorted")

    # End Step 2 
    # Step 3 

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only
    

        annotation_output = pwd + "/intermediates/" + "annotation_" + time_string + ".anno"
        peaks_command = [
            "annotatePeaks.pl",
            temp_peaks_file,
            "hg38",
            ">",
            annotation_output
        ]
        os.system(" ".join(peaks_command))

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only

        enhancer_bed = pwd + "/enhancer-gene/processed/unionAHSIL"
        intersect_output = pwd + "/intermediates/" + "annotation_" + time_string + ".bed"
        bed_command = [
            "bedtools",
            "intersect",
            "-a",
            enhancer_bed,
            "-b",
            temp_peaks_file,
            "-wb",
            ">",
            intersect_output
        ]
        os.system(" ".join(bed_command))

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer

        enhancer_bed = pwd + "/enhancer-gene/processed/unionAHSIL"
        intersect_output = pwd + "/intermediates/" + "annotation_" + time_string + ".bed"
        bed_command = [
            "bedtools",
            "intersect",
            "-a",
            enhancer_bed,
            "-b",
            temp_peaks_file,
            "-wb",
            ">",
            intersect_out
        ]

        annotation_output = pwd + "/intermediates/" + "annotation_" + time_string + ".anno"

        peaks_command = [
            "annotatePeaks.pl",
            temp_peaks_file,
            "hg38",
            ">",
            annotation_output
        ]

        os.system(" ".join(bed_command))
        os.system(" ".join(peaks_command))

    else:
        print("invalid query, neither promotor nor enhancer selected")
        return None


    print("||||||||||||||||||peaks annotated")

    # End Step 3
    # Step 4

    tg_table = create_empty_table(all_genes, all_tfs)

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only
        parse_promoter(user_params, annotation_output, tg_table)

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only
        parse_enhancer(user_params, intersect_output, tg_table)

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer
        parse_promoter(user_params, annotation_output, tg_table)
        parse_enhancer(user_params, intersect_output, tg_table)

    else:
        print("invalid query, neither promotor nor enhancer selected")
        return None

    tg_write_file = pwd + "/intermediates/" + "tgtable_" + time_string + ".tgtable"
    write_dict_tsv(tg_table, all_genes, all_tfs, tg_write_file)
    os.system("gzip " + tg_write_file)
    zipped_tg_table = tg_write_file + ".gz"
    solo_zipped_filename = zipped_tg_table.split("/")[-1]

    # End Step 4
    # Step 5

    send_from = "ChIPBaseApp@gmail.com"
    password = "chipbase"
    send_to = user_params["email"]
    subject = "Your output from ChIP-IO"
    text = "Here is your transcription factor - gene interaction table"
    server = 'smtp.gmail.com'
    send_mail(send_from, send_to, subject, text, zipped_tg_table, solo_zipped_filename, server, send_from, password)
    print(user_params["email"])

    # End Step 5
    # Clear up intermediate files
    

    print("END PIPELINE *************************************************************************************************************")

def write_dict_tsv(tg_table, all_genes, all_tfs, table_write):
    '''
    Writes dictionary table to a tsv file
    '''
    with open(table_write, "a") as tw:

        tw.write("Genes\t" + "\t".join(all_tfs) + "\n")

        for gene in all_genes:
            # if gene not in tg_table:
            #     continue
            cur_gene_string = gene + "\t"
            for tf in all_tfs:
                # if tf not in tg_table[gene]:
                #     continue
                cur_gene_string += (str(tg_table[gene][tf]) + "\t")
            tw.write(cur_gene_string.rstrip("\t")+"\n")

def create_empty_table(gene_list, tf_list):
    '''
    Creates an empty table (dictionary[tfs]->dictionary[genes]) initializing all values
    to 0 based upon a list of genes and tfs
    '''
    table = {}
    for gene in gene_list:
        tf_dict = {}
        for tf in tf_list:
            tf_dict[tf] = 0
        table[gene] = tf_dict
    
    return table

def constraints_met(data, user_params, constraints_type):
    '''
    Given some data in list form(a line in a bed file), the user input parameters, and the "type"
    of constraints corresponding to the given data, assesses whether the data "passes"
    '''
    if constraints_type == "peaks":
        # print(data)
        # print(user_params)
        if (float(data[7]) > user_params["log_p"] and
            float(data[8]) > user_params["fold_enrichment"]
        ):
            # print("Peaks constraints satisfied")
            return True
        else:
            return False
    elif constraints_type == "annotations":
        try:
            int(data[9])
        except:
            print("TSS DISTANCE NOT INTABLE: ", data[9])
            return False
        if ((int(data[9]) > (-1)*user_params["dist_tss_upstream"]) and (int(data[9]) < user_params["dist_tss_downstream"])
        ):
            print("TSS DISTANCE: " + data[9] + " < " + str(user_params["dist_tss"]))
            return True
        else:
            return False

def parse_enhancer(user_params, anno_file, tf_gene_table):
    '''
    Updates tf-gene table with enhancer annotation results which pass constraints
    '''
    line_count = 0
    with open(anno_file, "r") as anno:
        # print(tf_gene_table["PAX7"].keys())
        for line in anno:
            print(line)
            p_l= line.rstrip("\n").split("\t")
            # print(p_l, "annotation line")
            # print(line_count)
            if line_count == 0:
                line_count += 1
                continue


            peak_id = int(p_l[10])
            gene_id = p_l[6]
            ori_peak_tf = Peaks.query.filter(Peaks.id==peak_id)[0].transcription_factors
            print(ori_peak_tf, "ori peak")
            if gene_id in tf_gene_table:

                # print("CHECK gene")
                if ori_peak_tf in tf_gene_table[gene_id]:
                    # print("CHECK 2, added")
                    tf_gene_table[gene_id][ori_peak_tf] += 1
            else:
                print("Did not pass Check 1 ", gene_id)
            line_count += 1


def parse_promoter(user_params, anno_file, tf_gene_table):
    '''
    Updates tf-gene table with promoter annotation results which pass constraints
    '''
    line_count = 0
    with open(anno_file, "r") as anno:
        # print(tf_gene_table["PAX7"].keys())
        for line in anno:
            p_l= line.rstrip("\n").split("\t")
            # print(p_l, "annotation line")
            # print(line_count)
            if line_count == 0:
                line_count += 1
                continue

            if constraints_met(p_l, user_params, "annotations"):
                peak_id = int(p_l[0])
                ori_peak_tf = Peaks.query.filter(Peaks.id==peak_id)[0].transcription_factors
                print(ori_peak_tf, "ori peak")
                if p_l[15] in tf_gene_table:

                    print("CHECK 1")
                    if ori_peak_tf in tf_gene_table[p_l[15]]:
                        print("CHECK 2, added")
                        tf_gene_table[p_l[15]][ori_peak_tf] += 1
                else:
                    print("Did not pass Check 1 ", p_l[15])
            line_count += 1

def send_mail(send_from, send_to, subject, text, file_path, file_name, server, email_user, email_password):

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Subject'] = subject

    body = 'Hi there, sending this email from Python!'
    msg.attach(MIMEText(body,'plain'))

    filename=file_path
    attachment=open(filename,'rb')

    part = MIMEBase('application','octet-stream')
    part.set_payload((attachment).read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition',"attachment; filename= "+file_name)

    msg.attach(part)
    text = msg.as_string()
    server = smtplib.SMTP(server,587)
    server.starttls()
    server.login(email_user,email_password)


    server.sendmail(email_user,send_to,text)
    server.quit()

def convert_query_to_file(columns, query_result, user_params, file_path):
    '''
    Converts the output of a db query to a tsv file and saves the file
    '''
    peaks_list = []
    with open(file_path, "a") as f:
        x = 0
        for item in query_result:
            outputs = [str(getattr(item,x)) for x in columns]
            if constraints_met(outputs, user_params, "peaks"):
                peaks_list.append(outputs)
                writeable = "\t".join(outputs) + "\n"
                # print(writeable, "WRITE LINE")
                f.write(writeable)

    return peaks_list
    
def parse_input(in_string, delim, field, all_possible):

    if in_string.strip() == "ALL":
        return all_possible[field]

    return(in_string.strip().split(delim))

def build_query_hist(form):

    query_data = Query_History(
        log_p = form.log_p.data, 
        fold_enrichment = form.fold_enrichment.data, 
        promoter = False, 
        enrichment = True)

    return query_data

# View Functions

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/construct_query')
def construct_query():
    return render_template('construct_query.html')

@app.route('/promoter_form', methods=['GET', 'POST'])
def promoter_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL")
    print(form)
    if form.validate_on_submit():
        print("valid")
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": True, 
            "enhancer": False,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,

            "dist_tss": form.distance_from_TSS_upstream.data,
            "dist_tss": form.distance_from_TSS_upstream.data,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
    return render_template('promoter_form.html', form = form)

@app.route('/enhancer_form', methods=['GET', 'POST'])
def enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL", distance_from_TSS_upstream = 1, distance_from_TSS_downstream = 1) 
    if form.validate_on_submit():
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": False, 
            "enhancer": True,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,

            "dist_tss": 1,
            "dist_tss": 1,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
    return render_template('enhancer_form.html', form = form)

@app.route('/promoter_enhancer_form', methods=['GET', 'POST'])
def promoter_enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL")
    if form.validate_on_submit():
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": True, 
            "enhancer": True,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,

            "dist_tss": form.distance_from_TSS_upstream.data,
            "dist_tss": form.distance_from_TSS_upstream.data,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
    return render_template('promoter_enhancer_form.html', form = form)

@app.route('/output')
def output():
    return render_template('output.html', query_data = data)

@app.route('/complete')
def complete():
    return render_template('complete.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/downloads')
def downloads():
    # send_file("/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/peaks/ENCSR480LIS_peaks.xls", as_attachment=True)
    return render_template('downloads.html', download_files = download_files)

@app.route("/download_file/<file_path>", methods=['GET', 'POST'])
def download_file(file_path):
    print("downloading,", file_path)
    file_path_proper = "/".join(file_path.split("#"))
    if file_path is None:
        print("path is None")
    try:
        return send_file(file_path_proper, as_attachment=True)
    except Exception as e:
        print("problem with path")
    # return render_template('downloads.html', download_files = download_files)

@app.route('/contact')
def contact():
    return render_template('contact.html')


all_possible = {
"transcription_factors": list(set([x.transcription_factors for x in ChIP_Meta.query.all()])),
"tissue_types": list(set([x.tissue_types for x in ChIP_Meta.query.all()]))
}

download_files = DownloadFiles()
download_files.collect_peaks("peaks")
# print(download_files.peak_files)


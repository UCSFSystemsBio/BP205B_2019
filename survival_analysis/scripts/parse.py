import gzip
import csv
import io
import glob
import os
import pandas as pd
import numpy as np
files = glob.glob("./panc_expression/**/*.gz")

metalist = glob.glob("./panc_expression/**/*")
os.path.basename("./gdc_sample_sheet.2019-03-05.tsv")
name2tcga = "./gdc_sample_sheet.2019-03-05.tsv"
tcga = []
dirname=os.path.dirname
with open(name2tcga, 'r') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        for meta in metalist:
            if os.path.basename(meta) == row[1]:
                tcga.append([row[6], os.path.basename(dirname(meta))])

tcga2patientdata = "./paad_tcga_clinical_data.tsv"
tcga2patient = []
with open(tcga2patientdata, 'r') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    for row in spamreader:
        for tc in tcga:
            if tc[0][:-1] == row[2]:
                tmp = [tc[1], row[29], row[28], row[86]] #Disease Free (Months)	Disease \t Free Status \t sex
                tcga2patient.append(tmp)
df = pd.DataFrame(tcga2patient)
len(np.unique(df[2]))
df.to_csv("patient_metadata.tsv", sep='\t')

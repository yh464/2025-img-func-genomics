#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk

This script analyses if cortical size factor is specifically associated with a particular subtype of ADHD
'''

import numpy as np
import pandas as pd
import seaborn as sns
import os

plink = '/rds/project/rds-Nl99R8pHODQ/toolbox/plink'
tmpdir = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/temp'
outdir = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/test'
os.makedirs(tmpdir, exist_ok = True)
os.makedirs(outdir, exist_ok = True)

# BUHMBOX analysis
sumstats = f'{outdir}/buhmbox_input_ce.txt'
sumstats_temp = f'{tmpdir}/buhmbox_input_ce.txt'
if not os.path.exists(sumstats):
    ce_clump = pd.read_table('/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/clump/structural_factors/cortical_expansion_5e-08.clumped', usecols = ['SNP'])
    ce_gwas = pd.read_table('/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gwa/structural_factors/cortical_expansion.fastGWA', usecols = ['SNP','A1','AF1','BETA'])
    ce_clumped = pd.merge(ce_clump, ce_gwas, on = 'SNP', how = 'inner')
    ce_clumped['OR'] = np.exp(ce_clumped['BETA']) # BUHMBOX only accepts OR as input
    ce_clumped[['SNP','A1','AF1','OR']].to_csv(sumstats, sep = '\t', index = False)
    ce_clumped[['SNP','A1']].to_csv(sumstats_temp, sep = '\t', index = False)

cases_bed = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/params/bed_spark/alpha_omega_adhd'
cases_recode = f'{outdir}/adhd_recode.txt'
if not os.path.exists(cases_recode):
    os.system(f'{plink} --bfile {cases_bed} --extract {sumstats_temp} --reference-allele {sumstats_temp} --recode A --out {cases_recode[:-4]}')

ctrls_bed = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/params/bed_spark/alpha_omega_noadhd'
ctrls_recode = f'{outdir}/noadhd_recode.txt'
if not os.path.exists(ctrls_recode): 
    os.system(f'{plink} --bfile {ctrls_bed} --extract {sumstats_temp} --reference-allele {sumstats_temp} --recode A --out {ctrls_recode[:-4]}')

pc_file = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/params/bed_spark/alpha_omega.eigenvec'

buhmbox_script = f'{tmpdir}/bhumbox.R'
if not os.path.exists(buhmbox_script):
    os.system(f'wget https://software.broadinstitute.org/mpg/buhmbox/data/buhmbox_v0.38/buhmbox.R -O {buhmbox_script}')

buhmbox_out = f'{outdir}/buhmbox_ce_adhd'
os.system(f'Rscript {buhmbox_script} {sumstats} {cases_recode}.raw {ctrls_recode}.raw YY N Y {buhmbox_out} {pc_file}')
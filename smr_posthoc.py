import pandas as pd
import numpy as np

gene_sets = pd.read_table('./braun_2023_s4.txt', index_col = 0)
ref = pd.read_table('/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt', index_col = 'LABEL')
neuronal_like_ipc_genes = [gene for gene in gene_sets.index[:138] if gene in ref.index]
rg_like_ipc_genes = [gene for gene in gene_sets.index[138:] if gene in ref.index]
neuronal_like_ipc_genes = ref.loc[neuronal_like_ipc_genes, 'GENE'].values
rg_like_ipc_genes = ref.loc[rg_like_ipc_genes, 'GENE'].values

sig_tf = ['MEF2C','NFIA','NFIB','TCF4','TCF12','ZBTB20','LHX2','HMGA2','PAX6']
sig_tf = ref.loc[sig_tf, 'GENE'].values

cnmf_k = 16; cnmf_dt = '0_15'
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]

smr_files = [
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/structural_factors/cortical_expansion.smr',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/disorders/adhd2022.smr',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/annot/smr/disorders/adhd2025.smr'
]

smr_df = [ pd.read_table(f) for f in smr_files ]

def extract_smr(df, gene_set):
    gene_set = pd.Index(gene_set)
    genes = gene_set.intersection(df['probe'])
    out = df.loc[df['probe'].isin(genes),:]
    print(out)
    return out

def extract_sig(df, p_threshold = 0.05):
    df = df.loc[df.qtl == 'psychencode_eqtl',:]
    df = df.loc[(df.q < p_threshold) & (df.p_heidi > .01),:]
    return df

for smr_file, df in zip(smr_files, smr_df):
    print(smr_file)
    print('Neuronal-like IPC genes:')
    extract_smr(df, neuronal_like_ipc_genes)
    print('RG-like IPC genes:')
    extract_smr(df, rg_like_ipc_genes)
    for programme in ['F9','F2','F1','F8']:
        print(f'Top 200 of {programme}:')
        top_200_genes = cnmf_weights.nlargest(200, programme).index
        extract_smr(df, top_200_genes)
    print('Significant transcription factors:')
    extract_smr(df, sig_tf)

smr_sig = [ extract_sig(df, 0.1) for df in smr_df ]
genes_sig = np.union1d(
    np.intersect1d(smr_sig[0]['gene'].values, smr_sig[1]['gene'].values),
    np.intersect1d(smr_sig[0]['gene'].values, smr_sig[2]['gene'].values)
)

smr_overlap = pd.concat([ df.loc[(df['gene'].isin(genes_sig)) & (df.qtl == 'psychencode_eqtl'),['beta', 'gene']
    ].set_index('gene').rename(columns = {'beta': df.iloc[0,0]}) for df in smr_df ], axis = 1)
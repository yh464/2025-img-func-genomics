# Gene set enrichment tests for the cNMF programmes
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control as fdr
rng = np.random.default_rng(19260817)

cnmf_k = 16; cnmf_dt = '0_15'
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]
replication_k = 10; replication_dt = '0_1'
replication_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/polioudakis_2019/neocx_wang_2025_genes/'+
    f'neocx_wang_2025_genes.spectra.k_{replication_k}.dt_{replication_dt}.consensus.txt',
    index_col = 0
).T
replication_weights.columns = [f'F{i+1}' for i in range(replication_k)]

gene_sets = pd.read_table('./braun_2023_s4.txt', index_col = 0)
ref = pd.read_table('/home/yh464/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/genes_ref.txt', index_col = 'LABEL')
neuronal_like_ipc_genes = [gene for gene in gene_sets.index[:138] if gene in ref.index]
rg_like_ipc_genes = [gene for gene in gene_sets.index[138:] if gene in ref.index]
neuronal_like_ipc_genes = ref.loc[neuronal_like_ipc_genes, 'GENE'].values
rg_like_ipc_genes = ref.loc[rg_like_ipc_genes, 'GENE'].values

def gene_set_enrichment(weights, gene_set, n_perm = 100000):
    n_overlap = weights.index.intersection(gene_set).size
    mean_weight = weights.loc[weights.index.intersection(gene_set),:].mean(axis = 0).values
    perm_weights = np.stack([
        weights.sample(n = n_overlap, replace = False, random_state = rng).mean(axis = 0).values
        for _ in range(n_perm)
    ], axis = 0)
    p_values = (mean_weight > perm_weights).mean(axis = 0)
    p_values = np.minimum(p_values, 1 - p_values) * 2  # two-sided
    out = pd.DataFrame({'mean_weight': mean_weight, 'perm_mean': perm_weights.mean(axis = 0), 'p_value': p_values}, index = weights.columns)
    print(out)
    return out

cnmf_gsea_neuronal_like_ipc = gene_set_enrichment(cnmf_weights, neuronal_like_ipc_genes)
cnmf_gsea_rg_like_ipc = gene_set_enrichment(cnmf_weights, rg_like_ipc_genes)
replication_gsea_neuronal_like_ipc = gene_set_enrichment(replication_weights, neuronal_like_ipc_genes)
replication_gsea_rg_like_ipc = gene_set_enrichment(replication_weights, rg_like_ipc_genes)

gnova = pd.read_table(
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gcorr/gnova/disorders.structural_factors/'+
    'disorders_adhd2022.structural_factors_cortical_expansion.wang_2025.eregulons.gnova.txt'
)
gnova['eregulon'] = gnova['annot_name'].str.replace('eregulons.','').str.replace('_-','').str.replace('_+','')
gnova['gene'] = [ref.loc[x, 'GENE'] if x in ref.index else np.nan for x in gnova['eregulon']]
gnova_max = pd.concat([
    gnova.loc[gnova['rho_corrected']<0, ['gene','rho_corrected']].groupby('gene').min(),
    gnova.loc[gnova['rho_corrected']>0, ['gene','rho_corrected']].groupby('gene').max()
])
n_overlap = cnmf_weights.index.intersection(gnova_max.index).size
gnova_corr = cnmf_weights.corrwith(gnova_max['rho_corrected'], axis = 0)
gnova_corr

gnova_spearman = cnmf_weights.corrwith(gnova_max['rho_corrected'], axis = 0, method = 'spearman')
gnova_spearman

tf_adhd = pd.read_table(
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/magma_gsea/disorders/adhd2022/adhd2022.ENSG_10kb.wang_2025_eregulons.gsa.out', 
    index_col = 0, comment = '#', sep = '\\s+')
tf_adhd['q'] = fdr(tf_adhd['P'].values)
tf_expansion = pd.read_table(
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/magma_gsea/structural_factors/cortical_expansion/cortical_expansion.ENSG_10kb.wang_2025_eregulons.gsa.out', 
    index_col = 0, comment = '#', sep = '\\s+'
)
tf_expansion['q'] = fdr(tf_expansion['P'].values)
overlap_q005 = tf_adhd.loc[tf_adhd['q'] < 0.05, :].index.intersection(
    tf_expansion.loc[tf_expansion['q'] < 0.05, :].index
)
overlap_q005 = [x for x in overlap_q005.str.replace('eregulons.','').str.replace('_-','').str.replace('_+','') if x in ref.index]
overlap_q005 = ref.loc[overlap_q005, 'GENE'].values
print(overlap_q005)
overlap_005_enrichment = gene_set_enrichment(cnmf_weights, overlap_q005)

overlap_q01 = tf_adhd.loc[tf_adhd['q'] < 0.1, :].index.intersection(
    tf_expansion.loc[tf_expansion['q'] < 0.1, :].index
)
overlap_q01 = [x for x in overlap_q01.str.replace('eregulons.','').str.replace('_-','').str.replace('_+','') if x in ref.index]
overlap_q01 = ref.loc[overlap_q01, 'GENE'].values
print(overlap_q01)
overlap_01_enrichment = gene_set_enrichment(cnmf_weights, overlap_q01)
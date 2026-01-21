# Gene set enrichment tests for the cNMF programmes
import pandas as pd

cnmf_k = 16; cnmf_dt = '0_15'
cnmf_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_neocortex/' +
    f'wang_2025_neocortex.spectra.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',index_col = 0).T
cnmf_weights.columns = [f'F{i+1}' for i in range(cnmf_k)]
replication_k = 14; replication_dt = '0_1'
replication_weights = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/polioudakis_2019/neocx_wang_2025_genes/'+
    f'neocx_wang_2025_genes.spectra.k_{replication_k}.dt_{replication_dt}.consensus.txt',
    index_col = 0
).T
replication_weights.columns = [f'F{i+1}' for i in range(replication_k)]

import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os

cnmf_k = 16; cnmf_dt = '0_15'; factor_order = [f'F{i}' for i in [12,15,5,9,14,2,1,8,4,6,7,11,13,3,10,16]]
cnmf_scores = pd.read_table(
    '/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/'+
    f'wang_2025_neocortex/wang_2025_neocortex.usages.k_{cnmf_k}.dt_{cnmf_dt}.consensus.txt',
    index_col = 0)
cnmf_scores.columns = [f'F{i+1}' for i in range(cnmf_k)]
cnmf_scores = cnmf_scores[factor_order]

scdrs_scores = [
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/scdrs/disorders/adhd2022/wang_2025/adhd2022.wang_2025_neocortex.scdrs.score.txt',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/scdrs/disorders/adhd2025/wang_2025/adhd2025.wang_2025_neocortex.scdrs.score.txt',
    '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/sc/scdrs/structural_factors/cortical_expansion/wang_2025/cortical_expansion.wang_2025_neocortex.scdrs.score.txt'
]

adata = sc.read_h5ad('/rds/project/rds-Nl99R8pHODQ/multiomics/scdrs/wang_2025/wang_2025_neocortex.h5ad', 'r')
cell_types = adata.obs['Type_updated']
cell_type_order = [
    'RG-vRG','RG-tRG','RG-oRG',
    'IPC-EN','EN-Newborn','EN-IT-Immature','EN-L2_3-IT','EN-L4-IT','EN-L5-IT','EN-L6-IT','EN-Non-IT-Immature','EN-L5-ET','EN-L5_6-NP','EN-L6-CT','EN-L6b',
    'IN-NCx_dGE-Immature','IN-CGE-Immature','IN-CGE-SNCG','IN-CGE-VIP','IN-Mix-LAMP5','IN-MGE-Immature','IN-MGE-PV','IN-MGE-SST',
    'Tri-IPC','Astrocyte-Immature','Astrocyte-Protoplasmic','Astrocyte-Fibrous','OPC','Oligodendrocyte-Immature','Oligodendrocyte',
    'Cajal-Retzius cell','Microglia','Vascular','Unknown'
]
output_dir = '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/test'

def test_scdrs_corr(input_file):
    scdrs_score = pd.read_table(input_file, index_col = 0, usecols = [0,2])
    out = []
    corr = cnmf_scores.corrwith(scdrs_score['norm_score']).to_frame(name = 'All')
    out.append(corr)
    for ct in cell_types.unique():
        cells = cell_types[cell_types == ct].index
        corr = cnmf_scores.loc[cells, :].corrwith(scdrs_score.loc[cells, 'norm_score']).to_frame(name = ct)
        out.append(corr)
    result = pd.concat(out, axis = 1)
    result = result.loc[factor_order, ['All'] + cell_type_order]
    result.index = result.index.str.replace('F', '')
    out_prefix = f'{output_dir}/{os.path.basename(os.path.dirname(os.path.dirname(input_file)))}_cnmf_k{cnmf_k}_corr'

    fig = plt.figure(figsize = (12, 6))
    sns.heatmap(result, cmap = 'vlag', center = 0, yticklabels = True, xticklabels = True, square = True)
    fig.savefig(f'{out_prefix}.pdf', bbox_inches = 'tight')
    result.to_csv(f'{out_prefix}.txt', sep = '\t')
    plt.close()

for file in scdrs_scores:
    print(file)
    test_scdrs_corr(file)
import os
import pandas as pd
import statsmodels.formula.api as smf
diff_efficiency = pd.read_table('./jerber_2021_s5.txt', index_col = 0).dropna()
diff_efficiency = diff_efficiency.loc[diff_efficiency['in_study'] == 'succeeded', :]
metadata = pd.read_table('./jerber_2021_s10.txt', index_col = 0)
metadata.index.name = None; diff_efficiency.index.name = None

hipsci_pca = pd.read_table('/rds/project/rds-Nl99R8pHODQ/hipsci/merged/hipsci.merged.eigenvec', sep = '\\s+', header = None)
hipsci_pca.columns = ['FID','IID'] + [f'PC{i+1}' for i in range(hipsci_pca.shape[1]-2)]
hipsci_pca.index = hipsci_pca['FID'] + '_' + hipsci_pca['IID']
hipsci_pca = hipsci_pca.drop(columns = ['FID','IID'])

adhd2022_prs = pd.read_table('/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/prs/prs_score_hipsci/disorders/adhd2022.txt')
adhd2022_prs.index = adhd2022_prs['FID'] + '_' + adhd2022_prs['IID']
diff_efficiency['diff_efficiency'].corr(adhd2022_prs['score_norm'])

cortical_expansion_prs = pd.read_table('/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/prs/prs_score_hipsci/structural_factors/cortical_expansion.txt')
cortical_expansion_prs.index = cortical_expansion_prs['FID'] + '_' + cortical_expansion_prs['IID']
diff_efficiency['diff_efficiency'].corr(cortical_expansion_prs['score_norm'])

def test_prs(prs, test_var = None):
    if test_var is None: test_var = diff_efficiency['diff_efficiency']
    df = pd.concat([test_var, metadata[['Donor_Sex','Avg_Donor_Age']], hipsci_pca.iloc[:, :20], prs[['FID','score_norm']]], axis = 1).dropna()
    # model = smf.mixedlm('diff_efficiency ~ score_norm + ' + ' + '.join(hipsci_pca.columns[:20]), data = df, groups = df['FID']).fit()
    # print(model.summary())
    model = smf.ols(f'{test_var.name} ~ score_norm + Donor_Sex + Avg_Donor_Age + ' + ' + '.join(hipsci_pca.columns[:20]), data = df).fit()
    if model.pvalues['score_norm'] < 0.05: print(model.summary())
    # model_individual_only = smf.mixedlm('diff_efficiency ~ score_norm', data = df, groups = df['FID']).fit()
    # print(model_individual_only.summary())
    model_4pc = smf.ols(f'{test_var.name} ~ score_norm + Donor_Sex + Avg_Donor_Age + ' + ' + '.join(hipsci_pca.columns[:4]), data = df).fit()
    if model_4pc.pvalues['score_norm'] < 0.05: print(model_4pc.summary())
    model_0pc = smf.ols(f'{test_var.name} ~ score_norm + Donor_Sex + Avg_Donor_Age' , data = df).fit()
    if model_0pc.pvalues['score_norm'] < 0.05: print(model_0pc.summary())
    return model

test_prs(adhd2022_prs)
test_prs(cortical_expansion_prs)

# model with cNMF programmes
from sklearn.decomposition import NMF
import scanpy as sc
cnmf_programmes = pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/programmes/cnmf/wang_2025/wang_2025_magma_genes/wang_2025_magma_genes.spectra.k_16.dt_0_1.consensus.txt', index_col = 0)
for d in [11, 30, 52]:
    adata = f'/rds/project/rds-Nl99R8pHODQ/multiomics/raw/jerber_2021/jerber_2021_d{d}.h5ad'
    adata = sc.read_h5ad(adata,'r')
    out_file = f'/rds/project/rds-Nl99R8pHODQ/multiomics/raw/jerber_2021/d{d}_wang_2025_k16_activity.txt'
    if not os.path.isfile(out_file):
        overlap = cnmf_programmes.columns.intersection(adata.var_names)
        model = NMF(n_components = 16)
        model.components_ = cnmf_programmes[overlap].values
        adata = adata[:, overlap].to_memory()
        adata.X = adata.X.astype(float)
        sc.pp.normalize_total(adata, target_sum = 1e6)
        sc.pp.scale(adata, zero_center = False)
        activity = model.transform(adata.X)
        activity = pd.DataFrame(activity, index = adata.obs_names, columns = [f'F{i}' for i in range(1, 17)])
        activity.to_csv(out_file, sep = '\t')
    else: activity = pd.read_table(out_file, index_col = 0)
    for factor in [7, 3, 2, 8]:
        tmp = pd.concat([adata.obs, activity[[f'F{factor}']]], axis = 1).dropna()
        tmp = tmp[['donor_id',f'F{factor}']].groupby('donor_id').mean()
        test_prs(cortical_expansion_prs, test_var = tmp[f'F{factor}'])
        test_prs(adhd2022_prs, test_var = tmp[f'F{factor}'])
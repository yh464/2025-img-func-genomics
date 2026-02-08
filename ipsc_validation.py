import pandas as pd
import statsmodels.formula.api as smf
diff_efficiency = pd.read_table('./jerber_2021_s5.txt', index_col = 0).dropna()
diff_efficiency = diff_efficiency.loc[diff_efficiency['in_study'] == 'succeeded', :]
metadata = pd.read_table('./jerber_2021_s10.txt', index_col = 0)
metadata.index.name = None; diff_efficiency.index.name = None
diff_efficiency = pd.concat([diff_efficiency, metadata.loc[diff_efficiency.index, :]], axis = 1)

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

def test_prs(prs):
    df = pd.concat([diff_efficiency[['diff_efficiency', 'Donor_Sex','Avg_Donor_age']], hipsci_pca.iloc[:, :20], prs[['FID','score_norm']]], axis = 1).dropna()
    # model = smf.mixedlm('diff_efficiency ~ score_norm + ' + ' + '.join(hipsci_pca.columns[:20]), data = df, groups = df['FID']).fit()
    # print(model.summary())
    model = smf.ols('diff_efficiency ~ score_norm + Donor_Sex + Avg_Donor_age + ' + ' + '.join(hipsci_pca.columns[:20]), data = df).fit()
    print(model.summary())
    # model_individual_only = smf.mixedlm('diff_efficiency ~ score_norm', data = df, groups = df['FID']).fit()
    # print(model_individual_only.summary())
    model_4pc = smf.ols('diff_efficiency ~ score_norm + Donor_Sex + Avg_Donor_age + ' + ' + '.join(hipsci_pca.columns[:4]), data = df).fit()
    print(model_4pc.summary())
    model_0pc = smf.ols('diff_efficiency ~ score_norm + Donor_Sex + Avg_Donor_age' , data = df).fit()
    print(model_0pc.summary())
    return model

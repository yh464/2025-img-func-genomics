import pandas as pd
diff_efficiency = pd.read_table('./jerber_2021_s5.txt', index_col = 0).dropna()

adhd2022_prs = pd.read_table('/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/prs/prs_score_hipsci/disorders/adhd2022.txt')
adhd2022_prs.index = adhd2022_prs['FID'] + '_' + adhd2022_prs['IID']
diff_efficiency['diff_efficiency'].corrwith(adhd2022_prs['score_norm'])
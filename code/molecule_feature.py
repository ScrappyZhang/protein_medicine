import numpy as np
import pandas as pd

df_molecule = pd.read_csv('../input/df_molecule.csv') # 小分子信息，存在缺失值
# 指纹信息处理
molecule_fingerprint = df_molecule['Fingerprint'].str.split(',')
molecule_fingerprint = molecule_fingerprint.tolist()
molecule_fingerprint = pd.DataFrame(molecule_fingerprint)
molecule_fingerprint.columns = ['finger_{}'.format(i) for i in range(1, 168)]
molecule_fingerprint = molecule_fingerprint.astype(int)

df_molecule = df_molecule.drop(['Fingerprint'], axis=1)
# 缺失值处理
df_molecule_feat = df_molecule.drop(['Molecule_ID'], axis=1)

# df_molecule_feat= df_molecule_feat.fillna(0)
# df_molecule_feat= df_molecule_feat.fillna(df_molecule_feat.mean())

df_molecule_feat['cyp_3a4'] = df_molecule_feat['cyp_3a4'].fillna(df_molecule_feat['cyp_3a4'].mean())
df_molecule_feat['hia'] = df_molecule_feat['hia'].fillna(df_molecule_feat['hia'].mode()[0])
df_molecule_feat['ames_toxicity'] = df_molecule_feat['ames_toxicity'].fillna(0.87)
df_molecule_feat['tetrahymena_pyriformis_toxicity'] = df_molecule_feat['tetrahymena_pyriformis_toxicity'].fillna(df_molecule_feat['tetrahymena_pyriformis_toxicity'].mode()[0])
df_molecule_feat['honey_bee'] = df_molecule_feat['honey_bee'].fillna(df_molecule_feat['honey_bee'].quantile(0.75))
df_molecule_feat['logP'] = df_molecule_feat['logP'].fillna(df_molecule_feat['logP'].mean())
df_molecule_feat['CLtotal'] = df_molecule_feat['CLtotal'].fillna(df_molecule_feat['CLtotal'].mode()[0])
df_molecule_feat['NOAEL'] = df_molecule_feat['NOAEL'].fillna(df_molecule_feat['NOAEL'].mean())
df_molecule_feat['p_glycoprotein_inhibition'] = df_molecule_feat['p_glycoprotein_inhibition'].fillna(df_molecule_feat['p_glycoprotein_inhibition'].mean())
df_molecule_feat['fathead_minnow_toxicity'] = df_molecule_feat['fathead_minnow_toxicity'].fillna(df_molecule_feat['fathead_minnow_toxicity'].mode()[0])
df_molecule_feat['bbb'] = df_molecule_feat['bbb'].fillna(df_molecule_feat['bbb'].mode()[0])
df_molecule_feat['Vdd'] = df_molecule_feat['Vdd'].fillna(df_molecule_feat['Vdd'].mode()[0])
df_molecule_feat['biodegradation'] = df_molecule_feat['biodegradation'].fillna(df_molecule_feat['biodegradation'].mode()[0])
df_molecule_feat['cyp_2d6'] = df_molecule_feat['cyp_2d6'].fillna(df_molecule_feat['cyp_2d6'].quantile(0.75))

# 合并理化性质和指纹信息
df_molecule_feat['Molecule_ID'] = df_molecule['Molecule_ID']
df_molecu = pd.concat([df_molecule, molecule_fingerprint],axis=1)
df_molecu.to_csv('../output/molecule_feature.csv', index=False)
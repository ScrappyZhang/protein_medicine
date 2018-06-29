import os
import pandas as pd

data_path = '../input/'
os.chdir(data_path)
print(os.getcwd())

# 1.读取数据
df_protein_train = pd.read_csv('df_protein_train.csv')
df_protein_test = pd.read_csv('df_protein_test.csv')

protein_concat = pd.concat([df_protein_train, df_protein_test])
# print(protein_concat.iloc[336])
print(protein_concat.head())


# print(protein_concat['Sequence'])

def getP(P):
    AA = "ACDEFGHIKLMNPQRSTVWY"
    number = {k: P.count(k) for k in AA}
    P_length = len(P)
    XA = number['A'] / P_length
    XV = number['V'] / P_length
    XI = number['I'] / P_length
    XL = number['L'] / P_length
 
    Ai = 100 * (XA + 2.9 * XV + 3.9 * (XI + XL))
    return round(Ai,8)

# 构建蛋白质的脂肪系数
charc = pd.DataFrame([getP(i.upper()) for i in protein_concat['Sequence']])
charc.columns = ['Ai']
charc['Ai'] = charc.Ai.astype(float)
charc['Protein_ID'] = pd.Series(protein_concat['Protein_ID'].values)
print(charc.head())

charc.to_csv('../output/protein_get_Ai.csv', sep=',', header=True, index=False)

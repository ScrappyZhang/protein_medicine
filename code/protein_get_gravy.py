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
    GRAVY= {'A': 1.800, 'C': 2.50, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2,
                  'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6,
                  'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}

    P_length = len(P)
    number = {k: P.count(k) for k in AA}
    a = 0
    for key, value in number.items():
        a += GRAVY[key] * value
    a = a / P_length

    return round(a,8)

# 构建蛋白质的相对分子质量、Pi值 、 晓光系数
charc = pd.DataFrame([getP(i.upper()) for i in protein_concat['Sequence']])
# print charc
charc.columns = ['gravy']
# charc.columns = ['pi', 'xx']
charc['gravy'] = charc.gravy.astype(float)
charc['Protein_ID'] = pd.Series(protein_concat['Protein_ID'].values)
print(charc.head())
# charc.to_csv('protein_getp.csv', sep=',', header=True, index=False)
# charc.to_csv('protein_getp_no_drop.csv', sep=',', header=True, index=False)
charc.to_csv('../output/protein_get_gravy.csv', sep=',', header=True, index=False)

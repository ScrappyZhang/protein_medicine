import os
import pandas as pd

data_path = '../input/'
os.chdir(data_path)
print(os.getcwd())

# 1.读取数据
df_protein_train = pd.read_csv('df_protein_train.csv')
df_protein_test = pd.read_csv('df_protein_test.csv')

protein_concat = pd.concat([df_protein_train, df_protein_test])
# print(protein_concat.head())


def getP(P):
    AA = "ACDEFGHIKLMNPQRSTVWY"
    AA_residue = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766, 'G': 57.0519, 'H': 137.1411,
                  'I': 113.1594, 'K': 128.1741, 'L': 113.1594, 'M': 131.1926, 'N': 114.1038, 'P': 97.1167,
                  'Q': 128.1307, 'R': 156.1875, 'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.176}
    # AA_residue = {'A': 71.037114, 'C': 103.00919, 'D': 115.02694, 'E': 129.4259, 'F': 147.06841, 'G': 57.021464, 'H': 137.05894,
    #               'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04048, 'N': 114.04293, 'P': 97.052764,
    #               'Q': 128.05858, 'R': 156.10111, 'S': 87.032029, 'T': 101.04768, 'V': 99.068414, 'W': 186.07931, 'Y': 163.176}
    pI_e = {'C': 9.0, 'D': 4.0, 'E': 4.5, 'H': 6.4, 'K': 10.4, 'R': 12.0, 'Y': 10.0}
    COOH = "CDEY" # 带负电残基 
    NH2 = "HKR"  # 带正电残基
    number = {k: P.count(k) for k in AA}
    a = 0
    # 相对分子质量
    for k in AA:
        a = a + number[k] * AA_residue[k]
    # 晓光系数
    E = (number['Y'] * 1490 + number['W'] * 5500 + number['C'] / 2 * 125) / a
    # Pi迭代计算
    def f(x):
        b = 0
        c = 0
        for m in COOH:
            # 负电计算
            b += (number[m] * (10 ** x)) / (10 ** x + 10 ** pI_e[m])
        for n in NH2:
            c += (number[n] * 10 ** pI_e[n]) / (10 ** x + 10 ** pI_e[n])
        return b + 10 ** x / (10 ** x + 10 ** 3.2) - 10 ** 8.2 / (10 ** x + 10 ** 8.2) - c  # 3.2为末端羧基， 8.2为末端氨基

    r = 3.2 # 负电中最小
    s = 12.0 # 正电中最大
    x = (r + s) / 2
    for i in range(11):
        if f(x) > 0:
            s = x
            x = (r + s) / 2
        elif f(x) < 0:
            r = x
            x = (r + s) / 2
    return round((a+18.01524)/1000,8), round(x, 8), round(E, 8)


# 构建蛋白质的相对分子质量、Pi值 、 晓光系数
charc = pd.DataFrame([getP(i.upper()) for i in protein_concat['Sequence']])
# print charc
charc.columns = ['xdfz', 'pi', 'xx']
# charc.columns = ['pi', 'xx']
charc['xdfz'] = charc.xdfz.astype(float)
charc['pi'] = charc.pi.astype(float)
charc['xx'] = charc.xx.astype(float)
charc['Protein_ID'] = pd.Series(protein_concat['Protein_ID'].values)
print(charc.head())
# charc.to_csv('protein_getp.csv', sep=',', header=True, index=False)
# charc.to_csv('protein_getp_no_drop.csv', sep=',', header=True, index=False)
charc.to_csv('../output/protein_getp.csv', sep=',', header=True, index=False)

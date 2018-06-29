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

def getC(P):
    # AA_name = "APGCSNQDERKHYWFMLIVT"
    AA_name = "APSTCIMVLFWYHQRKNEDG"
    AA_position = {AA_name[0]: (1.00000000e+00, 0, 1),
                   AA_name[1]: (9.51056516e-01, 3.09016994e-01, 2),
                   AA_name[2]: (8.09016994e-01, 5.87785252e-01, 3),
                   AA_name[3]: (5.87785252e-01, 8.09016994e-01, 4),
                   AA_name[4]: (3.09016994e-01, 9.51056516e-01, 5),
                   AA_name[5]: (6.12323400e-17, 1.00000000e+00, 6),
                   AA_name[6]: (-3.09016994e-01, 9.51056516e-01, 7),
                   AA_name[7]: (-5.87785252e-01, 8.09016994e-01, 8),
                   AA_name[8]: (-8.09016994e-01, 5.87785252e-01, 9),
                   AA_name[9]: (-9.51056516e-01, 3.09016994e-01, 10),
                   AA_name[10]: (-1.00000000e+00, 1.22464680e-16, 11),
                   AA_name[11]: (-9.51056516e-01, -3.09016994e-01, 12),
                   AA_name[12]: (-8.09016994e-01, -5.87785252e-01, 13),
                   AA_name[13]: (-5.87785252e-01, -8.09016994e-01, 14),
                   AA_name[14]: (-3.09016994e-01, -9.51056516e-01, 15),
                   AA_name[15]: (-1.83697020e-16, -1.00000000e+00, 16),
                   AA_name[16]: (3.09016994e-01, -9.51056516e-01, 17),
                   AA_name[17]: (5.87785252e-01, -8.09016994e-01, 18),
                   AA_name[18]: (8.09016994e-01, -5.87785252e-01, 19),
                   AA_name[19]: (9.51056516e-01, -3.09016994e-01, 20)
                   }
    # x, y, z
    positon_x = 0
    positon_y = 0
    positon_z = 0
    P_length = len(P)
    for i, AA in enumerate(P):
        positon_x = positon_x + AA_position[AA][0] * (0.5 * (1 - 0.5 ** (P_length - i)) / 0.5)
        positon_y = positon_y + AA_position[AA][1] * (0.5 * (1 - 0.5 ** (P_length - i)) / 0.5)
        # positon_z = positon_z + i + 1
        # positon_z = positon_z + (i + 1) * (0.5 * (1 - 0.5 ** (P_length - i)) / 0.5)
        positon_z = positon_z + AA_position[AA][2] * (0.5 * (1 - 0.5 ** (P_length - i)) / 0.5)
        # positon_z = positon_z + AA_position[AA][2]
    E_x = positon_x / P_length
    E_y = positon_y / P_length
    # E_z = positon_z / P_length / P_length
    E_z = positon_z / P_length
    # E_z = positon_z

    return round(E_x, 8), round(E_y, 8), round(E_z, 8)



# 构建一阶中心距
charc = pd.DataFrame([getC(i.upper()) for i in protein_concat['Sequence']])
# print charc
charc.columns = ['E_1x', 'E_1y', 'E_z']
print(charc.head())
charc['E_1x'] = charc.E_1x.astype(float)
charc['E_1y'] = charc.E_1y.astype(float)
charc['E_z'] = charc.E_z.astype(float)
charc['Protein_ID'] = pd.Series(protein_concat['Protein_ID'].values)

charc.to_csv('../output/protein_getC.csv', sep=',', header=True, index=False)

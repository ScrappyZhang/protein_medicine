import re
import sys
import time
import os
import datetime

import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.metrics import mean_squared_error



# 1.读取数据
# 1.1 读取亲和力数据

data_path = '../input/'
os.chdir(data_path)
print(os.getcwd())


# data = pd.read_csv('data_not_nan_not_protein.csv')
# print(data.head())

# # 分子指纹
# feat = pd.read_csv('data_fingerprint.csv')
# print(feat.head())
# data = data.merge(feat, on='Molecule_ID', how='left')

df_affinity_train = pd.read_csv('df_affinity_train.csv')
df_affinity_test = pd.read_csv('df_affinity_test_toBePredicted.csv')
df_affinity_test['Ki'] = -11
data = pd.concat([df_affinity_train, df_affinity_test])

data_path = '../output/'
os.chdir(data_path)
print(os.getcwd())

# 1.1 药物数据
df_molecule = pd.read_csv('molecule_feature.csv')
data = data.merge(df_molecule, on='Molecule_ID', how='left')


# 1.2 蛋白质的特征


# 相对分子质量、pi、晓光系数
data_feat = pd.read_csv('protein_getp.csv')
print(data_feat.head())


feat = pd.read_csv('../../output/protein_get_II.csv') # 不稳定系数
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

feat = pd.read_csv('../../output/protein_get_Ai.csv') # 脂肪系数
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

feat = pd.read_csv('../../output/protein_get_gravy.csv') # GRAVY系数
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

#
# # 基于氨基酸理化性质的9大性质的APSTCIMVLFWYHQRKNEDG计算一阶中心距
feat = pd.read_csv('protein_getC.csv')
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

# 氨基酸主成分AAC 20维向量
# feat = pd.read_csv('protein_getAAC.csv')
# print(feat.head())
# data_feat = data_feat.merge(feat, on='Protein_ID', how='left')

data = data.merge(data_feat, on='Protein_ID', how='left')
del data_feat
print(data.shape)


# 3. lgb
# 3.1 缺失值0
train_feat = data[data['Ki'] > -11] # .fillna(0)
testt_feat = data[data['Ki'] <= -11] # .fillna(0)

# 线下验证集
offline_choice_list = [1, 2]
idx = train_feat['Protein_ID'].isin(offline_choice_list)  # 挑选线下测试集


offline_test = train_feat[idx ]  # 线下测试集
train_xy = train_feat[(~idx)]  # 训练集
# train_xy = train_feat

label_train = train_xy['Ki']
label_test = offline_test['Ki']
# 3.2 剔除个案标记列 id
# train_xy 训练集  testt_feat 线上测试集  offline_test 线下测试集
submission = testt_feat[['Protein_ID', 'Molecule_ID']]
train_xy = train_xy.drop('Ki', axis=1)
train_xy = train_xy.drop('Protein_ID', axis=1)
train_xy = train_xy.drop('Molecule_ID', axis=1)

testt_feat = testt_feat.drop('Ki', axis=1)
testt_feat = testt_feat.drop('Protein_ID', axis=1)
testt_feat = testt_feat.drop('Molecule_ID', axis=1)
#
offline_test = offline_test.drop('Ki', axis=1)
offline_test = offline_test.drop('Protein_ID', axis=1)
offline_test = offline_test.drop('Molecule_ID', axis=1)

print('训练集', train_xy.shape)
print('线下测试集', offline_test.shape)
print('线上预测集', testt_feat.shape)

# 3.3 lgb算法
# 3.3.1 导入数据集
train = lgb.Dataset(train_xy, label=label_train)
test = lgb.Dataset(offline_test, label=label_test, reference=train)

# 3.3.2 训练配置
params = {
    'task': 'predict',
    'boosting_type': 'gbdt',
    'objective': 'regression_l2',
    # 'metric': ['l2', 'auc'],
    'metric': ['l2'],
    # 'min_child_weight': 3,
    # 'max_depth': 4,
    'num_leaves': 85,
    'is_unbalance': True,
    'lambda_l2': 5,
    # 'feature_fraction':0.8,
    # 'subsample': 0.7,
    # 'colsample_bytree': 0.7,
    'learning_rate': 0.05,
    'seed': 2000,
    'nthread': 4,
}
print('当前参数:', params)

# 3.3.3 训练
num_round = 8000
gbm = lgb.train(params,
                train,
                # num_boost_round=boost_rounds,
                num_round,
                verbose_eval=50,
                valid_sets=[train, test]
                # valid_sets=[train]
                )


# 3.3.4 预测


preds_offline = gbm.predict(offline_test)

print('线下mse ：', mean_squared_error(label_test, preds_offline))

preds_sub = gbm.predict(testt_feat)
# #
# # # 4 保存结果
nowTime = datetime.datetime.now().strftime('%m%d%H%M')  # 现在
name = '../submit/lgb_' + nowTime + '.csv'
submission['Ki'] = preds_sub
submission.to_csv(name, index=False)



import time
import os
import datetime

import numpy as np
import pandas as pd
import lightgbm as lgb
# from opt_func import import_data

from sklearn.model_selection import KFold, StratifiedKFold

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error
from contextlib import contextmanager
from lightgbm import LGBMRegressor
import gc

@contextmanager
def timer(title):
    t0 = time.time()
    yield
    print("{} - done in {:.0f}s".format(title, time.time() - t0))


from bayes_opt import BayesianOptimization

space = {"max_depth": (5, 8),
         "num_leaves": (20, 100),
         # "learning_rate": (0.01, 0.02), 
         "reg_alpha": (0, 5),
         "reg_lambda": (0, 15),
         "feature_fraction":(0.85, 1)
         }

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


feat = pd.read_csv('../output/protein_get_II.csv') # 不稳定系数
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

feat = pd.read_csv('../output/protein_get_Ai.csv') # 脂肪系数
print(feat.head())
data_feat = data_feat.merge(feat, on='Protein_ID', how='left')
del feat

feat = pd.read_csv('../output/protein_get_gravy.csv') # GRAVY系数
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
train_df = data[data['Ki'] > -11] # .fillna(0)
test_df = data[data['Ki'] <= -11] # .fillna(0)

del data

def kfold_lightgbm(max_depth, num_leaves,reg_alpha, reg_lambda, feature_fraction):
    # Divide in training/validation and test data

    gc.collect()
    # Cross validation model
    stratified = False
    num_folds = 5
    debug= False
    if stratified:
        folds = StratifiedKFold(n_splits= num_folds, shuffle=True, random_state=1001)
    else:
        folds = KFold(n_splits= num_folds, shuffle=True, random_state=1001)
    # Create arrays and dataframes to store results
    oof_preds = np.zeros(train_df.shape[0])
    sub_preds = np.zeros(test_df.shape[0])
    feature_importance_df = pd.DataFrame()
    feats = [f for f in train_df.columns if f not in ['Ki','Protein_ID','index', 'Molecule_ID']]
    
    for n_fold, (train_idx, valid_idx) in enumerate(folds.split(train_df[feats], train_df['Ki'])):
        train_x, train_y = train_df[feats].iloc[train_idx], train_df['Ki'].iloc[train_idx]
        valid_x, valid_y = train_df[feats].iloc[valid_idx], train_df['Ki'].iloc[valid_idx]

        # LightGBM parameters
        clf = LGBMRegressor(
            nthread=3,
            n_estimators=10000,
            learning_rate=0.08,
            num_leaves=int(num_leaves),
#             colsample_bytree=0.9497036,
#             subsample=0.8715623,
            max_depth=int(max_depth),
            reg_alpha=reg_alpha,
            reg_lambda=reg_lambda,
#             min_split_gain=0.0222415,
#             min_child_weight=39.3259775,
            feature_fraction=feature_fraction,
            silent=-1,
            verbose=-1, )

        clf.fit(train_x, train_y, eval_set=[(train_x, train_y), (valid_x, valid_y)], 
            eval_metric= 'mse', verbose= 100, early_stopping_rounds= 200)

        oof_preds[valid_idx] = clf.predict(valid_x, num_iteration=clf.best_iteration_)
        sub_preds += clf.predict(test_df[feats], num_iteration=clf.best_iteration_)  / folds.n_splits

        fold_importance_df = pd.DataFrame()
        fold_importance_df["feature"] = feats
        fold_importance_df["importance"] = clf.feature_importances_
        fold_importance_df["fold"] = n_fold + 1
        feature_importance_df = pd.concat([feature_importance_df, fold_importance_df], axis=0)
        del clf, train_x, train_y, valid_x, valid_y
        gc.collect()
    result = mean_squared_error(train_df['Ki'], oof_preds)

    # Write submission file 
    if not debug:
        nowTime = datetime.datetime.now().strftime('%m%d%H%M')  # 现在
        submission_file_name = "../output/tap_test_regress_pre"+ nowTime +".csv"
        test_df['Ki'] = sub_preds
        test_df[['Protein_ID', 'Ki']].to_csv(submission_file_name, index= False)

# with timer("Run LightGBM with kfold"):
#     b_lgb = BayesianOptimization(kfold_lightgbm, space, random_state=2000)
#     b_lgb.maximize()
#     print(b_lgb.res['max']) # 获取最大化时的值和参数

with timer("Run LightGBM with kfold"):
    kfold_lightgbm(8, 100,0, 0, 0.85)


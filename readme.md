本仓库为《基于人工智能的药物分子筛选》比赛三等奖源码， B榜第五 分数为1.35198

比赛链接：http://www.dcjingsai.com/common/cmpt/%E5%9F%BA%E4%BA%8E%E4%BA%BA%E5%B7%A5%E6%99%BA%E8%83%BD%E7%9A%84%E8%8D%AF%E7%89%A9%E5%88%86%E5%AD%90%E7%AD%9B%E9%80%89_%E7%AB%9E%E8%B5%9B%E4%BF%A1%E6%81%AF.html

相关数据可从链接处报名下载。

## 参考文献：

- 《蛋白质等电点的计算》 乐茂华  湛江师范学院 2009
- 《蛋白质等电点的理论计算》 陈安和 渝州大学学报 1991年第三期
- 《蛋白质分类问题的特征提取算法研究》 张振慧 2006 理学博士论文 
- 《蛋白质序列的数学描述及其应用》张艳萍 浙江理工大学 2010 
- 《蛋白质序列离散灰色模型及其在药物开发中的应用研究》 林卫中 2013 
- 《蛋白质序列一级结构图形构造及相似性分析》 许时超 浙江理工大学 2015 
- 《基于2D分子指纹和非平衡数据集药物与受体交互作用预测研究》 闵建亮 2014 
- 《Correlation between stability of a protein and its dipeptide composion: a novel approach for predicting in vivo stability of a  protein from its primary sequence》 Kunchur Guruprasad, B.V.Bhasker Reddy  1990 
- 《生物信息学手册》郝柏林 张淑誉  2000 
- 《生物信息学札记》 樊龙江  浙江大学生物科学研究所、沃森基因组科学研究院 IBM生物计算实验室 
-  ExPASY ProtParam tool —— Bioinformatics Rescurce Portal

## 比赛思路

根据预测集发现，预测集和训练集的蛋白质ID不重合，药物分子有重合，因此，精力集中于蛋白质特征的挖掘。

#### 对药物分子特征：

缺失值填充、分子指纹直接字符串切断生成指纹矩阵作为特征

#### 对于蛋白质特征：

根据氨基酸序列计算并构建蛋白质理化性质：脂肪系数、grand average of hydropathicity (GRAVY)、不稳定系数、相对分子质量、Pi值 、 晓光系数。

根据计算机图形学原理，将一维氨基酸序列抽象拉取为三维空间向量特征。

#### 模型：

由于数据集分布原因， 笔者在此次比赛中没有在模型中下大力气，仅仅采用lightGBM手工调参而已。

## 软件包说明：

python 3.6.3
numpy 1.13.3
pandas 0.20.3
sklearn 0.19.1
lightgbm 2.1.1

## 文件夹说明

.
├── code  源代码文件夹
│   ├── molecule_feature.py      药物特征提取
│   ├── protein_getC.py          蛋白质理化性质特征提取
│   ├── protein_getp.py          蛋白质二维图形化表示特征提取
│   └── yaowu_data_for_train.py  模型训练与提交结果生成
├── input  竞赛输入文件夹
│   ├── df_affinity_test_toBePredicted.csv
│   ├── df_affinity_train.csv
│   ├── df_molecule.csv
│   ├── df_protein_test.csv
│   └── df_protein_train.csv
├── output 源代码中间生成的特征数据文件夹
│   ├── molecule_feature.csv  药物特征
│   ├── protein_getC.csv      蛋白质三种理化性质数据
│   └── protein_getp.csv      蛋白质三种二维图形化表示特征数据
├── readme.txt       代码运行和文件夹说明
└── submit           提交结果文件所在目录
  └── lgb_******.csv

## 代码运行说明：

需要先运行code文件夹中的三个特征获取代码文件：
molecule_feature.py    会获取药物特征  生成文件到output文件夹的molecule_feature.csv
protein_getC和protein_getp 会获取蛋白质相关特征 生成文件到output文件夹的protein_getC.csv和protein_getp.csv

待以上output文件生成后，
运行code文件夹中的yaowu_data_for_train.py文件，进行模型训练和预测；
生成submit文件夹下的以lgb开头的文件即为提交文件。



#!/Users/huangshi/anaconda3/envs/qiime2-2020.2/bin/python

__author__ = "Shi Huang"
# coding: utf-8

# # Microbiome age prediction in the new datasets
# ### May 2020
# ## Motivation 
# We published the age-prediction models on **healthy** gut, oral and skin microbiomes using Random Forest regression analyses.
# This jupyter notebook attempted to apply our microbiome age models to new datasets using Q2 API and potentially link the microbiome age to more phenotypes.
# 
# ### Reference
# Huang S, Haiminen N, Carrieri A-P, Hu R, Jiang L, Parida L, Russell B, Allaband C, Zarrinpar A, Vázquez-Baeza Y, Belda-Ferre P, Zhou H, Kim H-C, Swafford AD, Knight R, Xu ZZ. 2020. Human skin, oral, and gut microbiomes predict chronological age. mSystems 5:e00630-19. https://doi.org/10.1128/mSystems.00630-19.
# 
# ### Qiita study IDs involved in the trained model: 
# * Gut microbiota:
# 
# | QIITA Study ID | EBI accession ID | Project name | Publication(s) | # of samples involved |
# | ------------------ | ------------------ | ------------------ |------------------ | ------------------ |
# |[10317](https://qiita.ucsd.edu/study/description/10317)| ERP012803 | American Gut Project | [American Gut: an Open Platform for Citizen Science Microbiome Research](https://msystems.asm.org/content/3/3/e00031-18) | 2770 |
# |[11757](https://qiita.ucsd.edu/study/description/11757)| PRJEB18535 | GGMP regional variation | [Regional variation greatly limits application of healthy gut microbiome reference ranges and disease models](https://www.nature.com/articles/s41591-018-0164-x)| 1609 |
# 
# 
# ## About the implementation
# We re-trained the microbiome-age model using `q2-sample-classifer` which generated a Q2 artifact `SampleEstimator[Regressor]` for your applications. 

import numpy as np
import pandas as pd
import qiime2 as q2
import os
from optparse import OptionParser

from biom import Table
from qiime2 import Artifact
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from qiime2.plugins.sample_classifier.actions import predict_regression, regress_samples, scatterplot

usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--test_data_fp",
                  metavar="FILE", action="store", type="string", help="The file path of test feature table (qza)")
parser.add_option("-b", "--test_sample_md_fp",
                  metavar="FILE", action="store", type="string", help="The file path of test metadata (tsv)")    
parser.add_option("-c", "--test_prefix", default="test", 
                  metavar="FILE", action="store", type="string", help="The prefix of test data (str) [default: %default]")  
parser.add_option("-d", "--test_target_field",
                  action="store", type="string", help="The target field in the test metadata (str)")
parser.add_option("-o", "--test_outdir",
                  action="store", default="Test_out", type="string", help="The output directory for the age predictions [default: %default]")
parser.add_option("-A", "--train_data_fp",
                  metavar="FILE", action="store", type="string", help="The file path of train feature table (qza)")
parser.add_option("-B", "--train_sample_md_fp",
                  metavar="FILE", action="store", type="string", help="The file path of train metadata (tsv)")    
parser.add_option("-C", "--train_prefix", default="train", 
                  metavar="FILE", action="store", type="string", help="The prefix of train data (str) [default: %default]")  
parser.add_option("-D", "--train_target_field", type="string", 
                  action="store", help="The target field in the train metadata (str)")
parser.add_option("-T", "--retrain", 
                  action="store_true", default="True", help="If the model should be retained: [default: %default]")
parser.add_option("-O", "--train_outdir", type="string", 
                  action="store", default="Train_out", help="The output directory for the trained model [default: %default]")         
opts, args = parser.parse_args()


print(""" source activate qiime2-2020.2
python q2--MicrobiomeAgePrediction4NewDatasets.py \
--test_data_fp Test_data/CIAO/100nt/78846_feature-table.qza \
--test_sample_md_fp Test_data/CIAO/metadata_192.tsv \
--test_target_field agevisit \
--test_prefix CIAO \
--test_outdir Test_data/CIAO/ \
--train_data_fp Train_data/16S-100nt-gut_4434/gut_4434.qza \
--train_sample_md_fp Train_data/16S-100nt-gut_4434/gut_4434_map.txt \
--train_target_field age \
--retrain True \
--train_outdir Regressor/16S-100nt-gut_4434/ \
""")

# ## Input the test data table and metadata
test_data_fp = opts.test_data_fp #'Test_data/CIAO/100nt/78846_feature-table.qza'
test_sample_md_fp = opts.test_sample_md_fp #'Test_data/CIAO/metadata_192.tsv'
test_prefix = opts.test_prefix #'CIAO'
test_target_field = opts.test_target_field # 'agevisit'
OUTDIR= opts.test_outdir #'Test_data/CIAO/'

# ## Input the train data table, metadata, and a prebuilt Q2 RandomForestRegressor

train_data_fp= opts.train_data_fp#'Train_data/16S-100nt-gut_4434/gut_4434.qza' #'Train_data/shotgun-gotu-Finrisk/gotu.shared.feature-table.qza' 
train_sample_md_fp= opts.train_sample_md_fp #'Train_data/16S-100nt-gut_4434/gut_4434_map.txt'#'Train_data/shotgun-gotu-Finrisk/gotu.shared.metadata.txt'
train_target_field= opts.train_target_field #'age' #'FIBER_TOTAL' 
retrain=opts.retrain # False
trained_model_fp=opts.train_outdir #'Regressor/16S-100nt-gut_4434/' #'Regressor/shotgun-gotu-Finrisk-FIBER_TOTAL-ht/' 


try:
    if trained_model_fp is None:
        print("The output directory of trained model should be set!")
except NameError:
    print ("This output directory can not be defined")

try:
    if not os.path.exists(trained_model_fp):
        os.mkdir(trained_model_fp)
        print("The output directory of trained model has been created!")
except NameError:
    print ("The output directory of trained model has already existed!")

# ## Load the test data files
# test_X = q2.Artifact.load(test_data_fp).view(pd.DataFrame)
test_X_q2 = q2.Artifact.load(test_data_fp)
test_metadata=pd.read_csv(test_sample_md_fp, sep='\t', index_col=0)
try:
    test_y=test_metadata[test_target_field]
except NameError:
    print ("The test_target_field is not found!")
test_metadata_q2=q2.Metadata(test_metadata)
test_y_q2=test_metadata_q2.get_column(test_target_field)

# ## Load the train data files
train_X_q2 = q2.Artifact.load(train_data_fp) #q2 feature table
train_metadata=pd.read_csv(train_sample_md_fp, sep='\t', index_col=0) # read the metadaata file
try:
    train_y=train_metadata[train_target_field]
except NameError:
    print("The train_target_field is not found!")
train_metadata_q2=q2.Metadata(train_metadata) # q2 metadata
train_y_q2=train_metadata_q2.get_column(train_target_field)

# ## Load the Q2 RandomForestRegressor
# ### Option 1: load the pre-built model 
# ### Option 2: re-train the model using the train data table using `regress-samples`

if retrain==True:
    out=regress_samples(train_X_q2, train_y_q2, cv=5, n_jobs=8, n_estimators=500, parameter_tuning=False)
    trained_model_q2=out.sample_estimator
    out.sample_estimator.save(trained_model_fp+'sample_estimator.qza')
    out.feature_importance.save(trained_model_fp+'feature_importance.qza')
    out.predictions.save(trained_model_fp+'predictions.qza')
    out.model_summary.save(trained_model_fp+'model_summary.qzv')
    out.accuracy_results.save(trained_model_fp+'accuracy_results.qzv')

else:
    trained_model_q2=q2.Artifact.load(trained_model_fp+'sample_estimator.qza')


# ### You can view the model performance based on the outputs from `q2-sample-classifier`

# ## The essential preprocessing steps for the test table
# ## (1) The normalization of the ASV feature format
# ### Problem:
# The train data only contains 100-nt sequence features.
# 
# For example, the test data contains 150-nt sequence features or others, which will prevent this data from microbiome age prediction based on the train data.
# ### Solution: 
# We will truncate the 150-nt sequences into 100 nt ones and collapse all counts of identical ASVs after this step.

def trim_asvs_to_length(test_data, start=0, end=100):
    '''
    Parameters
    -------
        x: str
        The file path specify a qza file that contains sequence-like features in the columns
    Return
    -------
        x_dedup: q2 artifact 
        A table that contain sequence-like features with desired length
    Examples
    -------
    x=pd.DataFrame({'atcttc':[1, 3, 1, 3], 'ttcttc':[1, 3, 3, 1], 
                    'aatttc':[2, 5, 3, 1], 'ttcttc':[2, 5, 3, 1],
                    'aattcc':[2, 5, 3, 1], 'aatatc':[2, 0, 0, 1]})

    '''
    x=test_data.view(pd.DataFrame)
    ids=x.columns.tolist()
    all_length_equal_to_100=all([len(i)==100 for i in ids])
    if(all_length_equal_to_100):
        x_dedup=x
    else:            
        new_ids=[i[start:end] for i in ids]
        x.columns=new_ids
        def checkIfDuplicates(listOfElems):
            ''' Check if given list contains any duplicates '''
            if len(listOfElems) == len(set(listOfElems)):
                return False
            else:
                return True
        if(checkIfDuplicates(new_ids)):
            x_dedup=x.sum(axis=1, level=0)
        else:
            x_dedup=x
    x_dedup_qza=q2.Artifact.import_data('FeatureTable[Frequency]', x_dedup)

    return x_dedup_qza

test_X_q2=trim_asvs_to_length(test_X_q2)

# ## (2) The alignment of the ASV features from the train and test datasets
# ### Problem:
# The test data usually will not have a identical set of ASV features as that in the train data. 
# ### Solution: 
# We will only keep the test features consistent with those in the train data, and pad other train features with zeros in the test table to ensure the test table has the same columns with train data finally.

def pad_features_by_qza(train_data, test_data):
    '''
    Parameters
    ----------
    train_datafile : Q2 feature-table artifact i.e. 'FeatureTable[Frequency]'
        The train data table, 
    test_datafile : Q2 feature-table artifact i.e. 'FeatureTable[Frequency]'
        The test data table, 
    Returns
    -------
    new_b_qza: 'FeatureTable[Frequency]'
        A updated test data table with equal number of
        feature as the train table.
    '''
    a=train_data.view(pd.DataFrame)
    b=test_data.view(pd.DataFrame)
    a_feature_ids=a.columns.values.tolist()
    b_feature_ids=b.columns.values.tolist()
    print("The # of features in the train data: ", len(a_feature_ids))
    print("The # of features in the original test data: ", len(b_feature_ids))
    a_uniq_f=list(set(a_feature_ids)-set(b_feature_ids))
    ab_shared_f=set(a_feature_ids).intersection(set(b_feature_ids))
    print("The # of features with all zeros in the new test data: ", len(a_uniq_f))
    print("The # of shared features kept in the new test data: ", len(ab_shared_f))
    b_padding_matrix = pd.DataFrame(0, index=b.index, columns=a_uniq_f)
    new_b=pd.concat([b[ab_shared_f], b_padding_matrix], axis=1)
    new_b=new_b[a_feature_ids]    
    print("The shape of new test data: ", new_b.shape)
    new_b_qza=q2.Artifact.import_data('FeatureTable[Frequency]', new_b)
    return new_b_qza

# ## Align features in the test dataset with those in the train data

test_X_padding_qza=pad_features_by_qza(train_X_q2, test_X_q2)

# ##  Microbiome age prediction using `predict_regression`

predictions=predict_regression(test_X_padding_qza, trained_model_q2).predictions
predictions.save(OUTDIR+'test_predictions.qza')

test_pred_df=predictions.view(pd.Series)

result = pd.concat([test_metadata, test_pred_df], axis=1, sort=False)
result.to_csv(OUTDIR+'test_predictions_metadata.tsv',sep='\t')


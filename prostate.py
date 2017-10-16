from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, accuracy_score, roc_curve
from sklearn.linear_model import LogisticRegression
%matplotlib inline

def plot_roc(y_test, y_pred_prod, name):
    '''
    Using sklearn roc_curve plot roc curve
    INPUT:
    y_test: Array of true labels
    y_pred_prod: Array of probabilities of target variable
    OUTPUT:
    None
    '''
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_prod)
    plt.plot(fpr, tpr, label = name)
    plt.rcParams['font.size'] = 12
    plt.plot([0,1],[0,1],color='black',ls='--')
    #plt.title('ROC curve for Churn Classifier')
    plt.xlabel('Rate of False Diagnosis of Cancer')
    plt.ylabel('Rate of True Diagnosis of Cancer')
    plt.grid(False)

def report_metrics(y_test, y_pred_lbl):
    '''
    Return Precision/Recall/accuracy
    INPUT: y_test: Array of true labels
       y_pred_lbl: Array of predicted labels
    OUTPUT: Return precision, recall and accuracy score values
    '''
    precision = precision_score(y_test, y_pred_lbl)
    recall = recall_score(y_test, y_pred_lbl)
    accuracy = accuracy_score(y_test, y_pred_lbl)
    return precision, recall, accuracy

logistic_params = {'logistic__penalty': ['l1','l2'],
                    'logistic__C':[0.05, 0.1, 0.2, 0.5,1],
                    'logistic__random_state': [1]}

df = pd.read_csv('tcga.csv')
df = df[df['organ'] == 'Prostate']

df2 = pd.read_csv('gtex_sample_expression.csv')
df2 = df2[df2['tissue'] == 'prostate']
df1_test = df.groupby('gene_id').median().reset_index()
df2_test = df2.groupby('gene_id').median().reset_index()
gene_df = df2_test.merge(df1_test, how='inner', on = ['gene_id'])
gene_df['half_rpkm'] = gene_df['rpkm_expression']*.5
workable_df = gene_df[(gene_df['fpkm_expression'] < gene_df['half_rpkm']) | (gene_df['fpkm_expression'] > gene_df['rpkm_expression'])]
workable_df['final'] = workable_df['half_rpkm']
mask = (workable_df['fpkm_expression']-workable_df['half_rpkm']) > (workable_df['rpkm_expression']-workable_df['fpkm_expression'])
mask2 = (workable_df['fpkm_expression']-workable_df['half_rpkm']) < (workable_df['rpkm_expression']-workable_df['fpkm_expression'])
fpkm = workable_df[mask][['gene_id','fpkm_expression','sample_number']]

rpkm = workable_df[mask2][['gene_id','rpkm_expression','sample_number']]
fpkm['label'] = 1
rpkm['label'] = 0
fpkm['final'] = fpkm['fpkm_expression']
rpkm['final'] = rpkm['rpkm_expression']
final_df = pd.concat([fpkm,rpkm])
final_df = final_df[['final', 'sample_number','gene_id','label']]

matrix = defaultdict(list)

def matrix_transform(final_df):
    for row in final_df.values:
        matrix[row[1]].append((row[0],row[3]))
    len(sorted(matrix.values(), key=lambda v: len(v), reverse = True)[0])

def sort_matrix(final_df):
    for key,value in matrix.items():
        matrix[key] = sorted(value,key = lambda x: x[0])

labels = []
for value in matrix.values():
    labels.append(value[0][1])

for i,j in matrix.iteritems():
    if len(j)<17:
        num_zeros = 17-len(j)
        toAppend = [0]*num_zeros
        j.extend(toAppend)

for key, value in matrix.items():
    matrix[key] = [val[0] for val in value]

X = np.array(X)
X.shape

X_train, X_test, y_train, y_test = train_test_split(X,yl)

l = LogisticRegression()
l.fit(X_train, y_train)
predictions = l.predict(X_test)

predict_label = l.predict_proba(X_test)[:,1]
plot_roc(y_test,predict_label,'l')

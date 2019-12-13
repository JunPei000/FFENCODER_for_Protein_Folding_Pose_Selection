################################################################
#                                                              # 
#   Set up Random Forest -- grid search for protein folding    #
#   system                                                     #
################################################################

import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score
from random import shuffle
from sklearn import metrics

direct2 = PATH OF FINAL DESCRIPTOR FILES
direct3 = 'important_f_amber_ff14/'
impo = open(direct3+'amber_ff14_with_oop_500.txt','r')
imp = impo.readlines()
important_f = []
for line in imp:
    m = line.split()
    for i in m:
        if i not in important_f:
            important_f.append(i)
print (len(important_f))
print (important_f)
feature_to_drop = []
for column in [str(i) for i in range(0,19512)]:
    if column != '19511' and column not in important_f:
        feature_to_drop.append(column)
#feature_to_drop.append('Unnamed: 0.1')
print (len(feature_to_drop))


#filenames = shuffle([i for i in os.listdir(direct2) if '.csv' in i])
filenames = [i for i in os.listdir(direct2) if '.csv' in i and 'rosetta' not in i]
shuffle(filenames)
trainfiles, testfiles = train_test_split(filenames, train_size = 0.8, random_state = 42) 
print (testfiles)
x_train = pd.DataFrame()
y_train = pd.DataFrame()
x_test = pd.DataFrame()
y_test = pd.DataFrame()
train_data = pd.DataFrame()
test_data = pd.DataFrame()
for tfile in trainfiles:
    train_data = train_data.append(pd.read_csv(direct2+tfile).drop(feature_to_drop,axis=1))
if 'Unnamed: 0.1' in train_data.columns:
    train_data = train_data.drop(['Unnamed: 0.1'], axis=1)
x_train = train_data.drop(['19511'], axis=1).iloc[:,1:]
y_train = train_data['19511']
y_train=y_train.astype('float')
print(x_train)
print(y_train)

for tfile in testfiles:
    test_data = test_data.append(pd.read_csv(direct2+tfile).drop(feature_to_drop,axis=1))
    temp = pd.read_csv(direct2+tfile)[1:]
    for row in temp.index:
        if float(temp.loc[row, '19511']) != 0 and float(temp.loc[row, '19511']) != 1:
            print (tfile, row, temp.loc[row, '19511'])
if 'Unnamed: 0.1' in test_data.columns:
    test_data = test_data.drop(['Unnamed: 0.1'], axis=1)
x_test = test_data.drop(['19511'], axis=1).iloc[:,1:]
y_test = test_data['19511']
print(x_test)
print(y_test)
y_test=y_test.astype('float')

clf =  RandomForestClassifier(criterion='gini', class_weight = {0:0.5, 1:0.5}, n_jobs=30)
parameters = {
       'n_estimators':(1000,5000),
       'max_depth':(5,100),
       'min_samples_split':(2,5),
       'min_samples_leaf':(1,2) }
grid_search = GridSearchCV(clf, parameters, cv=5,verbose=1,scoring='accuracy')
grid_search.fit(x_train, y_train)
print ('Best Training score: %0.3f' % grid_search.best_score_)
print ('Best parameters set:')
best_parameters = grid_search.best_estimator_.get_params()
for param_name in sorted(parameters.keys()):
    print ('\t%s: %r' % (param_name, best_parameters[param_name]))
predictions = grid_search.predict(x_test)
print ("Testing accuracy:",round(accuracy_score(y_test,predictions),4))
print ("\nComplete report of Testingdata\n",classification_report(y_test, predictions))
print ("\n\nRandom Forest Grid Search- Test Confusion Matrix\n\n",pd.crosstab( y_test, predictions,rownames = ["Actuall"],colnames = ["Predicted"]))
probs = grid_search.predict_proba(x_test)
preds = probs[:,1]
fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
roc_auc = metrics.auc(fpr, tpr)
print (type(fpr), type(tpr), type(roc_auc))
newfile1 = open(direct2+'fpr_IMP_500_with_all.txt', 'w')
newfile2 = open(direct2+'tpr_IMP_500_with_all.txt', 'w')
newfile3 = open(direct2+'auc_IMP_500_with_all.txt', 'w')
for i in range(len(fpr)):
    newfile1.write(str(fpr[i])+'\n')
    newfile2.write(str(tpr[i])+'\n')
newfile3.write(str(roc_auc))
newfile1.close()
newfile2.close()
newfile3.close()
test_data = pd.DataFrame()
x_test = pd.DataFrame()

for tfile in testfiles:
    print (tfile)
    test_data = pd.read_csv(direct2+tfile)
    print (test_data)
    if 'Unnamed: 0' in test_data.columns:
        test_data = test_data.drop('Unnamed: 0', axis=1)
    x_test = test_data.drop(['19511','Unnamed: 0.1'], axis=1).astype('float').drop(feature_to_drop,axis=1)
    print (x_test)
    names = []
    scores = []
    rows = [2*i for i in range(0, int(len(x_test.index)/2))]
    scn = 0
    for row in rows:
        nums = []
        nums = [i for i in rows if i != row]
        x_testd = pd.DataFrame(columns= x_test.columns, index=[d1 for d1 in range(0, len(rows))])
        x_testd1 = pd.DataFrame(columns= x_test.columns, index=[d1 for d1 in range(0, len(rows))])
        x_testd2 = pd.DataFrame(columns= x_test.columns, index=[d1 for d1 in range(0, len(rows))])

        sc = 0
        if row+1 in x_test.index:
            predictionsn = grid_search.predict(x_test.loc[row+1].values.reshape(1, -1))
            scn += predictionsn
            print ('scn=', scn)
        names.append(test_data.loc[row, 'Unnamed: 0.1'])
        nums = []
        nums = [i for i in rows if i != row]
        print (len(rows), len(nums))
        x_testd1.loc[[d1 for d1 in range(0, len(rows))],:] = x_test.loc[row].values
        x_testd2.loc[[d1 for d1 in range(0, len(rows)-1)],:] = x_test.loc[[d2 for d2 in nums],:].values
        x_testd2 = x_testd2.fillna(0)
        x_testd=x_testd1.sub(x_testd2)
        predictions = grid_search.predict(x_testd)
        sc = sum(predictions)
        print (sc)
        scores.append(sc)
    scores.append(scn)
    names.append('native')
    print ('native ranking is : ', len(scores)-sorted(scores).index(scores[-1]))
    print (len(scores))
    print (scores)
    print (names)
    print ('best decoy is : ', names[:-1][scores[:-1].index(max(scores[:-1]))])


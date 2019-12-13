#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:23:59 2018

@author: peijun
"""

import os
import pandas as pd

direct1 = PATH OF ORIGINAL OUTPUT FILE FROM FFENCODER
direct_oop = PATH OF ORIGINAL OOP OUTPUT FILE FROM FFENCODER_OOP
direct2 = PATH OF THE FINAL DESCRIPTOR FILE 
if not os.path.isdir(direct2):
    os.makedirs(direct2)
aa = 'lmds_v2'


filenames = [i for i in os.listdir(direct1) if '.csv' in i and '.DS' not in i and aa in i
             and 'descriptort' in i and 'unk' not in i and 'l30' not in i and 'casp3_' not in i]
#filenames = [i for i in os.listdir(direct1) if '.csv' in i and '.DS' not in i and 'descriptort' in i]
for filename in filenames:
    print (filename)
    torsionpart = pd.read_csv(direct1+filename)
    nonbondpart = pd.read_csv(direct1+filename.replace('descriptort', 'descriptorn'))
    torsion_oop = pd.read_csv(direct_oop+filename.replace('.csv', '_oop.csv'))
    sysname = filename.replace(aa+'_', '')
    #runnum = sysname.replace('.csv', '')[-1]
    sysname = sysname[:sysname.index('descri')-1]
    if '_' in sysname:
        sysname = sysname.replace('_', '-')
    if sysname == '1shf':
        sysname = '1shf-A'
    #sysname = 'rosetta_natives'
    #sysname = 'native'
    #nativename = 'model1.pdb'
    nativename = sysname+'.pdb'
    #nativename = [i for i in torsionpart.columns if 'native' in i][0]   
    #nativesys = [i for i in torsionpart.columns if '.out' in i]
    newdescriptort = pd.DataFrame()
    newdescriptorn = pd.DataFrame()
    newdesctor_oop = pd.DataFrame()
    torsionpart = torsionpart.fillna(0)
    nonbondpart = nonbondpart.fillna(0)
    torsion_oop = torsion_oop.fillna(0)
    #torsionpart['avg'] = torsionpart[nativesys].mean(axis = 1)
    #nonbondpart['avg'] = nonbondpart[nativesys].mean(axis = 1)
    for column in torsionpart.columns:
        if 'Pair_name' in column:
            newdescriptort[column] = torsionpart[column].map(lambda x: x+'t')
            for row in nonbondpart.index:
                if nonbondpart.loc[row, column] == 0:
                    nonbondpart.loc[row, column] = 'sasa'
            newdescriptorn[column] = nonbondpart[column].map(lambda x: x+'n')
            newdesctor_oop[column] = torsion_oop[column].map(lambda x: x+'o')
        elif 'ref' in column:
            continue
         #   newdescriptort[column] = torsionpart[column]
          #  newdescriptorn[column] = nonbondpart[column]
        else:
            if column != nativename and 'Unnamed' not in column:
            #if '.out' not in column and 'Unnamed' not in column and 'avg' not in column:
                zero_y = column.replace('.pdb', '000')
                one_y = column.replace('.pdb','111')
                newdescriptort[zero_y] = torsionpart[column]- torsionpart[nativename]
                newdescriptort[one_y] = torsionpart[nativename] - torsionpart[column]
                newdescriptorn[zero_y] = nonbondpart[column]- nonbondpart[nativename]
                newdescriptorn[one_y] = nonbondpart[nativename] - nonbondpart[column]
                newdesctor_oop[zero_y] = torsion_oop[column] - torsion_oop[nativename]
                newdesctor_oop[one_y] = torsion_oop[nativename] - torsion_oop[column]
                #newdescriptort[zero_y] = torsionpart[column]- torsionpart['avg']
                #newdescriptort[one_y] = torsionpart['avg'] - torsionpart[column]
                #newdescriptorn[zero_y] = nonbondpart[column]- nonbondpart['avg']
                #newdescriptorn[one_y] = nonbondpart['avg'] - nonbondpart[column]
    for column in torsionpart.columns:
        if column not in nonbondpart.columns:
            print (column)
        if column not in torsion_oop.columns:
            print (column)
    for column in nonbondpart.columns:
        if column not in torsionpart.columns:
            print (column)
        if column not in torsion_oop.columns:
            print (column)
    newdescriptort = newdescriptort.fillna(0)
    newdescriptorn = newdescriptorn.fillna(0)
    newdesctor_oop = newdesctor_oop.fillna(0)
    fp1 = open(direct2+'/'+filename.replace('descriptort', 'detailed_atomp_gasphase'),'w')
    m  = pd.DataFrame()
    m = newdescriptort.append(newdescriptorn, ignore_index=True).append(newdesctor_oop, ignore_index=True)
    lr = len(m.index)
    for column in m.columns:
        if column[-3:] == '111':
            m.loc[lr, column] = 1
        if column[-3:] == '000':
            m.loc[lr, column] = 0
        if 'Pair_name' in column:
            m.loc[lr,column] = 'class'
    m = pd.DataFrame.transpose(m)
    print (len(m.columns))
    m.to_csv(fp1)
    fp1.close()


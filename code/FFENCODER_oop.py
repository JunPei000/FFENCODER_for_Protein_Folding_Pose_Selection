#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 16:00:25 2019

@author: peijun
"""

from protein_info_H import protein_info
import pandas as pd
import os
from insert import insert
import json
from out_of_plane import out_of_plane_calc_ff94, out_of_plane_calc_ff14

direct1 = PATH OF PROTEIN-FOLDING DECOY SET
newfile1 = open(direct1+'descriptort_amber_oop_ff94.csv','w')
newfile2 = open(direct1+'descriptort_amber_oop_ff14.csv','w')

filenames = [i for i in os.listdir(direct1) if '.pdb' in i and '.swp' not in i and 'inserted' not in i]

t1 = pd.read_csv('template/oop_ff94.csv')
t2 = pd.read_csv('template/oop_ff14.csv')
t1 = t1.drop(['Unnamed: 0'], axis=1)
t2 = t2.drop(['Unnamed: 0'], axis=1)

for pdbname in filenames:
    print (pdbname)
    insert(direct1, direct1, pdbname)
    pro_file = open(direct1+'inserted_'+pdbname, 'r').readlines()
    protein_oop_energy_ff94 = {}
    protein_oop_energy_ff14 = {}
    protein_oop_energy_ff94 = out_of_plane_calc_ff94(pro_file)
    protein_oop_energy_ff14 = out_of_plane_calc_ff14(pro_file)
    t1[pdbname] = t1['Pair_name'].map(protein_oop_energy_ff94)
    t2[pdbname] = t2['Pair_name'].map(protein_oop_energy_ff14)
t1 = t1.fillna(0)
t2 = t2.fillna(0)
t1.to_csv(newfile1)
t2.to_csv(newfile2)
newfile1.close()
newfile2.close()

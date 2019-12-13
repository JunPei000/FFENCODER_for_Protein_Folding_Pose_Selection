#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 16:00:25 2019

@author: peijun
"""

from protein_info_H import protein_info
from protein_info_finder import torsion_info_ff14
from tor_nonb_calc_ff14 import tor_calc
from tor_nonb_calc_ff14 import nonb_calc
from ref_des_build import ref_des_ff14
from change_key_to_standard import change_key_to_standard_amber_torsion as ckts_torsion
from change_key_to_standard import change_key_to_standard_amber_nonbond as ckts_nonbond
import pandas as pd
import os
from insert import insert
import json

direct1 = PATH OF THE PROTEIN FOLDING DECOY SET 
newfile1 = open(direct1+'descriptort_amber_ff14.csv','w')
newfile2 = open(direct1+'descriptorn_amber_ff14.csv','w')
filenames = [i for i in os.listdir(direct1) if '.pdb' in i and '.swp' not in i and 'inserted' not in i]


[dictor, dicnon] = ref_des_ff14()
t = pd.read_csv('template/descriptort_ff14.csv')
n = pd.read_csv('template/descriptorn_ff14.csv')
t = t.drop(['Unnamed: 0'], axis=1)
n = n.drop(['Unnamed: 0'], axis=1)
for pdbname in filenames:
    print (pdbname)
    insert(direct1, direct1, pdbname)
    pro_file = open(direct1+'inserted_'+pdbname, 'r').readlines()
    pro_info = protein_info(pro_file)
    pro_tor = torsion_info_ff14(pro_info[0], pro_info[1])
    dicT = tor_calc(pro_tor)
    dicN = nonb_calc(pro_info[3])
    protein_tor_energy = {}
    protein_nonb_energy = {}
    protein_tor_energy = ckts_torsion(dictor, dicT)
    protein_nonb_energy = ckts_nonbond(dicnon, dicN)
    t[pdbname] = t['Pair_name'].map(protein_tor_energy)
    n[pdbname] = n['Pair_name'].map(protein_nonb_energy)
t = t.fillna(0)
n = n.fillna(0)
t.to_csv(newfile1)
n.to_csv(newfile2)
newfile1.close()
newfile2.close()

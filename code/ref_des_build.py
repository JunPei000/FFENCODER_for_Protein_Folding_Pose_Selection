#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:13:33 2018

@author: peijun
"""
import numpy as np
import pandas as pd
import json

def ref_des_ff94():
    torsion = {}
    t = {}
    fpt = open('Amber_ff94_parameter/final_torsion_type.json','r')
    torsion = json.load(fpt)
    for key in torsion:
        if not key in t:
            t[key] = 1.0
    descriptort = pd.DataFrame(list(t.items()), columns=['Pair_name', 'ref'])
    #descriptort = pd.DataFrame.from_dict(t.items(), orient='index', columns=['Pair_name', 'ref']) 

    nonb_file = open('Amber_ff94_parameter/final_nonb_type.json','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if not key in n:
            n[key] = 1.0
    descriptorn = pd.DataFrame(list(n.items()), columns=['Pair_name', 'ref'])
    #descriptorn = pd.DataFrame.from_dict(n.items(), orient='index', columns=['Pair_name', 'ref'])
    #return descriptort, descriptorn,t,n
    return t, n

def ref_des_ff14():
    torsion = {}
    t = {}
    fpt = open('Amber_ff14_parameter/final_torsion_type.json','r')
    torsion = json.load(fpt)
    for key in torsion:
        if not key in t:
            t[key] = 1.0
    #descriptort = pd.DataFrame(t.items(),columns=['Pair_name', 'ref'])
    
    nonb_file = open('Amber_ff14_parameter/final_nonb_type.json','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if not key in n:
            n[key] = 1.0
    #descriptorn = pd.DataFrame(n.items(), columns=['Pair_name', 'ref'])
    #return descriptort, descriptorn,t,n
    return t, n






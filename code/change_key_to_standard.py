#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:52:44 2018

@author: peijun
"""
import json
from function import name
def change_key_to_standard(dic1,dic2):
    dic3 = {}
    for key in dic2:
        if dic1.has_key(key):
            if not dic3.has_key(key):    
                dic3[key] = dic2[key]
            elif dic3.has_key(key): 
                dic3[key] += dic2[key]
        elif not dic1.has_key(key):
            A = []
            A = name(key)
            #if 'OXT' in A[0]:
            #    A[0] = A[0].replace('OXT','O')
            #if 'OXT' in A[1]:
            #    A[1] = A[1].replace('OXT','O')
            ind0 = A[0].index('-')
            ind1 = A[1].index('-')
            if 'HIE' in A[0]:
                A[0] = A[0].replace('HIE','HID')
            elif 'HIP' in A[0]:
                A[0] = A[0].replace('HIP','HID') 
            elif 'HIS' in A[0]:
                A[0] = A[0].replace('HIS','HID')
            if 'HIE' in A[1]:
                A[1] = A[1].replace('HIE','HID')
            elif 'HIP' in A[1]:
                A[1] = A[1].replace('HIP','HID')
            elif 'HIS' in A[1]:
                A[1] = A[1].replace('HIS','HID')
            if len(A[0][ind0+1:]) == 3:
                resi1 = A[0][ind0+1:] 
            elif len(A[0][ind0+1:]) == 4:
                resi1 = A[0][ind0+2:] 
            if len(A[1][ind1+1:]) == 3:
                resi2 = A[1][ind1+1:] 
            elif len(A[1][ind1+1:]) == 4:
                resi2 = A[1][ind1+2:] 
            key_old = A[0].replace(A[0][ind0+1:],resi1)+'__'+A[1].replace(A[1][ind1+1:],resi2)
            key_new = A[1].replace(A[1][ind1+1:],resi2)+'__'+A[0].replace(A[0][ind0+1:],resi1)
            if dic1.has_key(key_old):
                if not dic3.has_key(key_old):    
                    dic3[key_old] = dic2[key]
                elif dic3.has_key(key_old): 
                    dic3[key_old] += dic2[key]
            elif dic1.has_key(key_new):
                if not dic3.has_key(key_new):    
                    dic3[key_new] = dic2[key]
                elif dic3.has_key(key_new): 
                    dic3[key_new] += dic2[key]
            elif not dic1.has_key(key_new) and not dic1.has_key(key_old):
                print ('undefined atom pair!', key, key_new, key_old, A, resi1, resi2)
    return dic3


def change_key_to_standard_amber_torsion(dic1,dic2):
    dic3 = {}
    for key in dic2:
        #if key in dic1:
        #    if key not in dic3:
        #        dic3[key] = dic2[key]
        #    elif key in dic3:
        #        dic3[key] += dic2[key]
        #elif key not in dic1:
        newlist = []
        newlist = key.split('_')
        #print newlist
        newkey1 = '_'.join([newlist[0], newlist[1], newlist[2], newlist[3]])
        newkey2 = '_'.join([newlist[0], newlist[2], newlist[1], newlist[3]])
        newkey3 = '_'.join([newlist[3], newlist[1], newlist[2], newlist[0]])
        newkey4 = '_'.join([newlist[3], newlist[2], newlist[1], newlist[0]])
        #print (newkey1)
        if newkey1 in dic1:
            if newkey1 not in dic3:
                dic3[newkey1] = dic2[key]
            elif newkey1 in dic3:
                dic3[newkey1] += dic2[key]
            continue
        elif newkey2 in dic1:
            if newkey2 not in dic3:
                dic3[newkey2] = dic2[key]
            elif newkey2 in dic3:
                dic3[newkey2] += dic2[key]
            continue
        elif newkey3 in dic1:
            if newkey3 not in dic3:
                dic3[newkey3] = dic2[key]
            elif newkey3 in dic3:
                dic3[newkey3] += dic2[key]
            continue
        elif newkey4 in dic1:
            if newkey4 not in dic3:
                dic3[newkey4] = dic2[key]
            elif newkey4 in dic3:
                dic3[newkey4] += dic2[key]
            continue
    return dic3

def change_key_to_standard_amber_nonbond(dic1,dic2):
    dic3 = {}
    for key in dic1:
        dic3[key] = 0
    for key in dic2:
        newlist = []
        newlist = key.split('__')
        key1 = newlist[1]+'__'+newlist[0]
        if key != key1:
            if key in dic1 and key1 not in dic1:
                dic3[key] += dic2[key]
            elif key not in dic1 and key1 in dic1:
                dic3[key1] += dic2[key]
            elif key not in dic1 and key1 not in dic1:
                print ('error! missing nonbond pair!', key, key1)
            elif key  in dic1 and key1 in dic1:
                print ('error! same nonbond pair exists!', key, key1)
        elif key == key1:
            if key in dic1:
                dic3[key] += dic2[key]
            elif key not in dic1:
                print ('error! missing nonbond atom pair!', key, key1)
    return dic3
                
            
    
                
                
                
                
                
                
                
                
                
                
                
                
                
        


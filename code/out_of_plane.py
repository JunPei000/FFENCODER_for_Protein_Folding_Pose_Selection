#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:17:18 2019

@author: peijun
"""
from __future__ import division
import os
from protein_info_H import protein_info
import json
from math import cos, pi
import numpy as np

def dihedral(p):
    """formula from Wikipedia article on "Dihedral angle"; formula was removed
    from the most recent version of article (no idea why, the article is a
    mess at the moment) but the formula can be found in at this permalink to
    an old version of the article:
    https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=689165217#Angle_between_three_vectors
    uses 1 sqrt, 3 cross products"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

    y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)
    return np.degrees(np.arctan2(y, x))

def modify_name(A):
    newlist = []
    if 'OXT' in A:
        A = A.replace('OXT', 'O')
    if 'HIE' in A:
        A = A.replace('HIE', 'HID')
    if 'HIP' in A:
        A = A.replace('HIP', 'HID')
    if 'HIS' in A:
        A = A.replace('HIS', 'HID')
    if '_' in A:
        newlist = A.split('_')
    elif '-' in A:
        newlist = A.split('-')
    if len(newlist[1]) == 4:
        if newlist[1][0] == 'C' or newlist[1][0] == 'N':
            newresi = newlist[1][1:]
        if '_' in A:
            return newlist[0]+'_'+newresi
        elif '-' in A:
            return newlist[0]+'-'+newresi
    else:
        return A

def out_of_plane_calc_ff94(pdbfile):
    resi_info = protein_info(pdbfile)[0]
    Bond_info = protein_info(pdbfile)[-1]
    resi_info_f = {}
    Bond_info_f = {}
    fp1 = open('Amber_ff94_parameter/out_of_plane_amber.json', 'r')
    out_of_plane = {}
    out_of_plane = json.load(fp1)
    fp2 = open('Amber_ff94_parameter/amber_atomtype_info.json', 'r')
    amber_atom_type = {}
    amber_atom_type = json.load(fp2)
    for i in range(len(resi_info)):
        for j in range(len(resi_info[i])):
            if resi_info[i][j][-1] not in resi_info_f:
                resi_info_f[resi_info[i][j][-1]] = {}
                resi_info_f[resi_info[i][j][-1]]['coor'] = [resi_info[i][j][1], resi_info[i][j][2], resi_info[i][j][3]]
                resi_info_f[resi_info[i][j][-1]]['KECSA_type'] = resi_info[i][j][0].replace('-', '_')
            elif resi_info[i][j][-1] in resi_info_f:
                print ('Error! Same atom exists in pdb! Please check: ', resi_info[i][j][-1])
    for i in range(len(Bond_info)):
        for j in range(len(Bond_info[i])):
            for k in range(len(Bond_info[i][j])):
                if Bond_info[i][j][k][0] not in Bond_info_f and Bond_info[i][j][k][1] not in Bond_info_f:
                    Bond_info_f[Bond_info[i][j][k][0]] = []
                    Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    Bond_info_f[Bond_info[i][j][k][1]] = []
                    Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                elif Bond_info[i][j][k][0] in Bond_info_f and Bond_info[i][j][k][1] not in Bond_info_f:
                    if Bond_info[i][j][k][1] not in Bond_info_f[Bond_info[i][j][k][0]]:
                        Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    Bond_info_f[Bond_info[i][j][k][1]] = []
                    Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                elif Bond_info[i][j][k][0] not in Bond_info_f and Bond_info[i][j][k][1] in Bond_info_f:
                    if Bond_info[i][j][k][0] not in Bond_info_f[Bond_info[i][j][k][1]]:
                        Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                    Bond_info_f[Bond_info[i][j][k][0]] = []
                    Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                elif Bond_info[i][j][k][0] in Bond_info_f and Bond_info[i][j][k][1] in Bond_info_f:
                    if Bond_info[i][j][k][1] not in Bond_info_f[Bond_info[i][j][k][0]]:
                        Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    if Bond_info[i][j][k][0] not in Bond_info_f[Bond_info[i][j][k][1]]:
                        Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                else:
                    print ('Error! Please check bond: ', Bond_info[i][j][k][1], Bond_info[i][j][k][0])
    out_of_plane_info = {}
    for key in Bond_info_f:
        if len(Bond_info_f[key]) >4:
            print ('Error! Please check bonds connect to ', key)
        elif len(Bond_info_f[key]) == 3:
            if modify_name(resi_info_f[Bond_info_f[key][0]]['KECSA_type']) == 'HE2_HID':
                continue
            if modify_name(resi_info_f[Bond_info_f[key][1]]['KECSA_type']) == 'HE2_HID':
                continue
            if modify_name(resi_info_f[Bond_info_f[key][2]]['KECSA_type']) == 'HE2_HID':
                continue
            if modify_name(resi_info_f[key]['KECSA_type']) == 'HE2_HID':
                continue
            amber_atom1 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][0]]['KECSA_type'])]['Amber_type']
            amber_atom2 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][1]]['KECSA_type'])]['Amber_type']
            amber_atom3 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][2]]['KECSA_type'])]['Amber_type']
            amber_atom4 = amber_atom_type[modify_name(resi_info_f[key]['KECSA_type'])]['Amber_type']
            amber_coor1 = resi_info_f[Bond_info_f[key][0]]['coor']
            amber_coor2 = resi_info_f[Bond_info_f[key][1]]['coor']
            amber_coor3 = resi_info_f[Bond_info_f[key][2]]['coor']
            amber_coor4 = resi_info_f[key]['coor']
            newkey1 = '_'.join([amber_atom1, amber_atom2, amber_atom4, amber_atom3])
            p1 = np.array([amber_coor1, amber_coor2, amber_coor4, amber_coor3])
            newkey2 = '_'.join([amber_atom2, amber_atom1, amber_atom4, amber_atom3])
            p2 = np.array([amber_coor2, amber_coor1, amber_coor4, amber_coor3])
            newkey3 = '_'.join([amber_atom1, amber_atom3, amber_atom4, amber_atom2])
            p3 = np.array([amber_coor1, amber_coor3, amber_coor4, amber_coor2])
            newkey4 = '_'.join([amber_atom3, amber_atom1, amber_atom4, amber_atom2])
            p4 = np.array([amber_coor3, amber_coor1, amber_coor4, amber_coor2])
            newkey5 = '_'.join([amber_atom3, amber_atom2, amber_atom4, amber_atom1])
            p5 = np.array([amber_coor3, amber_coor2, amber_coor4, amber_coor1])
            newkey6 = '_'.join([amber_atom2, amber_atom3, amber_atom4, amber_atom1])
            p6 = np.array([amber_coor2, amber_coor3, amber_coor4, amber_coor1])
            newkey11 = '_'.join(['X', amber_atom2, amber_atom4, amber_atom3])
            newkey21 = '_'.join(['X', amber_atom1, amber_atom4, amber_atom3])
            newkey31 = '_'.join(['X', amber_atom3, amber_atom4, amber_atom2])
            newkey41 = '_'.join(['X', amber_atom1, amber_atom4, amber_atom2])
            newkey51 = '_'.join(['X', amber_atom2, amber_atom4, amber_atom1])
            newkey61 = '_'.join(['X', amber_atom3, amber_atom4, amber_atom1])
            newkey12 = '_'.join(['X', 'X', amber_atom4, amber_atom3])
            newkey32 = '_'.join(['X', 'X', amber_atom4, amber_atom2])
            newkey52 = '_'.join(['X', 'X', amber_atom4, amber_atom1])
            barrier = 0; phase = 0; periodicity = 0; Etor = 0
            def check_key(somekey, somedic, somep, newdic):
                if somekey in somedic:
                    if somekey not in newdic:
                        newdic[somekey] = {}
                        newdic[somekey]['dihedral'] = []
                        newdic[somekey]['energy'] = []
                    degree = 0
                    degree = dihedral(somep)
                    barrier = float(somedic[somekey]['barrier'])
                    phase = float(somedic[somekey]['phase'])
                    periodicity = float(somedic[somekey]['periodicity'])
                    Etor = float(barrier)*(1+cos((float(periodicity)*degree-phase)/180*pi))
                    newdic[somekey]['dihedral'].append(degree)
                    newdic[somekey]['energy'].append(Etor)
                    return True
                elif somekey not in somedic:
                    return False
            ans1 = check_key(newkey1, out_of_plane, p1, out_of_plane_info)
            if ans1 == True:
                continue
            ans2 = check_key(newkey2, out_of_plane, p2, out_of_plane_info)
            if ans2 == True:
                continue
            ans3 = check_key(newkey3, out_of_plane, p3, out_of_plane_info)
            if ans3 == True:
                continue
            ans4 = check_key(newkey4, out_of_plane, p4, out_of_plane_info)
            if ans4 == True:
                continue
            ans5 = check_key(newkey5, out_of_plane, p5, out_of_plane_info)
            if ans5 == True:
                continue
            ans6 = check_key(newkey6, out_of_plane, p6, out_of_plane_info)
            if ans6 == True:
                continue
            ans7 = check_key(newkey11, out_of_plane, p1, out_of_plane_info)
            if ans7 == True:
                continue
            ans8 = check_key(newkey21, out_of_plane, p2, out_of_plane_info)
            if ans8 == True:
                continue
            ans9 = check_key(newkey31, out_of_plane, p3, out_of_plane_info)
            if ans9 == True:
                continue
            ans10 = check_key(newkey41, out_of_plane, p4, out_of_plane_info)
            if ans10 == True:
                continue
            ans11 = check_key(newkey51, out_of_plane, p5, out_of_plane_info)
            if ans11 == True:
                continue
            ans12 = check_key(newkey61, out_of_plane, p6, out_of_plane_info)
            if ans12 == True:
                continue
            ans13 = check_key(newkey12, out_of_plane, p1, out_of_plane_info)
            if ans13 == True:
                continue
            ans14 = check_key(newkey32, out_of_plane, p3, out_of_plane_info)
            if ans14 == True:
                continue
            ans15 = check_key(newkey52, out_of_plane, p5, out_of_plane_info)
            if ans15 == True:
                continue
    final_out_of_plane = {}
    for key in out_of_plane_info:
        final_out_of_plane[key] = 0
        final_out_of_plane[key] = sum(out_of_plane_info[key]['energy'])
    return final_out_of_plane

def out_of_plane_calc_ff14(pdbfile):
    resi_info = protein_info(pdbfile)[0]
    Bond_info = protein_info(pdbfile)[-1]
    resi_info_f = {}
    Bond_info_f = {}
    fp1 = open('Amber_ff14_parameter/out_of_plane_amber.json', 'r')
    out_of_plane = {}
    out_of_plane = json.load(fp1)
    fp2 = open('Amber_ff14_parameter/final_type.json', 'r')
    amber_atom_type = {}
    amber_atom_type = json.load(fp2)
    for i in range(len(resi_info)):
        for j in range(len(resi_info[i])):
            if resi_info[i][j][-1] not in resi_info_f:
                resi_info_f[resi_info[i][j][-1]] = {}
                resi_info_f[resi_info[i][j][-1]]['coor'] = [resi_info[i][j][1], resi_info[i][j][2], resi_info[i][j][3]]
                resi_info_f[resi_info[i][j][-1]]['KECSA_type'] = resi_info[i][j][0]
            elif resi_info[i][j][-1] in resi_info_f:
                print ('Error! Same atom exists in pdb! Please check: ', resi_info[i][j][-1])
    for i in range(len(Bond_info)):
        for j in range(len(Bond_info[i])):
            for k in range(len(Bond_info[i][j])):
                if Bond_info[i][j][k][0] not in Bond_info_f and Bond_info[i][j][k][1] not in Bond_info_f:
                    Bond_info_f[Bond_info[i][j][k][0]] = []
                    Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    Bond_info_f[Bond_info[i][j][k][1]] = []
                    Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                elif Bond_info[i][j][k][0] in Bond_info_f and Bond_info[i][j][k][1] not in Bond_info_f:
                    if Bond_info[i][j][k][1] not in Bond_info_f[Bond_info[i][j][k][0]]:
                        Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    Bond_info_f[Bond_info[i][j][k][1]] = []
                    Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                elif Bond_info[i][j][k][0] not in Bond_info_f and Bond_info[i][j][k][1] in Bond_info_f:
                    if Bond_info[i][j][k][0] not in Bond_info_f[Bond_info[i][j][k][1]]:
                        Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                    Bond_info_f[Bond_info[i][j][k][0]] = []
                    Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                elif Bond_info[i][j][k][0] in Bond_info_f and Bond_info[i][j][k][1] in Bond_info_f:
                    if Bond_info[i][j][k][1] not in Bond_info_f[Bond_info[i][j][k][0]]:
                        Bond_info_f[Bond_info[i][j][k][0]].append(Bond_info[i][j][k][1])
                    if Bond_info[i][j][k][0] not in Bond_info_f[Bond_info[i][j][k][1]]:
                        Bond_info_f[Bond_info[i][j][k][1]].append(Bond_info[i][j][k][0])
                else:
                    print ('Error! Please check bond: ', Bond_info[i][j][k][1], Bond_info[i][j][k][0])
    out_of_plane_info = {}
    for key in Bond_info_f:
        if len(Bond_info_f[key]) >4:
            print ('Error! Please check bonds connect to ', key)
        elif len(Bond_info_f[key]) == 3:
            if modify_name(resi_info_f[Bond_info_f[key][0]]['KECSA_type']) == 'HE2-HID':
                continue
            if modify_name(resi_info_f[Bond_info_f[key][1]]['KECSA_type']) == 'HE2-HID':
                continue
            if modify_name(resi_info_f[Bond_info_f[key][2]]['KECSA_type']) == 'HE2-HID':
                continue
            if modify_name(resi_info_f[key]['KECSA_type']) == 'HE2-HID':
                continue
            amber_atom1 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][0]]['KECSA_type'])]['amber_type']
            amber_atom2 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][1]]['KECSA_type'])]['amber_type']
            amber_atom3 = amber_atom_type[modify_name(resi_info_f[Bond_info_f[key][2]]['KECSA_type'])]['amber_type']
            amber_atom4 = amber_atom_type[modify_name(resi_info_f[key]['KECSA_type'])]['amber_type']
            amber_coor1 = resi_info_f[Bond_info_f[key][0]]['coor']
            amber_coor2 = resi_info_f[Bond_info_f[key][1]]['coor']
            amber_coor3 = resi_info_f[Bond_info_f[key][2]]['coor']
            amber_coor4 = resi_info_f[key]['coor']
            newkey1 = '_'.join([amber_atom1, amber_atom2, amber_atom4, amber_atom3])
            p1 = np.array([amber_coor1, amber_coor2, amber_coor4, amber_coor3])
            newkey2 = '_'.join([amber_atom2, amber_atom1, amber_atom4, amber_atom3])
            p2 = np.array([amber_coor2, amber_coor1, amber_coor4, amber_coor3])
            newkey3 = '_'.join([amber_atom1, amber_atom3, amber_atom4, amber_atom2])
            p3 = np.array([amber_coor1, amber_coor3, amber_coor4, amber_coor2])
            newkey4 = '_'.join([amber_atom3, amber_atom1, amber_atom4, amber_atom2])
            p4 = np.array([amber_coor3, amber_coor1, amber_coor4, amber_coor2])
            newkey5 = '_'.join([amber_atom3, amber_atom2, amber_atom4, amber_atom1])
            p5 = np.array([amber_coor3, amber_coor2, amber_coor4, amber_coor1])
            newkey6 = '_'.join([amber_atom2, amber_atom3, amber_atom4, amber_atom1])
            p6 = np.array([amber_coor2, amber_coor3, amber_coor4, amber_coor1])
            newkey11 = '_'.join(['X', amber_atom2, amber_atom4, amber_atom3])
            newkey21 = '_'.join(['X', amber_atom1, amber_atom4, amber_atom3])
            newkey31 = '_'.join(['X', amber_atom3, amber_atom4, amber_atom2])
            newkey41 = '_'.join(['X', amber_atom1, amber_atom4, amber_atom2])
            newkey51 = '_'.join(['X', amber_atom2, amber_atom4, amber_atom1])
            newkey61 = '_'.join(['X', amber_atom3, amber_atom4, amber_atom1])
            newkey12 = '_'.join(['X', 'X', amber_atom4, amber_atom3])
            newkey32 = '_'.join(['X', 'X', amber_atom4, amber_atom2])
            newkey52 = '_'.join(['X', 'X', amber_atom4, amber_atom1])
            barrier = 0; phase = 0; periodicity = 0; Etor = 0
            def check_key(somekey, somedic, somep, newdic):
                if somekey in somedic:
                    if somekey not in newdic:
                        newdic[somekey] = {}
                        newdic[somekey]['dihedral'] = []
                        newdic[somekey]['energy'] = []
                    degree = 0
                    degree = dihedral(somep)
                    barrier = float(somedic[somekey]['barrier'])
                    phase = float(somedic[somekey]['phase'])
                    periodicity = float(somedic[somekey]['periodicity'])
                    Etor = float(barrier)*(1+cos((float(periodicity)*degree-phase)/180*pi))
                    newdic[somekey]['dihedral'].append(degree)
                    newdic[somekey]['energy'].append(Etor)
                    return True
                elif somekey not in somedic:
                    return False
            ans1 = check_key(newkey1, out_of_plane, p1, out_of_plane_info)
            if ans1 == True:
                continue
            ans2 = check_key(newkey2, out_of_plane, p2, out_of_plane_info)
            if ans2 == True:
                continue
            ans3 = check_key(newkey3, out_of_plane, p3, out_of_plane_info)
            if ans3 == True:
                continue
            ans4 = check_key(newkey4, out_of_plane, p4, out_of_plane_info)
            if ans4 == True:
                continue
            ans5 = check_key(newkey5, out_of_plane, p5, out_of_plane_info)
            if ans5 == True:
                continue
            ans6 = check_key(newkey6, out_of_plane, p6, out_of_plane_info)
            if ans6 == True:
                continue
            ans7 = check_key(newkey11, out_of_plane, p1, out_of_plane_info)
            if ans7 == True:
                continue
            ans8 = check_key(newkey21, out_of_plane, p2, out_of_plane_info)
            if ans8 == True:
                continue
            ans9 = check_key(newkey31, out_of_plane, p3, out_of_plane_info)
            if ans9 == True:
                continue
            ans10 = check_key(newkey41, out_of_plane, p4, out_of_plane_info)
            if ans10 == True:
                continue
            ans11 = check_key(newkey51, out_of_plane, p5, out_of_plane_info)
            if ans11 == True:
                continue
            ans12 = check_key(newkey61, out_of_plane, p6, out_of_plane_info)
            if ans12 == True:
                continue
            ans13 = check_key(newkey12, out_of_plane, p1, out_of_plane_info)
            if ans13 == True:
                continue
            ans14 = check_key(newkey32, out_of_plane, p3, out_of_plane_info)
            if ans14 == True:
                continue
            ans15 = check_key(newkey52, out_of_plane, p5, out_of_plane_info)
            if ans15 == True:
                continue
    final_out_of_plane = {}
    for key in out_of_plane_info:
        final_out_of_plane[key] = 0
        final_out_of_plane[key] = sum(out_of_plane_info[key]['energy'])
    return final_out_of_plane


            
            
                
                
                
                    
        


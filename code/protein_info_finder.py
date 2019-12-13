# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:01:20 2017

@author: peijun
"""
from __future__ import division
from function import distance
from function import check
from function import distance
import json
import numpy as np

def torsion_info(residue, Torsion):
    fp1 = open('Amber_ff94_parameter/atom_num_type.json', 'r')
    dic1 = {}
    dic1 = json.load(fp1)
    resi_dic = {}
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            if residue[i][j][-1] not in resi_dic:
                resi_dic[residue[i][j][-1]] = {}
                resi_dic[residue[i][j][-1]]['KECSA2_type'] = residue[i][j][0]
                resi_dic[residue[i][j][-1]]['coor'] = [residue[i][j][1], residue[i][j][2], residue[i][j][3]]
            else:
                print ('same atom exists twice! please check:', residue[i][j])
    tor_dic = {}; tor_checklist = []; count = 0
    for i in range(len(Torsion)):
        for j in range(len(Torsion[i])):
            for k in range(len(Torsion[i][j])):
                count += 1
                newtor1 = '-'.join(Torsion[i][j][k][:4])
                newtor2 = '-'.join([Torsion[i][j][k][3], Torsion[i][j][k][2], Torsion[i][j][k][1], Torsion[i][j][k][0]])
                if newtor1 not in tor_checklist and newtor2 not in tor_checklist:
                    tor_checklist.append(newtor1)
                #else:
                    #print ('same torsion type exists twice! Please check:', Torsion[i][j][k])
    for tor in tor_checklist:
        [atom1, atom2, atom3, atom4] = tor.replace('-', ' ').split()
        atom1_num = atom1.replace('_', ' ').split()
        atom2_num = atom2.replace('_', ' ').split()
        atom3_num = atom3.replace('_', ' ').split()
        atom4_num = atom4.replace('_', ' ').split()
        atom1_amber_type = dic1[atom1_num[0]]
        atom2_amber_type = dic1[atom2_num[0]]
        atom3_amber_type = dic1[atom3_num[0]]
        atom4_amber_type = dic1[atom4_num[0]]
        newtor1 = '-'.join([atom1_amber_type, atom2_amber_type, atom3_amber_type, atom4_amber_type])
        newtor2 = '-'.join([atom4_amber_type, atom3_amber_type, atom2_amber_type, atom1_amber_type])
        p = []; degree = 0; dist = 0
        p.append(resi_dic[atom1]['coor'])
        p.append(resi_dic[atom2]['coor'])
        p.append(resi_dic[atom3]['coor'])
        p.append(resi_dic[atom4]['coor'])
        p = np.array(p)
        degree = dihedral(p)
        dist = distance(resi_dic[atom1]['coor'], resi_dic[atom4]['coor'])
        if newtor1 != newtor2:
            if newtor1 not in tor_dic and newtor2 not in tor_dic:
                tor_dic[newtor1] = {}
                tor_dic[newtor1]['dihedral'] = []
                tor_dic[newtor1]['dist'] = []
                tor_dic[newtor1]['1_4_info'] = []
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 in tor_dic and newtor2 not in tor_dic:
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 not in tor_dic and newtor2 in tor_dic:
                tor_dic[newtor2]['dihedral'].append(degree)
                tor_dic[newtor2]['dist'].append(dist)
                tor_dic[newtor2]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            else:
                print ('error! same torsion exists twice in tor_dic! Please check:', newtor1, newtor2)
        elif newtor1 == newtor2:
            if newtor1 not in tor_dic:
                tor_dic[newtor1] = {}
                tor_dic[newtor1]['dihedral'] = []
                tor_dic[newtor1]['dist'] = []
                tor_dic[newtor1]['1_4_info'] = []
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 in tor_dic:
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
    return tor_dic

def torsion_info_ff14(residue, Torsion):
    fp1 = open('Amber_ff14_parameter/ff14_num_nonb.json', 'r')
    dic1 = {}
    dic1 = json.load(fp1)
    resi_dic = {}
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            if residue[i][j][-1] not in resi_dic:
                resi_dic[residue[i][j][-1]] = {}
                resi_dic[residue[i][j][-1]]['KECSA2_type'] = residue[i][j][0]
                resi_dic[residue[i][j][-1]]['coor'] = [residue[i][j][1], residue[i][j][2], residue[i][j][3]]
            else:
                print ('same atom exists twice! please check:', residue[i][j])
    tor_dic = {}; tor_checklist = []; count = 0
    for i in range(len(Torsion)):
        for j in range(len(Torsion[i])):
            for k in range(len(Torsion[i][j])):
                count += 1
                newtor1 = '-'.join(Torsion[i][j][k][:4])
                newtor2 = '-'.join([Torsion[i][j][k][3], Torsion[i][j][k][2], Torsion[i][j][k][1], Torsion[i][j][k][0]])
                if newtor1 not in tor_checklist and newtor2 not in tor_checklist:
                    tor_checklist.append(newtor1)
                #else:
                    #print ('same torsion type exists twice! Please check:', Torsion[i][j][k])
    for tor in tor_checklist:
        [atom1, atom2, atom3, atom4] = tor.replace('-', ' ').split()
        atom1_num = atom1.replace('_', ' ').split()
        atom2_num = atom2.replace('_', ' ').split()
        atom3_num = atom3.replace('_', ' ').split()
        atom4_num = atom4.replace('_', ' ').split()
        atom1_amber_type = dic1[atom1_num[0]]['Amber_type']
        atom2_amber_type = dic1[atom2_num[0]]['Amber_type']
        atom3_amber_type = dic1[atom3_num[0]]['Amber_type']
        atom4_amber_type = dic1[atom4_num[0]]['Amber_type']
        newtor1 = '-'.join([atom1_amber_type, atom2_amber_type, atom3_amber_type, atom4_amber_type])
        newtor2 = '-'.join([atom4_amber_type, atom3_amber_type, atom2_amber_type, atom1_amber_type])
        p = []; degree = 0; dist = 0
        p.append(resi_dic[atom1]['coor'])
        p.append(resi_dic[atom2]['coor'])
        p.append(resi_dic[atom3]['coor'])
        p.append(resi_dic[atom4]['coor'])
        p = np.array(p)
        degree = dihedral(p)
        dist = distance(resi_dic[atom1]['coor'], resi_dic[atom4]['coor'])
        if newtor1 != newtor2:
            if newtor1 not in tor_dic and newtor2 not in tor_dic:
                tor_dic[newtor1] = {}
                tor_dic[newtor1]['dihedral'] = []
                tor_dic[newtor1]['dist'] = []
                tor_dic[newtor1]['1_4_info'] = []
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 in tor_dic and newtor2 not in tor_dic:
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 not in tor_dic and newtor2 in tor_dic:
                tor_dic[newtor2]['dihedral'].append(degree)
                tor_dic[newtor2]['dist'].append(dist)
                tor_dic[newtor2]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            else:
                print ('error! same torsion exists twice in tor_dic! Please check:', newtor1, newtor2)
        elif newtor1 == newtor2:
            if newtor1 not in tor_dic:
                tor_dic[newtor1] = {}
                tor_dic[newtor1]['dihedral'] = []
                tor_dic[newtor1]['dist'] = []
                tor_dic[newtor1]['1_4_info'] = []
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
            elif newtor1 in tor_dic:
                tor_dic[newtor1]['dihedral'].append(degree)
                tor_dic[newtor1]['dist'].append(dist)
                tor_dic[newtor1]['1_4_info'].append([atom1_num[0]+'_'+atom1_num[1], atom4_num[0]+'_'+atom4_num[1]])
    return tor_dic


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

    #if np.arctan2(y, x) < 0:
    #    return np.arctan2(y, x)+2*pi
    #else:
    #    return np.arctan2(y, x)
    return np.degrees(np.arctan2(y, x))





            

                        



        

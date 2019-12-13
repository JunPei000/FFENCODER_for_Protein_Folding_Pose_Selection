#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 16:47:35 2019

@author: peijun
"""
from __future__ import division
import json
from math import sqrt, exp, log, cos, pi

def tor_calc(torsion):
    fp1 = open('Amber_ff94_parameter/amber_num_info.json', 'r')
    dic_num_nonb = {}
    dic_num_nonb = json.load(fp1)
    fp2 = open('Amber_ff94_parameter/torsion_amber.json', 'r')
    dic_tor = {}
    dic_tor = json.load(fp2)
    fp3 = open('Amber_ff94_parameter/atom_num_KECSA.json', 'r')
    dic_num_KECSA = {}
    dic_num_KECSA = json.load(fp3)
    fp4 = open('Amber_ff94_parameter/ring_info.json', 'r')
    dic_ring_info = {}
    dic_ring_info = json.load(fp4)
    fp4 = open('Amber_ff94_parameter/six_ring_info.json', 'r')
    dic_6_ring_info = {}
    dic_6_ring_info = json.load(fp4)
    fp5 = open('Amber_ff94_parameter/final_type.json', 'r')
    dic_nonb_types = {}
    dic_nonb_types = json.load(fp5)
    dic_tor_energy = {}; energy_info = []
    for key_1 in torsion:
        newlist = [];barrier = []; divider = []; phase = []; periodicity = []
        newlist = key_1.split('-')
        #if newlist[0][0] == 'H' or newlist[1][0] == 'H' or newlist[2][0] == 'H' or newlist[3][0] == 'H':
        #    continue
        key_2 = '-'.join([newlist[3], newlist[2], newlist[1], newlist[0]])
        if key_1 not in dic_tor and key_2 not in dic_tor:
            key1_1 = '-'.join(['X', newlist[1], newlist[2], newlist[3]])
            key1_2 = '-'.join([newlist[0], newlist[1], newlist[2], 'X'])
            key1_3 = '-'.join(['X', newlist[2], newlist[1], newlist[0]])
            key1_4 = '-'.join([newlist[3], newlist[2], newlist[1], 'X'])
            if key1_1 not in dic_tor and key1_2 not in dic_tor and key1_3 not in dic_tor and key1_4 not in dic_tor:
                key2_1 = '-'.join(['X', newlist[1], newlist[2], 'X'])
                key2_2 = '-'.join(['X', newlist[2], newlist[1], 'X'])
                if key2_1 not in dic_tor and key2_2 not in dic_tor:
                    print ('Cannot find torsion! Please check :', key_1)
                elif key2_1 in dic_tor:
                    barrier = dic_tor[key2_1]['barrier']
                    divider = dic_tor[key2_1]['divider']
                    phase = dic_tor[key2_1]['phase']
                    periodicity = dic_tor[key2_1]['periodicity']
                elif key2_2 in dic_tor:
                    barrier = dic_tor[key2_2]['barrier']
                    divider = dic_tor[key2_2]['divider']
                    phase = dic_tor[key2_2]['phase']
                    periodicity = dic_tor[key2_2]['periodicity']
            elif key1_1 in dic_tor:
                barrier = dic_tor[key1_1]['barrier']
                divider = dic_tor[key1_1]['divider']
                phase = dic_tor[key1_1]['phase']
                periodicity = dic_tor[key1_1]['periodicity']
            elif key1_2 in dic_tor:
                barrier = dic_tor[key1_2]['barrier']
                divider = dic_tor[key1_2]['divider']
                phase = dic_tor[key1_2]['phase']
                periodicity = dic_tor[key1_2]['periodicity']
            elif key1_3 in dic_tor:
                barrier = dic_tor[key1_3]['barrier']
                divider = dic_tor[key1_3]['divider']
                phase = dic_tor[key1_3]['phase']
                periodicity = dic_tor[key1_3]['periodicity']
            elif key1_4 in dic_tor:
                barrier = dic_tor[key1_4]['barrier']
                divider = dic_tor[key1_4]['divider']
                phase = dic_tor[key1_4]['phase']
                periodicity = dic_tor[key1_4]['periodicity']
        elif key_1 in dic_tor:
            barrier = dic_tor[key_1]['barrier']
            divider = dic_tor[key_1]['divider']
            phase = dic_tor[key_1]['phase']
            periodicity = dic_tor[key_1]['periodicity']
        elif key_2 in dic_tor:
            barrier = dic_tor[key_2]['barrier']
            divider = dic_tor[key_2]['divider']
            phase = dic_tor[key_2]['phase']
            periodicity = dic_tor[key_2]['periodicity']
        if barrier == [] or divider == [] or phase == [] or periodicity == []:
            print ('Cannot find parameters for torsion! Please check: ', key_1)
        if len(torsion[key_1]['dist']) != len(torsion[key_1]['dihedral']) or len(torsion[key_1]['dist']) != len(torsion[key_1]['1_4_info']) or len(torsion[key_1]['1_4_info']) != len(torsion[key_1]['dihedral']):
            print ('Lost info for torsion! Please check: ', key_1)
        for i in range(len(torsion[key_1]['dist'])):
            dist = 0; dihedral = 0; charge_1 = 0; charge_4 = 0; Rmin_1 = 0; Rmin_4 = 4; epsilon_1 = 0; epsilon_4 = 0
            dist = round(float(torsion[key_1]['dist'][i]), 3)
            dihedral = round(float(torsion[key_1]['dihedral'][i]), 3)
            print 
            charge_1 = dic_num_nonb[torsion[key_1]['1_4_info'][i][0].split('_')[0]]['charge']
            charge_4 = dic_num_nonb[torsion[key_1]['1_4_info'][i][1].split('_')[0]]['charge']
            Rmin_1 = dic_num_nonb[torsion[key_1]['1_4_info'][i][0].split('_')[0]]['Rmin']
            Rmin_4 = dic_num_nonb[torsion[key_1]['1_4_info'][i][1].split('_')[0]]['Rmin']
            epsilon_1 = dic_num_nonb[torsion[key_1]['1_4_info'][i][0].split('_')[0]]['epsilon']
            epsilon_4 = dic_num_nonb[torsion[key_1]['1_4_info'][i][1].split('_')[0]]['epsilon']
            Etor = 0; Evdw = 0; Echarge = 0; Etotal = 0
            for j in range(len(barrier)):
                Etor += (float(barrier[j])/float(divider[j]))*(1+cos((float(periodicity[j])*dihedral-phase[j])/180*pi))
            Evdw = sqrt(epsilon_1*epsilon_4)*(((Rmin_1+Rmin_4)/dist)**12-2*((Rmin_1+Rmin_4)/dist)**6)
            Echarge = charge_1*charge_4*332.05/dist
            #print (key_1, barrier, divider, phase, periodicity, dist, dihedral, torsion[key_1], i, Etor, Evdw, Echarge)
            #Etotal = Etor+Evdw/2+Echarge/1.2
            KECSAtype_1 = dic_num_KECSA[torsion[key_1]['1_4_info'][i][0].split('_')[0]]
            KECSAtype_4 = dic_num_KECSA[torsion[key_1]['1_4_info'][i][1].split('_')[0]]
            newlist1 = []; newlist2 = []
            newlist1 = KECSAtype_1.split('-')
            newlist2 = KECSAtype_4.split('-')
            if newlist1[0] == 'OXT':
                newlist1[0] = 'O'
                continue
            if newlist2[0] == 'OXT':
                newlist2[0] = 'O'
                continue
            if len(newlist1[1]) == 4 and (newlist1[1][0] == 'N' or newlist1[1][0] == 'C'):
                continue
                newlist1[1] = newlist1[1][1:]
                if newlist1[0] == 'H1' or newlist1[0] == 'H2' or newlist1[0] == 'H3':
                    newlist1[0] = 'H'
                elif newlist1[0] == 'OXT':
                    newlist1[0] = 'O'
            if len(newlist2[1]) == 4 and (newlist2[1][0] == 'N' or newlist2[1][0] == 'C'):
                continue
                newlist2[1] = newlist2[1][1:]
                if newlist2[0] == 'H1' or newlist2[0] == 'H2' or newlist2[0] == 'H3':
                    newlist2[0] = 'H'
                elif newlist1[0] == 'OXT':
                    newlist1[0] = 'O'
            if newlist1[1] == 'HIE' or newlist1[1] == 'HIP':
                newlist1[1] = 'HID'
            if newlist2[1] == 'HIE' or newlist2[1] == 'HIP':
                newlist2[1] = 'HID'
            KECSAtype_1 = newlist1[0]+'-'+newlist1[1]
            KECSAtype_4 = newlist2[0]+'-'+newlist2[1]
            if KECSAtype_1 == 'HE2-HID' or KECSAtype_4 == 'HE2-HID':
                continue
            final_type1 = dic_nonb_types[KECSAtype_1]['final_type']
            final_type4 = dic_nonb_types[KECSAtype_4]['final_type']
            final_tor_type1 = '_'.join([final_type1, newlist[1], newlist[2], final_type4])
            final_tor_type2 = '_'.join([final_type1, newlist[2], newlist[1], final_type4])
            final_tor_type3 = '_'.join([final_type4, newlist[1], newlist[2], final_type1])
            final_tor_type4 = '_'.join([final_type4, newlist[2], newlist[1], final_type1])
            #print (key_1, barrier, divider, phase, periodicity, dist, dihedral, Etor, Evdw, Echarge, KECSAtype_1, KECSAtype_4)
            #print (torsion[key_1]['1_4_info'][i])
            for key in dic_ring_info:
                if torsion[key_1]['1_4_info'][i][0].split('_')[1] == torsion[key_1]['1_4_info'][i][1].split('_')[1]:
                    if torsion[key_1]['1_4_info'][i][0].split('_')[0] in dic_ring_info[key] and torsion[key_1]['1_4_info'][i][1].split('_')[0] in dic_ring_info[key]:
                        Etor = Etor
                        Evdw = 0
                        Echarge = 0
            for key in dic_6_ring_info:
                if torsion[key_1]['1_4_info'][i][0].split('_')[1] == torsion[key_1]['1_4_info'][i][1].split('_')[1]:
                    if torsion[key_1]['1_4_info'][i][0].split('_')[0] in dic_6_ring_info[key] and torsion[key_1]['1_4_info'][i][1].split('_')[0] in dic_6_ring_info[key]:
                        Etor = Etor
                        Evdw = Evdw/2
                        Echarge = Echarge/2
            Etotal = Etor+Evdw/2+Echarge/1.2
            energy_info.append([key_1, str(Etor), str(Evdw), str(Echarge), str(Etotal)])
            #newkey1 = KECSAtype_1+'__'+KECSAtype_4
            #newkey2 = KECSAtype_4+'__'+KECSAtype_1
            #if newkey1 != newkey2:
            #    if newkey1 not in dic_tor_energy and newkey2 not in dic_tor_energy:
            #        dic_tor_energy[newkey1] = []
            #        dic_tor_energy[newkey1].append(round(Etotal, 3))
            #    elif newkey1 in dic_tor_energy and newkey2 not in dic_tor_energy:
            #        dic_tor_energy[newkey1].append(round(Etotal, 3))
            #    elif newkey1 not in dic_tor_energy and newkey2 in dic_tor_energy:
            #        dic_tor_energy[newkey2].append(round(Etotal, 3))
            #    elif newkey1 in dic_tor_energy and newkey2 in dic_tor_energy:
            #        print ('Same torsion in dic_tor_energy! Please check : ', newkey1)
            #elif newkey1 == newkey2:
            #    if newkey1 not in dic_tor_energy:
            #        dic_tor_energy[newkey1] = []
            #        dic_tor_energy[newkey1].append(round(Etotal, 3))
            #    elif newkey1 in dic_tor_energy:
            #        dic_tor_energy[newkey1].append(round(Etotal, 3))
            #print (final_tor_type1, final_tor_type2, final_tor_type3, final_tor_type4)
            #Etotal = Etor+Evdw/2+Echarge/1.2
            if final_tor_type1 not in dic_tor_energy and final_tor_type2 not in dic_tor_energy and final_tor_type3 not in dic_tor_energy and final_tor_type4 not in dic_tor_energy:
                dic_tor_energy[final_tor_type1] = []
                dic_tor_energy[final_tor_type1].append(round(Etotal, 4))
            elif final_tor_type1 in dic_tor_energy:
                dic_tor_energy[final_tor_type1].append(round(Etotal, 4))
            elif final_tor_type2 in dic_tor_energy:
                dic_tor_energy[final_tor_type2].append(round(Etotal, 4))
            elif final_tor_type3 in dic_tor_energy:
                dic_tor_energy[final_tor_type3].append(round(Etotal, 4))
            elif final_tor_type4 in dic_tor_energy:
                dic_tor_energy[final_tor_type4].append(round(Etotal, 4))
            
    dic_tor_final = {}
    for key in dic_tor_energy:
        if key not in dic_tor_final:
            dic_tor_final[key] = round(sum(dic_tor_energy[key]), 4)
        elif key in dic_tor_final:
            print ('Same torsion in dic_tor_energy! Please check : ', key)
    return dic_tor_final

def nonb_calc(nonbond):
    fp1 = open('Amber_ff94_parameter/KECSA_num.json', 'r')
    dic_KECSA_num = {}
    dic_KECSA_num = json.load(fp1)
    fp2 = open('Amber_ff94_parameter/amber_num_info.json', 'r')
    dic_num_nonb = {}
    dic_num_nonb = json.load(fp2)
    fp3 = open('Amber_ff94_parameter/final_type.json', 'r')
    dic_nonb_types = {}
    dic_nonb_types = json.load(fp3)
    dic_nonb_energy = {}
    dic_nonb_final = {}
    dic_nonb_final_ambertype = {}
    energy_info = []
    for key in nonbond:
        newlist = []
        newlist = key.split('__')
        #if newlist[0][0] == 'H' or newlist[1][0] == 'H':
        #    continue
        dic_nonb_energy[key] = []
        dic_nonb_final[key] = 0
        charge_1 = 0; charge_2 = 0; Rmin_1 = 0; Rmin_2 = 0; epsilon_1 = 0; epsilon_2 = 0
        charge_1 = dic_num_nonb[dic_KECSA_num[newlist[0]]]['charge']
        charge_2 = dic_num_nonb[dic_KECSA_num[newlist[1]]]['charge']
        Rmin_1 = dic_num_nonb[dic_KECSA_num[newlist[0]]]['Rmin']
        Rmin_2 = dic_num_nonb[dic_KECSA_num[newlist[1]]]['Rmin']
        epsilon_1 = dic_num_nonb[dic_KECSA_num[newlist[0]]]['epsilon']
        epsilon_2 = dic_num_nonb[dic_KECSA_num[newlist[1]]]['epsilon']
        for i in range(len(nonbond[key])):
            dist = 0; Etotal = 0
            dist = round(nonbond[key][i], 3)
            Evdw, Echarge = E_calc(epsilon_1, epsilon_2, Rmin_1, Rmin_2, charge_1, charge_2, dist)
            Etotal = Evdw + Echarge
            #print key, epsilon_1, epsilon_2, Rmin_1, Rmin_2, charge_1, charge_2, dist, Evdw, Echarge
            #print key, dic_num_nonb[dic_KECSA_num[newlist[0]]]['Amber_typw'], dic_num_nonb[dic_KECSA_num[newlist[1]]]['Amber_typw'], Evdw, Echarge
            energy_info.append([str(Evdw), str(Echarge)])
            dic_nonb_energy[key].append(round(Etotal, 4))
    for key in dic_nonb_energy:
        newlist =[]
        newlist = key.split('__')
        newlist1 = []; newlist2 = []
        newlist1 = newlist[0].split('-')
        newlist2 = newlist[1].split('-')
        if newlist1[0] == 'OXT':
            continue
            newlist1[0] = 'O'
        if newlist2[0] == 'OXT':
            continue
            newlist2[0] = 'O'
        if len(newlist1[1]) == 4 and (newlist1[1][0] == 'N' or newlist1[1][0] == 'C'):
            continue
            newlist1[1] = newlist1[1][1:]
            if newlist1[0] == 'H1' or newlist1[0] == 'H2' or newlist1[0] == 'H3':
                newlist1[0] = 'H'
            elif newlist1[0] == 'OXT':
                newlist1[0] = 'O'
        if len(newlist2[1]) == 4 and (newlist2[1][0] == 'N' or newlist2[1][0] == 'C'):
            continue
            newlist2[1] = newlist2[1][1:]
            if newlist2[0] == 'H1' or newlist2[0] == 'H2' or newlist2[0] == 'H3':
                newlist2[0] = 'H'
            elif newlist1[0] == 'OXT':
                newlist1[0] = 'O'
        if newlist1[1] == 'HIE' or newlist1[1] == 'HIP':
            newlist1[1] = 'HID'
        if newlist2[1] == 'HIE' or newlist2[1] == 'HIP':
            newlist2[1] = 'HID'
        newlist[0] = newlist1[0]+'-'+newlist1[1]
        newlist[1] = newlist2[0]+'-'+newlist2[1]
        if newlist[0] == 'HE2-HID' or newlist[1] == 'HE2-HID':
            continue
        key1 = dic_nonb_types[newlist[0]]['final_type'] + '__' + dic_nonb_types[newlist[1]]['final_type']
        key2 = dic_nonb_types[newlist[1]]['final_type'] + '__' + dic_nonb_types[newlist[0]]['final_type']
        if key1 == key2:
            if key1 not in dic_nonb_final_ambertype:
                dic_nonb_final_ambertype[key1] = round(sum(dic_nonb_energy[key]), 4)
            elif key1 in dic_nonb_final_ambertype:
                dic_nonb_final_ambertype[key1] += round(sum(dic_nonb_energy[key]), 4)
        elif key1 != key2:
            if key1 not in dic_nonb_final_ambertype and key2 not in dic_nonb_final_ambertype:
                dic_nonb_final_ambertype[key1] = round(sum(dic_nonb_energy[key]), 4)
            elif key1 in dic_nonb_final_ambertype and key2 not in dic_nonb_final_ambertype:
                dic_nonb_final_ambertype[key1] += round(sum(dic_nonb_energy[key]), 4)
            elif key1 not in dic_nonb_final_ambertype and key2 in dic_nonb_final_ambertype:
                dic_nonb_final_ambertype[key2] += round(sum(dic_nonb_energy[key]), 4)
            elif key1 in dic_nonb_final_ambertype and key2 in dic_nonb_final_ambertype:
                print ('error! same nonbond pair exists twice!', key1, key2)
        dic_nonb_final[key] = round(sum(dic_nonb_energy[key]), 4)
    return dic_nonb_final_ambertype
            
            
            
def E_calc(epsilon_1, epsilon_4, Rmin_1, Rmin_4, charge_1, charge_4, dist):
    Evdw = 0; Echarge = 0
    Evdw = sqrt(epsilon_1*epsilon_4)*(((Rmin_1+Rmin_4)/dist)**12-2*((Rmin_1+Rmin_4)/dist)**6)
    Echarge = charge_1*charge_4*332.05/dist
    return Evdw, Echarge

                
            
            
            
        
            

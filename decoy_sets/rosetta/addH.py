import os
import time


direct1 = '/mnt/home/peijun/Documents/aminoacidstru/proten_folding/rosetta/structure_add_H/1acf/'
filenames = [i for i in os.listdir(direct1) if 'inserted_' in i and '.pdb' in i and '.swp' not in i and i.count('inserted') == 1]
print (len(filenames))
for filename in filenames:
    print (filename)
    os.chdir(direct1)
    input1 = filename+' '+filename.replace('.pdb', '_addH.pdb')
    os.system('/opt/software/schrodinger/suite_2018-4/utilities/prepwizard -fillsidechains -rehtreat -fix %s' % input1)
    while filename.replace('.pdb', '_addH.pdb') not in os.listdir(direct1):
        time.sleep(5)




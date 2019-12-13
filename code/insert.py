def insert(direct1,direct2,name):
    fp1 = open(direct1+name,'r').readlines()
    fp2 = open(direct2+'inserted_'+name,'w')
    for line in fp1:
        newlist = []
        newlist = line.split()
        if len(newlist) > 4 and 'ATOM' in line:
            if newlist[3]+'     '+newlist[4] in line:
                if len(newlist[4]) == 1:
                    line = line.replace(newlist[3]+'     '+newlist[4], newlist[3]+' T   '+newlist[4])
                elif len(newlist[4]) == 2:
	                line = line.replace(newlist[3]+'     '+newlist[4], newlist[3]+' T  '+newlist[4])
                elif len(newlist[4]) == 3:
	                line = line.replace(newlist[3]+'     '+newlist[4], newlist[3]+' T '+newlist[4])
            elif newlist[3]+'    '+newlist[4] in line:
	            if len(newlist[4]) == 1:
	                line = line.replace(newlist[3]+'    '+newlist[4], newlist[3]+' T   '+newlist[4])
	            elif len(newlist[4]) == 2:
	                line = line.replace(newlist[3]+'    '+newlist[4], newlist[3]+' T  '+newlist[4])
	            elif len(newlist[4]) == 3:
	                line = line.replace(newlist[3]+'    '+newlist[4], newlist[3]+' T '+newlist[4])
            elif newlist[3]+'   '+newlist[4] in line:
	            if len(newlist[4]) == 1:
	                line = line.replace(newlist[3]+'   '+newlist[4], newlist[3]+' T   '+newlist[4])
	            elif len(newlist[4]) == 2:
	                line = line.replace(newlist[3]+'   '+newlist[4], newlist[3]+' T  '+newlist[4])
	            elif len(newlist[4]) == 3:
	                line = line.replace(newlist[3]+'   '+newlist[4], newlist[3]+' T '+newlist[4])
            fp2.write(line)
    fp2.close()

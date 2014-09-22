import re, sys, os

RBX = ["HO2M12","HO2M10","HO2M08","HO2M06","HO2M04","HO2M02",
       "HO1M12","HO1M10","HO1M08","HO1M06","HO1M04","HO1M02",
       "HO001","HO002","HO003","HO004","HO005","HO006","HO007","HO008","HO009","HO010","HO011","HO012",
       "HO1P02","HO1P04","HO1P06","HO1P08","HO1P10","HO1P12",
       "HO2P02","HO2P04","HO2P06","HO2P08","HO2P10","HO2P12"]

ofileName ="HOPedestalSetting_new.xml"
os.system("rm "+ofileName)
ofile = open(ofileName,"w")

versionWritten = 0
ChangeTagName = 0
TagName = "NewPedTag"

for j in range (0,len(RBX)):
    name = RBX[j]+".xml"
    #print name
    
    thisfile = open(name,"r")
    lines = thisfile.readlines()
    thisfile.close()
    
    for line in lines:
        #print line
        if len(re.split("xml version",line))>1 and versionWritten==1:
            continue
            #print "will write ", line
        if ChangeTagName == 1:
            if len(re.split("CREATIONTAG",line))>1:
                print "Orig", line
                elements = re.split("\"string\"\>",line)
                tag = re.split("<",elements[1])
                newline = elements[0]+"\"string\">"+TagName+"<"+tag[1]
                print "New2", newline
                line = newline
        ofile.write(line)
        versionWritten=1
    ofile.write("\n")
    ofile.write("\n")

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sys


# In[2]:

print(sys.argv)
pedestal_table_name = sys.argv[1]
missing_channels_name = sys.argv[2]
output_file_name = sys.argv[3]

pedestal_table = open(pedestal_table_name)
missing_channels = open(missing_channels_name)


# In[3]:


output = open(output_file_name,"w")


# In[4]:


had1 = "#U ADC  << this is the unit\n"
had2 = "#             eta             phi             dep             det     cap0     cap1     cap2     cap3 widthcap0 widthcap1 widthcap2 widthcap3      DetId\n"
output.write(had1)
output.write(had2)


# In[5]:


#make it list
lines = []
for line in pedestal_table:
        v = line.split()
        lines.append(v)
        
        


# In[6]:


miss = []
for line in missing_channels:
        v = line.split()
        miss.append(v)


# In[7]:


miss_array=np.array(miss[2:])


# In[8]:


sur = np.zeros([len(miss_array),9,4])

pm1 = [-1,-1,-1,0,0,0,1,1,1]
pm2 = [-1,0,1,-1,0,1,-1,0,1] 

for ii in range(len(miss_array)):
    for jj in [-1,0,1]:
        for kk in [-1,0,1]:
            for ll in range(9):
                sur[ii,ll,0] = int(miss_array[ii,0])+pm1[ll]
                sur[ii,ll,1] = int(miss_array[ii,1])+pm2[ll]
                sur[ii,ll,2] = int(miss_array[ii,2])
                if miss_array[ii,3] == "HE":
                    sur[ii,ll,3] = 0
                if miss_array[ii,3] == "HB":
                    sur[ii,ll,3] = 1
                if miss_array[ii,3] == "HF":
                    sur[ii,ll,3] = 2
                if miss_array[ii,3] == "HO":
                    sur[ii,ll,3] = 3
                else:
                    sur[ii,ll,3] == -999


# In[9]:


sur_list = sur.astype(int).astype(str).tolist()


# In[10]:


for ii in range(len(sur_list)):
    for jj in range(len(sur_list[ii])):
        if (sur_list[ii][jj][3] != "0")and(sur_list[ii][jj][3] != "1")and(sur_list[ii][jj][3] != "2")and(sur_list[ii][jj][3] != "3"):
            sur_list[ii][jj][3] = "-999"
        if sur_list[ii][jj][3] == "0":
            sur_list[ii][jj][3] = "HE"
        if sur_list[ii][jj][3] == "1":
            sur_list[ii][jj][3] = "HB"
        if sur_list[ii][jj][3] == "2":
            sur_list[ii][jj][3] = "HF"
        if sur_list[ii][jj][3] == "3":
            sur_list[ii][jj][3] = "HO"


# In[11]:


miss_list = miss[2:]


# In[12]:


data_sur_list = []

for ii in range(len(lines)):
    for jj in range(len(sur_list)):
        for kk in range(len(sur_list[jj])):
            if lines[ii][0:4] == sur_list[jj][kk]:
                data_sur_list.append([np.array(lines[ii][4:-1]), sur_list[jj][4]] )


# In[13]:


sur_mean = []
for ii in range(len(miss_list)):
    sum_B0 = np.zeros(8)
    i = 0 
    for jj in range(len(data_sur_list)):

        #key_1 = float(mis[ii][0]) + np.sqrt(2)*float(mis[ii][1])+np.sqrt(2)*float(mis[ii][2])
        #key_2 = float(B[jj][1][0]) + np.sqrt(2)*float(B[jj][1][1])+np.sqrt(2)*float(B[jj][1][2])
        
        if miss_list[ii] == data_sur_list[jj][1]:
            sum_B0 = sum_B0 + data_sur_list[jj][0].astype(float)
            i = i + 1
            #print("-----")
            #print(mis[ii])
            #print(B[jj][1])
            #print("-----")

    
    if i != 0:
        sur_mean.append([(1/i)*sum_B0, miss_list[ii]])
    #if i == 0:
        
  
   
        
        
    


# In[14]:


add_line = []
for ii in range(len(sur_mean)):
    lines = ""
    if len(sur_mean[ii][1][0]) == 1:
        lines = lines + "                "
    if len(sur_mean[ii][1][0]) == 2:
        lines = lines + "               "
    if len(sur_mean[ii][1][0]) == 3:
        lines = lines + "              "
        
    lines = lines + str(sur_mean[ii][1][0])
    
    if len(sur_mean[ii][1][1]) == 1:
        lines = lines + "               "
    if len(sur_mean[ii][1][1]) == 2:
        lines = lines + "              "
    lines = lines + str(sur_mean[ii][1][1])
    
    lines = lines + "               "
    
    lines = lines + str(sur_mean[ii][1][2])
    
    lines = lines + "              "
    lines = lines + sur_mean[ii][1][3]
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][0]),".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][1]) ,".5f")   
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][2]), ".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][3]), ".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][4]),".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][5]),".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][6]), ".5f")
    lines = lines + "  "
    lines = lines + format(float(sur_mean[ii][0][7]), ".5f")
    lines = lines + "   "
    lines = lines + "XXXXXXXX"
    add_line.append(lines)
    


# In[15]:


for ii in range(len(add_line)):
    output.write(add_line[ii] + "\n")
    
output.close()


    


# In[ ]:





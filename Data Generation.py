#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import sys
import time
import pandas as pd


# In[10]:


startT=time.time();
f0=50;
max_rocof=2.5*1.1;
cost=np.array([[0.070203195,	0.157639984,	0.001658217],
[0.070203195,	0.157639984,	0.001658217],
[0.070203195,	0.157639984,	0.001658217],
[0.101700083,	0.168704813,	0.001236125],
[0.189613203,	0.157282398,	0.000747425],
[0.189613203,	0.157282398,	0.000747425],
[0.161812093,	0.142950917,	0.000662586],
[0.169909152,	0.176088611,	0.000429804],
[0.169909152,	0.176088611,	0.000429804],
[0.169909152,	0.176088611,	0.000429804],
[0.882502764,	0.168566554,	0.000124427]]) #quadradic cost function a*p^2+b*p+c for generators [c,b,a]

inertia=np.array([1.75,1.75,1.75,1.73,2.16,1.88,2.10,2.10,2.10,2.10,6.50]) #inertia of units in seconds

mbase=np.array([5.4,5.4,5.4,6.3,9.4,9.6,15.75,14.5,14.5,14.5,26.82]) #base power in MW

pArray=np.array([
[0,2.5,3,3.5], #i1
[0,2.5,3,3.5], #i2
[0,2.5,3,3.5], #i3
[0,3,3.5,4], #i4
[0,3.5,4,4.5,5,5.5,6,6.5], #i5
[0,3.5,4,4.5,5,5.5,6,6.5], #i6
[0,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5], #i7
[0,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5], #i8
[0,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5], #i9
[0,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5], #i10
#[0,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21] #i11
               ],dtype=object)

#i11p=np.arange(5,21.5,0.5)


# In[11]:


combos= np.zeros((4*4*4*4*8*8*11*11*11*11,11))
delta = np.zeros(11,dtype=bool)
costSum=np.zeros(4*4*4*4*8*8*11*11*11*11)
thermalMax=36; #maximum total thermal generation in MWs
thermalMin=16; #minimum total thermal generation in MWs
j=0
for i1 in range(4):
    for i2 in range(4):
        for i3 in range(4):
            for i4 in range(4):
                for i5 in range(8):
                      for i6 in range(8):
                        for i7 in range(11):
                            for i8 in range(11):
                                for i9 in range(11):
                                    for i10 in range(11):
                    #                     for i11 in range(32):
                                            reserve=0;
                                            delta[:]=0;
                                            combos[j,0]=pArray[0][i1];
                                            if(pArray[0][i1]>0):
                                                reserve=reserve+pArray[0][-1]-pArray[0][i1]
                                                delta[0]=1
                                            combos[j,1]=pArray[1][i2];
                                            if(pArray[1][i2]>0):
                                                reserve=reserve+pArray[1][-1]-pArray[1][i2]
                                                delta[1]=1
                                            combos[j,2]=pArray[2][i3];
                                            if(pArray[2][i3]>0):
                                                reserve=reserve+pArray[2][-1]-pArray[2][i3]
                                                delta[2]=1
                                            combos[j,3]=pArray[3][i4];
                                            if(pArray[3][i4]>0):
                                                reserve=reserve+pArray[3][-1]-pArray[3][i4]
                                                delta[3]=1
                                            combos[j,4]=pArray[4][i5];
                                            if(pArray[4][i5]>0):
                                                reserve=reserve+pArray[4][-1]-pArray[4][i5]
                                                delta[4]=1
                                            combos[j,5]=pArray[5][i6];
                                            if(pArray[5][i6]>0):
                                                reserve=reserve+pArray[5][-1]-pArray[5][i6]
                                                delta[5]=1
                                            combos[j,6]=pArray[6][i7];
                                            if(pArray[6][i7]>0):
                                                reserve=reserve+pArray[6][-1]-pArray[6][i7]
                                                delta[6]=1
                                            combos[j,7]=pArray[7][i8];
                                            if(pArray[7][i8]>0):
                                                reserve=reserve+pArray[7][-1]-pArray[7][i8]
                                                delta[7]=1
                                            combos[j,8]=pArray[8][i9];
                                            if(pArray[8][i9]>0):
                                                reserve=reserve+pArray[8][-1]-pArray[8][i9]
                                                delta[8]=1
                                            combos[j,9]=pArray[9][i10];
                                            if(pArray[9][i10]>0):
                                                reserve=reserve+pArray[9][-1]-pArray[9][i10]
                                                delta[9]=1
                                            #combos[j,10]=i11p[i11];
                                            if (combos[j].sum()>=thermalMax or combos[j].sum()<=thermalMin):
                                                continue
                                                 
                                            for counter in range(10): 
                                                 if (reserve < 0.9*pArray[counter][-1]*delta[counter]
                                                     or 2*(np.sum(mbase*inertia*delta)-mbase[counter]*inertia[counter]*delta[counter])<f0*combos[j,counter]*delta[counter]/max_rocof):
                                                        #print(2*(np.sum(mbase*inertia*delta)-mbase[counter]*inertia[counter]*delta[counter]))
                                                        #print(f0*pArray[counter][-1]*delta[counter]/max_rocof)
                                                        #print(combos[j,counter])
                                                        costSum[j]=0
                                                        j=j-1
                                                        break
                                            if (combos[j,counter]>0):
                                                     costSum[j]=costSum[j]+cost[counter,0]+combos[j,counter]*cost[counter,1]+combos[j,counter]*combos[j,counter]*cost[counter,2]
                                            j=j+1;
print(j)


# In[12]:


posCombos=combos[:j,:]
costSum=costSum[:j]
posCombos=np.c_[posCombos,np.transpose(posCombos.sum(axis=1)),np.transpose(costSum)]
posCombos=posCombos[posCombos[:, 11].argsort()]
print(len(posCombos));
demandSteps=(thermalMax-thermalMin)*2-1;

count=np.zeros(demandSteps,dtype=np.uintc)
for i in range(demandSteps):
    count[i] = np.count_nonzero(posCombos[:,11] == (thermalMin+0.5+i/2));
print(count)

np.set_printoptions(threshold=sys.maxsize)

cheaperCombos=np.zeros((500*demandSteps,11))
cheaperCombosRes=np.zeros((500*demandSteps,11))


# In[13]:


startAt=0;
j=0;
for i in range(demandSteps):
    cheaperCombos[j:j+np.minimum(500,count[i]),:]=posCombos[startAt:startAt+np.minimum(500,count[i]),0:11];
    startAt=startAt+count[i];
    j=j+np.minimum(500,count[i]);
    
cheaperCombos=cheaperCombos[:j,:];
    
for i in range(len(cheaperCombos)):
    for k in range(10):
        if cheaperCombos[i,k]==0:
            cheaperCombosRes[i,k]=0;
        else:
            cheaperCombosRes[i,k]=pArray[k][-1]-cheaperCombos[i,k];   

cheaperCombosRes=cheaperCombosRes[:j,:];
endT=time.time();
print(f"Runtime of the program is {endT - startT}")


# In[14]:


with pd.ExcelWriter('All_cheap_combos_laPalma.xlsx') as writer:  
    pd.DataFrame(cheaperCombos).to_excel(writer,sheet_name='generation',index=False,header=False)
    pd.DataFrame(cheaperCombosRes).to_excel(writer,sheet_name='reserve',index=False,header=False)


# In[15]:


print(startAt)


# In[18]:


(15570*30/21)/3600


# In[ ]:





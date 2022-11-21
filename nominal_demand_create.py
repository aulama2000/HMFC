# -*- coding: utf-8 -*-

import numpy as np
from numpy import genfromtxt
from scipy.optimize import fsolve
from Auxiliary import create_buildings, itemp, eqw

def aggnompower(temp,draw,resnumber,comnumber,whnumber, solar = np.zeros(24*6)):
    setpoint = [23.0,23.0,46,0]
    deadband = [1.5,1.5,-5.0] 
    w=np.transpose(np.stack((np.repeat(temp,6),solar)))
    
    def eq(u,t,intemp,Ad,Bd,Gd,COP):
        return (np.matmul(Ad,intemp[t-1,])+Bd*u*COP+np.matmul(Gd,w[t-1,]))[0]-setpoint[0]
    
    def f_function(delta,alpha,beta,alphaa):
        return delta/beta/(1+abs((alphaa-alpha)/alpha))
    
    intemp=np.zeros((len(w),3))
    wtemp=np.zeros(len(w))
    buildings = create_buildings(resnumber=resnumber, comnumber=comnumber, whnumber=whnumber)
    types=np.zeros(resnumber + comnumber + whnumber)
    nom_power = np.zeros([len(temp)*6,resnumber+comnumber+whnumber])
    
    
    for i in range(nom_power.shape[1]):
        for t in range(nom_power.shape[0]+1):
            if buildings[i].b_type == 'residential' or buildings[i].b_type == 'commercial':
                if t==0:
                    intemp[t,]=np.array([setpoint[0], setpoint[0], setpoint[0]])
                else:
                    nom_power[t-1,i] = fsolve(eq,10,args=(t,intemp,buildings[i].Ad,buildings[i].Bd,buildings[i].Gd,buildings[i].COP))
                    if t!=nom_power.shape[0]:
                        intemp[t,]=itemp(nom_power[t-1,i],intemp[t-1,],buildings[i].Ad,buildings[i].Bd,buildings[i].Gd,buildings[i].COP, w[t,])
                if buildings[i].b_type == 'residential':
                    types[i]=0
                else:
                    types[i]=1
            else:
                if t==0:
                    wtemp[t]=setpoint[2]
                else:
                    nom_power[t-1,i] = fsolve(eqw,5,args=(t,wtemp[t-1],buildings[i],draw*8.35/2.205))
                    if t!=nom_power.shape[0]:
                        wtemp[t,] = eqw(nom_power[t-1,i],t,wtemp[t-1],buildings[i], draw*8.35/2.205)+setpoint[2]
                types[i]=2
    
    nom_power[:,types==2]=nom_power[:,types==2]*4.5
    agg_nom_power = np.sum(nom_power, axis=1)
    


    CC=np.zeros(len(buildings))
    powers = np.zeros(len(buildings))
    
    
    for k in range(len(CC)):
        if buildings[k].b_type != 'water_heater':
            CC[k]=deadband[0]
            #power3[i] = float(CC[i]/(buildings[i].Bd[0]*buildings[i].COP))
            powers[k] = f_function(deadband[0], 0.99, buildings[k].Bd[0]*buildings[k].COP, 0.99)
            
        else:
            CC[k] = -deadband[2]
            #power3[i] = CC[i]/buildings[i].tao*buildings[i].beta*np.sum(buildings[i].w)+buildings[i].tao*buildings[i].gamma
            #power3[i] = CC[i]/(buildings[i].Prated*0.99/6)
            powers[k] = f_function(-deadband[2],0.99,buildings[k].Prated*0.99/6,0.99)
            
    power3=np.sum(powers)###capacity

    temp1=np.zeros((len(w),resnumber+comnumber+whnumber))
    temp2=np.zeros((len(w),resnumber+comnumber+whnumber))
    maxlevel=np.zeros(len(w))
    minlevel=np.zeros(len(w))
    for t in range(len(w)):
        for k in range(resnumber+comnumber+whnumber):
            temp1[t,k] = nom_power[t,k]/powers[k]
            temp2[t,k] = (buildings[k].Prated-nom_power[t,k])/powers[k]
        minlevel[t]=min(temp1[t,:])*power3
        maxlevel[t]=min(temp2[t,:])*power3
        maxlevel[maxlevel<0]=0
        minlevel[maxlevel<0]=0
        
    minlevel = agg_nom_power - minlevel
    maxlevel = agg_nom_power + maxlevel
    
    minlevel = np.mean(minlevel.reshape(-1,6),1)
    maxlevel = np.mean(maxlevel.reshape(-1,6),1)
    agg_nom_power = np.mean(agg_nom_power.reshape(-1,6),1)

            
    return agg_nom_power, minlevel, maxlevel, power3/6, buildings












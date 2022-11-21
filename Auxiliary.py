# -*- coding: utf-8 -*-

import numpy as np
from control.matlab import *
import random


random.seed(19)
np.random.seed(19)


def MFC(u, y, e, Kp, tsamp=600, ts=2, alpha=-0.00010):
    L = ts*tsamp
    sigma = np.arange(0, L+1, tsamp)
    term_y = (L-2*sigma)*y
    term_u = alpha*sigma*(L-sigma)*u
    to_be_integrated = -6*(term_y + term_u)/(L**3)
    Phi = np.trapz(to_be_integrated)
    return -(Phi+Kp*e)/alpha


class Building:
    def __init__(self, b_type='residential'):

        self.temp_list = list()
        self.u_list = list()

        #self.R1 = 1/200*np.random.uniform(0.8,1.2)
        self.COP = 4
        self.b_type = b_type
        if b_type == 'residential':
            self.K1 = 0.8*np.random.uniform(0.8,1.2)*0.5
            self.K2 = 5.4*np.random.uniform(0.8,1.2)*0.5
            self.K3 = 0.2*np.random.uniform(0.8,1.2)*0.5
            self.K4 = 1.5*np.random.uniform(0.8,1.2)*0.5
            self.K5 = 1.1*np.random.uniform(0.8,1.2)*0.5
            self.C1 = 45_000*np.random.uniform(0.8,1.2)*0.35
            self.C2 = 140_000*np.random.uniform(0.8,1.2)*0.35
            self.C3 = 30_000*np.random.uniform(0.8,1.2)*0.35
            self.solar_coeff = 0.15
            #self.C1 = Air_Density*self.Vol_Air*Cv
            self.Prated = 4.1
            
            #self.number = number
        elif b_type == 'commercial':
            self.K1 = 2.4*np.random.uniform(0.8,1.2)*0.5
            self.K2 = 16.2*np.random.uniform(0.8,1.2)*0.5
            self.K3 = 0.7*np.random.uniform(0.8,1.2)*0.5
            self.K4 = 4.5*np.random.uniform(0.8,1.2)*0.5
            self.K5 = 3.4*np.random.uniform(0.8,1.2)*0.5
            self.C1 = 140_000*np.random.uniform(0.8,1.2)*0.35
            self.C2 = 440_000*np.random.uniform(0.8,1.2)*0.35
            self.C3 = 100_000*np.random.uniform(0.8,1.2)*0.35
            self.solar_coeff = 0.15
            #self.C1 = Air_Density*self.Vol_Air*Cv
            self.Prated = 11.0
            #self.number = number
            # A = -1/R/C1; equivalent to 1/(R*C1)

        elif b_type == 'water_heater':
            self.Prated = 4.5
            self.UA = 0.002 #kW/C0.07*np.random.uniform(0.8,1.2)  #Btu/min
            self.Tamb = 23.0
            self.cp = 0.00116 #kWh/(kg.C)1
            self.Tfresh = 22#C70
            self.m=8.35*66/2.205#kg #lb
            self.setpoint = 46
            #self.mdot=draw*8.35/2.205 

        if b_type != 'water_heater':
            self.A = np.array([[-(self.K1+self.K2+self.K3+self.K5)/self.C1,(self.K1+self.K2)/self.C1,self.K5/self.C1],[(self.K1+self.K2)/self.C2, -(self.K1+self.K2)/self.C2, 0],[self.K5/self.C3, 0, -(self.K4+self.K5)/self.C3]])
            self.B = np.array([[1/self.C1], [0], [0]])
            self.G = np.array([[self.K3/self.C1, self.solar_coeff*1/self.C1], [0, self.solar_coeff*1/self.C2], [self.K4/self.C3, 0]])
            self.C = np.array([[1, 0, 0]])
            self.D = np.array([[0, 0, 0]])
            self.sysc = ss(self.A,np.concatenate((self.B,self.G), axis = 1),self.C,self.D)
            self.sysd = c2d(self.sysc,600)
            self.a, self.b, self.c, self.d = ssdata(self.sysd)
            self.Ad = np.squeeze(np.asarray(self.a))
            self.Bd = np.squeeze(np.asarray(self.b))[:,0]
            self.Gd = np.squeeze(np.asarray(self.b))[:,1:]
            self.Cd = np.squeeze(np.asarray(self.c))
            self.Dd = np.squeeze(np.asarray(self.d))[0]


def create_buildings(resnumber, comnumber, whnumber, buildings=[]):
    buildings = [Building(b_type='residential') for _ in range(resnumber)]
    buildings = buildings + [Building(b_type='commercial') for _ in range(comnumber)]
    buildings = buildings + [Building(b_type='water_heater') for _ in range(whnumber)]
    return buildings


def itemp(u, intemp, Ad, Bd, Gd, COP, w):
    return np.matmul(Ad, intemp) + Bd*u*COP+np.matmul(Gd, w)

def eqw(u,t,intemp,WH,mdraw):
    #deltaT = (u*WH.Prated*1000*3.412/6*0.99-WH.UA*10*(intemp-WH.Tamb)-WH.mdot[t-1]*WH.cp*(intemp-WH.Tfresh))/(WH.cp*WH.m)
    deltaT = (u*WH.Prated/6*0.99-WH.UA/6*(intemp-WH.Tamb)-mdraw[t-1]*WH.cp*(intemp-WH.Tfresh))/(WH.cp*WH.m)
    tmp1=intemp+deltaT-WH.setpoint
    return tmp1


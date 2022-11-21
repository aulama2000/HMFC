import time
from copy import deepcopy
from datetime import datetime, timedelta
import pandas as pd
import pytz
import numpy as np
from numpy import genfromtxt
from Auxiliary import MFC, itemp, eqw
from weather import acquireweather
from nominal_demand_create import aggnompower
from price_demand import Parameter, Solvers, optimize
import pickle



class Experiment:
    def __init__(self,buildings,power,dt,nom_power, non_hvac, nhnom,forecasted,water_draw):
        self.buildings = buildings
        self.pratings = np.zeros(len(buildings))
        self.setpoints = np.zeros(len(buildings))
        self.deadbands = np.zeros(len(buildings))
        for i in range(len(self.buildings)):
            self.pratings[i] = self.buildings[i].Prated
            if self.buildings[i].b_type=='residential' or self.buildings[i].b_type=='commercial':
                self.setpoints[i] = 23
                self.deadbands[i] = 1.5
            else:
                self.setpoints[i] = 46
                self.deadbands[i] = 5  
        self.dt = dt
        self.ts = 2
        self.u = np.zeros((len(buildings),int(24*60/self.dt)))
        self.temps = np.zeros((len(buildings),int(24*60/self.dt)*3))
        self.power = np.repeat(power,int(60/self.dt))*1000
        self.nhpower = np.repeat(non_hvac,int(60/self.dt))*1000
        self.np = np.repeat(nom_power,int(60/self.dt))*1000
        self.nhnp = np.repeat(nhnom,int(60/self.dt))*1000
        self.outdoortemp = np.zeros((int(24*60/self.dt),2))
        self.forecastedtemp = np.repeat(forecasted,int(60/self.dt))
        self.draw = water_draw


    def main(self, t, otemp):
        self.outdoortemp[t,0] = otemp
        if t < 4:
            pi = np.argsort((self.temps[:, 3*t]-self.setpoints)/self.deadbands)
        else:
            error = self.temps[:, 3*t-3]-self.setpoints
            u_tmp = deepcopy(self.u).astype(int)
            pi = np.argsort(MFC(u_tmp[:, t-self.ts-2:t-self.ts+1],
                                self.temps[:, [3*t-9, 3*t-6, 3*t-3]],
                                error,
                                3.8,
                                tsamp=self.dt*60,
                                ts=self.ts,
                                alpha=-0.00010))

        best = np.zeros(len(self.temps)+1)
        for on in range(len(self.temps)+1):
            best[on] = np.sum(self.pratings[pi[:on]])-self.power[t]
            if on > 0:
                if np.abs(best[on]) > np.abs(best[on-1]):
                    on = on-1
                    break
        
        turn_on = on
        self.u[:, t] = False
        self.u[pi[:turn_on], t] = True
        
        if t==0:            
            for i in range(len(self.temps)):
                if self.buildings[i].b_type != 'water_heater':
                    self.temps[i, 3*t:3*t+3] = itemp(self.u[i, t]*self.buildings[i].Prated,
                                        np.array([np.random.uniform(self.setpoints[i]-self.deadbands[i], self.setpoints[i]+self.deadbands[i]),np.random.uniform(self.setpoints[i]-self.deadbands[i], self.setpoints[i]+self.deadbands[i]),np.random.uniform(self.setpoints[i]-self.deadbands[i], self.setpoints[i]+self.deadbands[i])]),
                                        self.buildings[i].Ad,
                                        self.buildings[i].Bd,
                                        self.buildings[i].Gd,
                                        self.buildings[i].COP,
                                        self.outdoortemp[t,:])
                else:
                    self.temps[i,3*t]=eqw(self.u[i, t],1,np.random.uniform(self.setpoints[i]-self.deadbands[i], self.setpoints[i]+self.deadbands[i]),self.buildings[i],[self.draw[t]*8.35/2.205])+self.setpoints[i]
        else:
            for i in range(len(self.temps)):
                if self.buildings[i].b_type != 'water_heater':
                    self.temps[i, 3*t:3*t+3] = itemp(self.u[i, t]*self.buildings[i].Prated,
                                        self.temps[i, 3*t-3:3*t],
                                        self.buildings[i].Ad,
                                        self.buildings[i].Bd,
                                        self.buildings[i].Gd,
                                        self.buildings[i].COP,
                                        self.outdoortemp[t,:])
                else:
                    self.temps[i,3*t]=eqw(self.u[i, t],1,self.temps[i, 3*t-3],self.buildings[i],[self.draw[t]*8.35/2.205])+self.setpoints[i]
            
if __name__ == "__main__":
    tz = pytz.timezone('US/Eastern')
    number_of_loads = 300
    lalo=[(35.96,-83.92),(36.01,-84.27),(36.12,-83.49)]  #latitudes and longitudes of the locations  ## add more tuples for more locations
    resnumber = comnumber = whnumber = np.zeros(len(lalo), dtype=int)
    for n in range(len(lalo)):
        rn = np.random.uniform(size=len(lalo))
        rn = np.around((number_of_loads/np.sum(rn)*rn)).astype(int)  #create heteregenous loads adding up to 300
        resnumber[n] = rn[0]
        comnumber[n] = rn[1]
        whnumber[n] = rn[2]
    dt = 10
    draw = np.sum(genfromtxt('waterdraw.csv').reshape(-1,10),1)
    aggregators=[]
    marginal_cost = pd.read_excel('Cost_marginal.xlsx', index_col=0).to_numpy().flatten()
    alpha = pd.read_excel('alpha.xlsx', index_col=0).to_numpy()
    non_hvac = pd.read_excel('Non_HVAC_load.xlsx', index_col=0).to_numpy()[:len(lalo),:]
    exp = []
    time_steps = []
    newday = False
    
    while True:
        tm_stp = datetime.now(tz)
        tm_day = tm_stp.day
        tm_stp = int(tm_stp.hour*60/dt+tm_stp.minute//dt)  # get current time step
        if tm_stp == 0:
            newday = True
            time_steps = []
        elif tm_stp in time_steps:
            time.sleep(1*60)
            continue
        time_steps.append(tm_stp) 
        print(time_steps)
        print(tm_stp)
        try:
            temps,otemp = acquireweather(lalo)
        except:
            print('api error'+str(tm_stp))
        if newday:  #a new day            
            nom_powers = np.zeros(temps.shape)
            minlevels = np.zeros(temps.shape)
            maxlevels = np.zeros(temps.shape)
            capacities = np.zeros(temps.shape[1])
            for n in range (temps.shape[1]):
                nom_powers[:,n], minlevels[:,n], maxlevels[:,n],  capacities[n], tmp = aggnompower(temps[:,n],draw,resnumber[n],comnumber[n],whnumber[n]) 
                aggregators.append(tmp)
            with open(f'results{tm_day-1}.txt', "wb") as fp:
                pickle.dump(exp, fp)
            exp = []
            para = Parameter(nom_powers, maxlevels, minlevels, capacities, marginal_cost, alpha, non_hvac)
            solver = Solvers()
            results=optimize(para, solver)
            power = np.zeros((len(lalo), 24))
            non_hvac_power = np.zeros((len(lalo), 24))
            for n in range(power.shape[0]):   
                for t in range(power.shape[1]):
                    power[n,t] = results.hr[n,t].value
                    non_hvac_power[n,t] = results.dr[n,t].value                    
                exp.append(Experiment(aggregators[n], power[n,:],dt, nom_powers[:,n],non_hvac_power[n,:], para.Dd[n,:],temps[:,n],draw))
            newday=False
            time.sleep(10*60)
        for n in range(len(lalo)):
            try:
                exp[n].main(tm_stp,otemp[n])
            except:
                print('waiting midnight')        

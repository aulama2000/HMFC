

from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.core import *
import pyomo.environ as pyo
import matplotlib.pyplot as plt
from nominal_demand_create import aggnompower
import pandas as pd
import numpy as np
from weather import acquireweather
from numpy import genfromtxt


    
class Solvers():
    def __init__(self):
        """ define solver options """
        self.ipopt_path = "/Users/4ka/miniconda3/bin/ipopt"
        #self.ipopt_path = "C:/Users/4ka/Downloads/Ipopt-3.13.3-win64-msvs2019-md/Ipopt-3.13.3-win64-msvs2019-md/bin/ipopt"
        self.opt = pyo.SolverFactory("ipopt", executable = self.ipopt_path, solver_io='nl')
        

class Parameter():
    def __init__(self, nom_powers, maxlevels, minlevels, capacities, marginal_cost, alpha, non_hvac):
        """ define parameters and read data """
        self.Time_period = nom_powers.shape[0]  # define time period
        self.Aggregator_num = nom_powers.shape[1] # define aggregator number
        self.theta = 5 # define fluctuation penalty
        self.omega = 1 # define overall satisfaction penalty 
        self.beta = 1 # battery dissipation rate
        self.BigM = 1000  # big number for MIP 
        self.SmallN = 0.0000001 # small number for output 
        self.BI = 0 # initial virtual battery level 
        self.PX = 0.12 # production value 
        self.Pmax = 0.18  # upper bound of price, $/mwi 
        self.Bmin = -capacities/1000 # minimum virtual battery 
        self.Bmax = capacities/1000 # maximum virtual battery 
        self.BTL = 0 # min virtual battery, end time 
        self.BTU = 0 # max virtual battery, end time 
                      
        'electricity purchase price'
        self.Cost = marginal_cost
        self.C = self.Cost/1000
        
        'nominal hvac load mw'
        self.Dh = np.transpose(nom_powers)/1000

        'nominal non hvac load MW'
        self.Dd = non_hvac
        self.Dd = self.Dd*np.sum(self.Dh)/np.sum(self.Dd)
               
        'satisfaction parameter'
        self.alpha = alpha
        self.w = self.PX/self.alpha
      
        'max load level for non hvac'
        self.ED_max = self.Dd*1.2
        
        'min load level for non hvac'
        self.ED_min = self.Dd*0.8
        self.sum_ed_min = np.sum(self.ED_min)
 
        'max load level for hvac'
        #self.EH_max = pd.read_excel('EHmax_heter_load.xlsx', index_col=0).to_numpy()
        self.EH_max = np.transpose(maxlevels)/1000
        'min load level for hvac'
        #self.EH_min = pd.read_excel('EHmin.xlsx', index_col=0).to_numpy()
        self.EH_min = np.transpose(minlevels)/1000
        self.sum_eh_min = np.sum(self.EH_min)




def optimize(para, solver):
    Model = pyo.ConcreteModel()
    Model.Aggregators = pyo.RangeSet(0, para.Aggregator_num-1)
    Model.timeStep = pyo.RangeSet(0, para.Time_period-1)
    Model.Aggregators_timeStep = Model.Aggregators * Model.timeStep
    'define all the variables for the Model'
    #Model.obj_utility = pyo.Var()
    Model.obj_profit = pyo.Var()
    Model.obj_revenue = pyo.Var(Model.Aggregators_timeStep, domain = Reals) # revenue of dso
    Model.obj_cost = pyo.Var(Model.Aggregators_timeStep, domain = Reals) # cost of dso
    Model.obj_aggreg = pyo.Var(Model.Aggregators, domain = Reals) # objective value of aggregators
    Model.obj_pay = pyo.Var(Model.Aggregators_timeStep, domain = Reals) # payment of aggregators
    Model.obj_saty = pyo.Var(Model.Aggregators_timeStep, domain = Reals) # satisfaction valueof aggregators
    Model.obj_ave_price = pyo.Var(Model.Aggregators) # averaged price for aggregators
    Model.p = pyo.Var(Model.timeStep, bounds = (0,None)) # optimized price
    Model.dl = pyo.Var(Model.Aggregators_timeStep, domain = NonNegativeReals) # optimized total load
    Model.dr = pyo.Var(Model.Aggregators_timeStep, domain = NonNegativeReals) # optimized non hvac looad
    Model.hr = pyo.Var(Model.Aggregators_timeStep, domain = NonNegativeReals) # optimized hvac load
    Model.dAve = pyo.Var(bounds = ((1/para.Time_period)*(para.sum_ed_min+para.sum_eh_min),None))
    Model.ee = pyo.Var()
    Model.rho = pyo.Var() 
    Model.b = pyo.Var(Model.Aggregators_timeStep) # virtual battery status
    
    Model.constraints = ConstraintList()
    '============================================= constraints list========================================'
    'define objective for DSO'
    Model.obj_utility = pyo.Objective(expr = Model.obj_profit + para.omega * sum(Model.obj_saty[n,t] for n in Model.Aggregators for t in Model.timeStep) \
                                      - para.theta * para.Time_period * Model.ee, sense=maximize)
    
    
    Model.constraints.add(Model.obj_profit == sum(Model.obj_revenue[n,t] - Model.obj_cost[n,t] for n in Model.Aggregators for t in Model.timeStep))
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.obj_revenue[n,t] == Model.p[t] * Model.dl[n,t]) 
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.obj_cost[n,t] == para.C[t] * Model.dl[n,t])
    
    for t in range(para.Time_period):
        Model.constraints.add(float(para.C[t]) <= Model.p[t])
        
    for t in range(para.Time_period):
        Model.constraints.add(Model.p[t] <= para.Pmax)     
    
    Model.constraints.add(Model.dAve == (1/para.Time_period) * sum(Model.dl[n,t] for n in Model.Aggregators for t in Model.timeStep))
      
    Model.constraints.add(Model.ee/Model.dAve == Model.rho)
    
    for t in range(para.Time_period):
        Model.constraints.add(Model.ee >= sum(Model.dl[n,t] for n in Model.Aggregators))  
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(Model.obj_aggreg[n] == sum(Model.obj_saty[n,t] - Model.obj_pay[n,t] for t in Model.timeStep))
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.obj_saty[n,t] == (para.Dd[n,t] + para.Dh[n,t]) * para.w[n,t] * (Model.dl[n,t]/(para.Dd[n,t] + para.Dh[n,t]))**para.alpha[n,t])
             
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.obj_pay[n,t] == Model.p[t] * Model.dl[n,t])
             
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.dl[n,t] == Model.dr[n,t] + Model.hr[n,t])
            
    for n in range(para.Aggregator_num):
        Model.constraints.add(Model.obj_ave_price[n] == sum(Model.p[t] * Model.dl[n,t] for t in Model.timeStep)/sum(Model.dl[n,t] for t in Model.timeStep))
            
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.dl[n,t] == ((Model.p[t]/(para.w[n,t] * para.alpha[n,t]))**(1/(para.alpha[n,t]-1))) * (para.Dd[n,t] + para.Dh[n,t]))
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(sum(Model.dr[n,t] for t in Model.timeStep) >= 0.9 * sum(para.Dd[n,t] for t in Model.timeStep))
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(sum(Model.dr[n,t] for t in Model.timeStep) <= 1.1 * sum(para.Dd[n,t] for t in Model.timeStep))
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.dr[n,t] >= para.ED_min[n,t])
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.dr[n,t] <= para.ED_max[n,t])
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.hr[n,t] >= para.EH_min[n,t])
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.hr[n,t] <= para.EH_max[n,t])                
    
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.b[n,t] >= para.Bmin[n])
               
    for n in range(para.Aggregator_num):
        for t in range(para.Time_period):
            Model.constraints.add(Model.b[n,t] <= para.Bmax[n])    
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(Model.b[n,23] >= para.BTL)
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(Model.b[n,23] <= para.BTU)    
    
    for n in range(para.Aggregator_num):
        Model.constraints.add(Model.b[n,0] == para.beta * (para.BI + Model.hr[n,0] - para.Dh[n,0])) 
            
    for n in range(para.Aggregator_num):
        for t in range(1,para.Time_period):
            Model.constraints.add(Model.b[n,t] == para.beta * (Model.b[n,t-1] + Model.hr[n,t] - para.Dh[n,t])) 
    
     
    '============================================= solving process ========================================'
    instance = Model
    results = solver.opt.solve(instance, tee=False, #stream the solver output
                            keepfiles=False, #print the MILP file for examination
                            symbolic_solver_labels=False, # use human readable names
                            load_solutions=True) #  results moved to instance 
    return Model


if __name__ == "__main__":

    ###### This part retrieves the necessary infromation from nominal_demand_create and weather scripts  #######
    lalo=[(35.96,-83.92),(36.01,-84.27),(36.12,-83.49)]  #latitudes and longitudes of the locations  ## add more tuples for more locations
    #lalo=[(35.96,-63.92),(36.01,-84.27),(36.12,-103.49)]  #latitudes and longitudes of the locations  ## add more tuples for more locations
    temps,_=acquireweather(lalo)  # temperature forecasts of next 24 hours
    nom_powers = np.zeros(temps.shape)
    minlevels = np.zeros(temps.shape)
    maxlevels = np.zeros(temps.shape)
    capacities = np.zeros(temps.shape[1])
    draw = np.sum(genfromtxt('waterdraw.csv').reshape(-1,10),1)
    for n in range (temps.shape[1]):
        rn = np.random.uniform(size=3)
        rn = (300/np.sum(rn)*rn).astype(int) 
        print(rn)
        nom_powers[:,n], minlevels[:,n], maxlevels[:,n],  capacities[n], _ = aggnompower(temps[:,n],draw, resnumber=rn[0],comnumber=rn[1],whnumber=rn[2]) 
        #nom_powers[:,n], minlevels[:,n], maxlevels[:,n],  capacities[n], _ = aggnompower(temps[:,n],draw, resnumber=100,comnumber=100,whnumber=100) 

    marginal_cost = pd.read_excel('Cost_marginal.xlsx', index_col=0).to_numpy().flatten()
    alpha = pd.read_excel('alpha.xlsx', index_col=0).to_numpy()
    non_hvac = pd.read_excel('Non_HVAC_load.xlsx', index_col=0).to_numpy()[:len(lalo),:]
    ###########################################################################################################        
    
    para = Parameter(nom_powers, maxlevels, minlevels, capacities, marginal_cost, alpha, non_hvac)
    solver = Solvers()
    
    
    results=optimize(para, solver)
    
    
    #     for var in instance.component_data_objects(Var):
    #         print(str(var), var.value)
    #         for t in Model.time_range:
    #             if str(var) == 'pTin[%s]'%(t):
    #                 print(str(var), var.value)       
    #     for t in Model.time_range:
    #         a = Model.pTin[t].value
    #         print(a)
    #         print(Model.pTin[t].value)
    print('************** obj_utility **************')
    v = results.obj_utility
    print(v.value())
    
    print('************** p **************')
    p = []
    for t in results.timeStep:
        v = results.p[t]
        p.append(v.value)
        print(t, v.value)
        
    print('************** b **************')   
    optimized_load = np.zeros((para.Aggregator_num, para.Time_period))
    optimized_load2 = np.zeros((para.Aggregator_num, para.Time_period))
    optimized_load3 = np.zeros((para.Aggregator_num, para.Time_period))
    
    
    for n in results.Aggregators:   
        for t in results.timeStep:
            v = results.b[n,t]
            optimized_load[n,t] = results.hr[n,t].value
            optimized_load2[n,t] = results.dr[n,t].value
            optimized_load3[n,t] = results.dl[n,t].value
    
    
    
    
    plt.figure(1)
    plt.ylabel('resulted price $/MWh')
    plt.xlabel('hours')
    plt.plot(p)
    plt.title('Price')

    plt.show()     
     
    plt.plot(np.transpose(optimized_load), label='optimized')
    plt.gca().set_prop_cycle(None)
    plt.plot(nom_powers/1000, '--', label='nominal')
    plt.legend()
    plt.title('TCL loads')
    plt.ylabel('Load (MW)')
    plt.xlabel('Hour')
    plt.show()
    
    plt.plot(np.transpose(np.sum(optimized_load3,axis=0)))
    plt.plot(np.transpose(np.sum(para.Dd+para.Dh,axis=0)),'--')
    plt.title('TCL + non-TCL loads for all aggregators')
    plt.ylabel('Load (MW)')
    plt.xlabel('Hour')
    plt.show()


    
    







        

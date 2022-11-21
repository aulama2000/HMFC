# -*- coding: utf-8 -*-


import requests
import json
import numpy as np

def acquireweather (lalo):
    temp = np.zeros((24,len(lalo)))
    otemp = np.zeros(len(lalo))
    for n in range(len(lalo)):
        response = requests.get("https://api.openweathermap.org/data/2.5/onecall?lat="+ str(lalo[n][0])+"&lon="+str(lalo[n][1])+"&units=metric&appid=xxxxx")#replace xxxxx with your app id
        for t in range(len(temp)):
            temp[t,n] = json.loads(response.content)['hourly'][t]['temp']
        otemp[n] = json.loads(response.content)['current']['temp']
    return (temp), otemp




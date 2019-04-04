# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 21:50:29 2019

@author: Administrator
@ Data Analysis for Real-time data
@ Data Courtesy-- Sinopec

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
import itertools
import matplotlib.dates as mdates
from matplotlib import dates
from pylab import rcParams
import statsmodels.api as sm


'''
from IPython.display import HTML
video = open("animation.mp4", "rb").read()
video_encoded = video.encode("base64")
video_tag = '<video controls alt="test" src="data:video/x-m4v;base64,{0}">'.format(video_encoded)
HTML(video_tag)
'''



os.chdir('D:/Backup/Downloads/real-time data analysis')

df = pd.read_csv('sinopec.csv')

# check the NA and drop Na
df.head()
df.tail()
df.isnull().values.any()
df.isnull().values.sum()
# df.isnull()
df.dropna(axis = 'rows', inplace = True)

df['Time & Date'] = df['Time'] + ' ' + df['Date']

df['Time & Date'] = pd.to_datetime(df['Time & Date'])

# data wraggling for the calculation
Key_Parameter = ['Time & Date',\
             'Bit Depth',\
         'TORQUE', \
         'ROP', \
         'RPM', \
         'SPP', \
         'WOB', \
         'WOH', \
         'Pump 1', \
         'FLW OUT']

# we first analysis the out flow rate from the date file
df_KValue = df[Key_Parameter]
del df_KValue['Time & Date']
df_KValue = df_KValue.astype(float)
df_KValue.index = df['Time & Date']
df_KValue = df_KValue.resample('T').mean()
# data aggrate for the minutes domain

'''###
df_flo = df[['FLW OUT']]
df_flo.index = df['Time & Date']
# check the value of every element
df_flo.values
df_flo = df_flo.astype(int)
y = df_flo['FLW OUT'].resample('T').mean()
# data visulization for the calculation
plt.style.use('ggplot')
y.plot()
'''###

f,axes=plt.subplots(9,1)

axes[0].plot((df_KValue['Bit Depth']), 'g-')
axes[1].plot((df_KValue['TORQUE']), 'g-')
axes[2].plot((df_KValue['ROP']), 'r-')
axes[3].plot((df_KValue['RPM']), 'r-')
axes[4].plot((df_KValue['SPP']), 'y-')
axes[5].plot((df_KValue['WOB']), 'y-')
axes[6].plot((df_KValue['Pump 1']), 'm-')
axes[7].plot((df_KValue['WOH']), 'm-')
axes[8].plot((df_KValue['FLW OUT']), 'b-')


axes[0].set_ylabel("Bit Depth", rotation=90)
axes[0].xaxis.set_major_locator(mdates.HourLocator(interval=3))
axes[0].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
axes[1].set_ylabel("Torque", rotation=90)
axes[2].set_ylabel("ROP", rotation=90)
axes[3].set_ylabel("RPM", rotation=90)
axes[4].set_ylabel("SPP", rotation=90)
axes[5].set_ylabel("WOB", rotation=90)
axes[6].set_ylabel("Pump 1", rotation=90)
axes[7].set_ylabel("WOH", rotation=90)
axes[8].set_ylabel("Flow out", rotation=90)

hour_fmt = dates.DateFormatter('%H:%M')

for i in range(9):
    axes[i].xaxis.set_major_formatter(hour_fmt)
    axes[i].xaxis.set_major_locator(mdates.HourLocator(interval=3))
    axes[i].xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    axes[i].set_xlim(pd.Timestamp('03-28-2019  8:00:10'), pd.Timestamp('03-29-2019  5:10:10'))
    axes[i].spines['top'].set_visible(False)
    axes[i].xaxis.set_visible(False)
# plt.tight_layout()

plt.xlabel('real-time lost circulation data March-29th-2019 (source: Sinopec)')

axes[0].fill_between(df_KValue.index,
                df_KValue.iloc[:, 0],
                df_KValue.iloc[:, 1], color='k', alpha=.2)

# invisible the figure in the matplot
axes[1].set_visible(False)
axes[0].xaxis.set_visible(False)
axes[0].spines['right'].set_visible(False)


# statistic analysis pump pressure response!
y = df_KValue['SPP']
from statsmodels.tsa.seasonal import seasonal_decompose
# decomposition = sm.tsa.seasonal_decompose(y, model='additive')
decomposition = seasonal_decompose(y.values, freq=168)
fig = decomposition.plot()
plt.show()


# ARIMA real-time machine learning model for the calculation
p = d = q = range(0, 2)
pdq = list(itertools.product(p, d, q))
seasonal_pdq = [(x[0], x[1], x[2], 12) for x in list(itertools.product(p, d, q))]
list(itertools.product(p, d, q))
print('Examples of parameter combinations for Seasonal ARIMA...')
print('SARIMAX: {} x {}'.format(pdq[1], seasonal_pdq[1]))
print('SARIMAX: {} x {}'.format(pdq[1], seasonal_pdq[2]))
print('SARIMAX: {} x {}'.format(pdq[2], seasonal_pdq[3]))
print('SARIMAX: {} x {}'.format(pdq[2], seasonal_pdq[4]))

for param in pdq:
    for param_seasonal in seasonal_pdq:
        try:
            mod = sm.tsa.statespace.SARIMAX(y,
                                            order=param,
                                            seasonal_order=param_seasonal,
                                            enforce_stationarity=False,
                                            enforce_invertibility=False)
            results = mod.fit()
            print('ARIMA{}x{}12 - AIC:{}'.format(param, param_seasonal, results.aic))
        except:
            continue
        
print(results.summary().tables[1]) 


'''
axes[1].plot((df_KValue['Time & Date']),(df_KValue['TORQUE']), 'g-')
axes[2].plot((df_KValue['Time & Date']),(df_KValue['ROP']), 'k-')
axes[3].plot((df_KValue['Time & Date']),(df_KValue['RPM']), 'k-')
axes[4].plot((df_KValue['Time & Date']),(df_KValue['SPP']), 'y-')
axes[5].plot((df_KValue['Time & Date']),(df_KValue['WOB']), 'y-')
axes[6].plot((df_KValue['Time & Date']),(df_KValue['Pump 1']), 'm-')
axes[7].plot((df_KValue['Time & Date']),(df_KValue['WOH']), 'm-')
axes[8].plot((df_KValue['Time & Date']),(df_KValue['FLW OUT']), 'b-')

'''

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 12:52:27 2022

@author: jstout
"""

import pandas as pd

dframe = pd.read_csv('conversion_table.csv', 
                     skiprows=1,)

def rem_trail_space(cell):
    if type(cell) is not str:
        return cell
    elif len(cell)==0:
        return cell
    elif cell[-1]==' ':
        return cell[:-1]
    else:
        return cell
    

for idx in dframe.keys():
    dframe[idx]=dframe[idx].apply(rem_trail_space)
    print(idx)
dframe.set_index('system', inplace=True)  

conv_dict=dframe.to_dict('index')

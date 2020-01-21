#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:27:50 2019

@author: rplazevedo
"""

#import matplotlib.pyplot as plt
import configparser
import reduce_data
import cleanup


config = configparser.ConfigParser()
config.read('input.ini')
all_phi = config['Parameters']['all_phi']

if input('Run cleanup? y/(n)\n') in ['y','Y']:
    cleanup.run()

run_init = input('Run initial sPRS? (y)/n \n')
run_sPRS_cont = input('Run full sPRS? (y)/n \n')
run_mPRS = input('Run mPRS? (y)/n \n')
plots = input('Plot data? y/(n) \n')

import sprs_init
import sprs_cont
import meshless_prs
import merge_v
import merge_phi
import plot_data

if run_init in ['','y','Y']:
    print('\n---Running initial sPRS---' )
    sprs_init.run()

if run_sPRS_cont in ['','y','Y']:
    print('\n---Running rest of sPRS---' )
    sprs_cont.run()
    merge_v.run()
    reduce_data.run()
    if all_phi in ['yes', 'YES', 'some', 'SOME']:
        merge_phi.run()

if run_mPRS in ['','y','Y']:
    print('\n---Running mPRS---')
    meshless_prs.run()
    merge_v.run()
    reduce_data.run()
    if all_phi in ['yes', 'YES', 'some', 'SOME']:
        merge_phi.run()
        
if plots in ['y','Y']:
    plot_data.run()

print('---Done!---')

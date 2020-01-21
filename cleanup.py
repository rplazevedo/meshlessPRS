#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 17:46:49 2019

@author: rplazevedo
"""

import os
import configparser

def run():    
    config = configparser.ConfigParser()
    config.read('input.ini')
    parts = int(config['PRS']['sPRS_parts'])+int(config['PRS']['mPRS_parts'])
    name = str(config['Parameters']['name'])
    
    if input('Delete all data? (y)/n \n') in ['', 'y','Y']:
        print('\n---Deleting all data files---')
        for p in range(1,parts+1):
            try:
                os.remove(str(name) + '_vdata' +str(p)+'.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_nwalls_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_exc_pts_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_tspan_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_nop_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_phi_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_d_phi_data' + str(p) + '.npy')
            except FileNotFoundError:
                pass
    
        try:         
            os.remove(str(name) + '_v.mymemmap')
        except FileNotFoundError:
            pass
        try:
            os.remove(str(name) + '_nwalls.mymemmap')
        except FileNotFoundError:
            pass
        try:         
            os.remove(str(name) + '_exc_pts.mymemmap')
        except FileNotFoundError:
            pass
        try:
            os.remove(str(name) + '_nop.mymemmap')
        except FileNotFoundError:
            pass
        
         
        number = 0
        while os.path.isfile(str(name) +'_'+str(number)+ '_v.txt'):
            os.remove(str(name) +'_'+str(number)+ '_v.txt')
            os.remove(str(name) +'_'+str(number)+ '_nwalls.txt')
            os.remove(str(name) +'_'+str(number)+ '_exc_pts.txt')
            try:
                os.remove(str(name) +'_'+str(number)+ '_nop.txt')
            except FileNotFoundError:
                pass
            number+=1    
    
        number = 0
        while os.path.isfile(str(name) +'_'+str(number)+ '_v.data'):
            os.remove(str(name) +'_'+str(number)+ '_v.data')
            os.remove(str(name) +'_'+str(number)+ '_nwalls.data')
            os.remove(str(name) +'_'+str(number)+ '_exc_pts.data')
            try:
                os.remove(str(name) +'_'+str(number)+ '_nop.data')
            except FileNotFoundError:
                pass
            number+=1
        try:
            os.remove('out.txt') 
        except FileNotFoundError:
            pass 
    
    else:
        delete_bulk = input('Delete bulk data? (y)/n \n')
        delete_mem_maps = input('Delete memory maps? (y)/n \n')
        delete_merged = input('Delete merged data? (y)/n \n')
    
        if delete_bulk in ['', 'y','Y']:
            for p in range(1,parts+1):
                try:
                    os.remove(str(name) + '_vdata' +str(p)+'.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_nwalls_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_tspan_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_exc_pts_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_nop_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_phi_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(str(name) + '_d_phi_data' + str(p) + '.npy')
                except FileNotFoundError:
                    pass
                    
        if delete_mem_maps in ['', 'y','Y']:
            try:         
                os.remove(str(name) + '_v.mymemmap')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_nwalls.mymemmap')
            except FileNotFoundError:
                pass
            try:         
                os.remove(str(name) + '_exc_pts.mymemmap')
            except FileNotFoundError:
                pass
            try:
                os.remove(str(name) + '_nop.mymemmap')
            except FileNotFoundError:
                pass
            
        if delete_merged in ['', 'y','Y']:
            number = 0
            while os.path.isfile(str(name) +'_'+str(number)+ '_v.txt'):
                os.remove(str(name) +'_'+str(number)+ '_v.txt')
                os.remove(str(name) +'_'+str(number)+ '_nwalls.txt')
                os.remove(str(name) +'_'+str(number)+ '_exc_pts.txt')
                try:
                    os.remove(str(name) +'_'+str(number)+ '_nop.txt')
                except FileNotFoundError:
                    pass
                number+=1
            number = 0
            while os.path.isfile(str(name) +'_'+str(number)+ '_v.data'):
                os.remove(str(name) +'_'+str(number)+ '_v.data')
                os.remove(str(name) +'_'+str(number)+ '_nwalls.data')
                os.remove(str(name) +'_'+str(number)+ '_exc_pts.data')
                try:
                    os.remove(str(name) +'_'+str(number)+ '_nop.data')
                except FileNotFoundError:
                    pass
                number+=1
           
    print('---Cleanup done!---')

# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 14:31:44 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

import os
#import shutil
#classes_dir = r'C:\Users\jlbus\OneDrive\Documents\mann\code\classes'
#classes_dir = r'/mnt/c/Users/jlbus/OneDrive/Documents/mann/code/classes' 
#os.chdir(classes_dir)   
from rotations import myRotations
from mastQuery import mastObj
from eleanorQuery import eleanorObj
from tesscutQuery import tesscutObj
from k2Query import k2obj

from datetime import date
#from functions import * 

class queryAll:
    def __init__(self,query_df,download_dir,lc_list = ['spoc','eleanor','cpm','k2','k2sff'], key_col = 'tic'):
        self.download_dir = download_dir
        self.query_df = query_df
        self.key_col = key_col
        self.lc_list = lc_list
        ### initialize storage locations and download storage diagnostic lists
        rot_result_cols = ['sector','LS_Per1','LS_Per2','LS_Per3','LS_Power1','LS_Power2','LS_Power3','ac_period']
        self.pdc_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.sap_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.ele_corr_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.ele_pca_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.cpm_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.k2_pdc_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.k2_sap_rots_df = pd.DataFrame(columns = rot_result_cols)
        self.k2sff_cor_rots_df = pd.DataFrame(columns = rot_result_cols)
        
        self.pdc_rot_fn = os.path.join(download_dir,'pdc_rots.csv')
        self.sap_rot_fn = os.path.join(download_dir,'sap_rots.csv')
        self.ele_corr_rot_fn = os.path.join(download_dir,'ele_corr_rots.csv')
        self.ele_pca_rot_fn = os.path.join(download_dir,'ele_pca_rots.csv')
        self.cpm_rot_fn = os.path.join(download_dir,'cpm_rots.csv')
        self.k2_pdc_rot_fn = os.path.join(download_dir,'k2_pdc_rots.csv')
        self.k2_sap_rot_fn = os.path.join(download_dir,'k2_sap_rots.csv')
        self.k2sff_cor_rot_fn = os.path.join(download_dir,'k2sff_cor_rots.csv')
        
        self.spoc_dir = os.path.join(download_dir,'spoc_LCs')
        self.ele_dir = os.path.join(download_dir,'eleanor_LCs')
        self.cpm_dir = os.path.join(download_dir,'cpm_LCs')
        self.k2_dir = os.path.join(download_dir,'k2_LCs')
        self.k2sff_dir = os.path.join(download_dir,'k2sff_LCs')
        
        self.mast_flux_type = ['PDCSAP_flux','SAP_flux']
        self.ele_flux_type = ['corr','pca']
        self.k2sff_flux_type = ['fcor']
        ### download folders and availability results columns
        results_cols = [key_col]
        if 'spoc' in lc_list:
            results_cols.append('spoc_avail')
            if os.path.exists(self.spoc_dir) == False:
                os.mkdir(self.spoc_dir)
        if 'eleanor' in lc_list:
            results_cols.append('ele_avail')
            if os.path.exists(self.ele_dir) == False:
                os.mkdir(self.ele_dir)
        if 'cpm' in lc_list:
            results_cols.append('tesscut_avail')
            if os.path.exists(self.cpm_dir) == False:
                os.mkdir(self.cpm_dir) 
        if 'k2' in lc_list:
            results_cols.append('k2_avail')
            if os.path.exists(self.k2_dir) == False:
                os.mkdir(self.k2_dir)
        if 'k2sff' in lc_list:
            results_cols.append('k2sff_avail')
            if os.path.exists(self.k2sff_dir) == False:
                os.mkdir(self.k2sff_dir)
        ### best period/power columns
        if 'spoc' in lc_list: results_cols = np.append(results_cols, ['pdc_per','pdc_power','sap_per','sap_power'])
        if 'eleanor' in lc_list: results_cols = np.append(results_cols, ['ele_corr_per','ele_corr_power','ele_corr_ac_per','ele_pca_per','ele_pca_power','ele_pca_ac_per'])
        if 'cpm' in lc_list: results_cols = np.append(results_cols, ['cpm_per','cpm_power','cpm_ac_per'])
        if 'k2' in lc_list: results_cols = np.append(results_cols, ['k2_pdc_per','k2_pdc_power','k2_pdc_ac_per','k2_sap_per','k2_sap_power','k2_sap_ac_per'])
        if 'k2sff' in lc_list: results_cols = np.append(results_cols, ['k2sff_cor_per','k2sff_cor_power','k2sff_ac_per'])
        
        self.results_cols = results_cols
        self.results_df = pd.DataFrame(columns = results_cols) 
        self.results_fn = os.path.join(download_dir,'results_' + str(date.today()) + '.csv')
        self.spoc_cols = ['spoc_avail','pdc_per','pdc_power','pdc_ac_per','sap_per','sap_power','sap_ac_per']
        self.ele_cols = ['ele_avail','ele_corr_per','ele_corr_power','ele_corr_ac_per','ele_pca_per','ele_pca_power','ele_pca_ac_per']
        self.cpm_cols = ['tesscut_avail','cpm_per','cpm_power','cpm_ac_per']
        self.k2_cols = ['k2_avail','k2_pdc_per','k2_pdc_power','k2_pdc_ac_per','k2_sap_per','k2_sap_power','k2_sap_ac_per']
        self.k2sff_cols = ['k2sff_avail','k2sff_cor_per','k2sff_cor_power','k2sff_ac_per']
        
    def get_the_dat(self):
        for i,row in self.query_df.iterrows():
            if str(row[self.key_col]) != 'nan':
                temp_results = pd.DataFrame(columns = self.results_cols)
                temp_results[self.key_col] = [row[self.key_col]]
                if ('spoc' in self.lc_list) & (str(row['tic']) != 'nan'): 
                    spoc_res_df = self.get_spoc(row=row,i=i)
                    for col in self.spoc_cols: temp_results[col] = spoc_res_df[col]
                if 'eleanor' in self.lc_list: 
                    ele_res_df = self.get_eleanor(row=row,i=i)
                    for col in self.ele_cols: temp_results[col] = ele_res_df[col]
                if 'cpm' in self.lc_list: 
                    cpm_res_df = self.get_cpm(row=row,i=i)
                    for col in self.cpm_cols: temp_results[col] = cpm_res_df[col]
                if 'k2' in self.lc_list: 
                    k2_res_df = self.get_k2(row=row,i=i,k2sff = False)
                    for col in self.k2_cols: temp_results[col] = k2_res_df[col]
                if 'k2sff' in self.lc_list: 
                    k2sff_res_df = self.get_k2sff(row=row,i=i,k2sff = True)
                    for col in self.k2sff_cols: temp_results[col] = k2sff_res_df[col]
                self.results_df = pd.concat([self.results_df,temp_results])
                self.results_df.to_csv(self.results_fn,index=False)
            else:
                print("Object " + str(i+1) + " has no key column identifier. Moving on...")

    def get_spoc(self,row,i):
        print("Working on SPOC LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
        mast_obj = mastObj(tic = row['tic'], ra = row['ra'], dec = row['dec'], download_dir = self.spoc_dir, products = ['LC'], rotations = False)
        mast_obj.download(keep_fits = False)
        if mast_obj.query_success == "success":
            spoc_avail = True
            print("Found SPOC lightcurve!")
            for flux in self.mast_flux_type:
                rot_obj = myRotations(id_label = 'tic', id_num = row['tic'], lc_df = mast_obj.lc_df, 
                                          flux_type = flux, flux_err_available = True)
                rot_obj.LS_adv()
                rot_obj.get_exo_acf()
                # append LS_results to list
                if rot_obj.LS_results_available == True:
                    temp_df = rot_obj.LS_results
                    if rot_obj.exo_acf_avail == True:
                        temp_df2 = rot_obj.acf_result
                        #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                        temp_df = temp_df.merge(right = temp_df2, on = ['tic','sector'], how = 'left')
                    #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                    #temp_df['type'] = label_type
                    #tic_label = np.repeat(a = row['tic'], repeats = len(temp_df))
                    #temp_df['tic'] = tic_label
                    if flux == 'PDCSAP_flux': 
                        self.pdc_rots_df = pd.concat([self.pdc_rots_df,temp_df])
                        self.pdc_rots_df.to_csv(self.pdc_rot_fn, index = False)
                        pdc_per,pdc_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False: 
                            pdc_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            pdc_ac_per = temp_df['ac_period'].iloc[0]
                    if flux == 'SAP_flux': 
                        self.sap_rots_df = pd.concat([self.sap_rots_df,temp_df])
                        self.sap_rots_df.to_csv(self.sap_rot_fn, index = False)
                        sap_per,sap_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False:
                            sap_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            sap_ac_per = temp_df['ac_period'].iloc[0]
                    print("Rotations added.")
        else:
            ### try to redownload, then try rotations. sometimes reorganization goes wrong/download gets interrupted
            print("Trying to redownload SPOC data.")
            mast_obj = mastObj(tic = row['tic'], ra = row['ra'], dec = row['dec'], download_dir = self.spoc_dir, products = ['LC'], rotations = False)
            mast_obj.download(keep_fits = False)
            if mast_obj.query_success == "success":
                spoc_avail = True
                for flux in self.mast_flux_type:
                    rot_obj = myRotations(id_label = 'tic', id_num = row['tic'], lc_df = mast_obj.lc_df, 
                                          flux_type = flux, flux_err_available = True)
                    rot_obj.LS_adv()
                    rot_obj.get_exo_acf()
                    # append LS_results to list
                    if rot_obj.LS_results_available == True:
                        temp_df = rot_obj.LS_results
                        if rot_obj.exo_acf_avail == True:
                            temp_df2 = rot_obj.acf_result
                            #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                            temp_df = temp_df.merge(right = temp_df2, on = ['tic','sector'], how = 'left')
                        #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                        #temp_df['type'] = label_type
                        #tic_label = np.repeat(a = row['tic'], repeats = len(temp_df))
                        #temp_df['tic'] = tic_label
                        if flux == 'PDCSAP_flux': 
                            self.pdc_rots_df = pd.concat([self.pdc_rots_df,temp_df])
                            self.pdc_rots_df.to_csv(self.pdc_rot_fn, index = False)
                            pdc_per,pdc_power,sect = rot_obj.best_period(rot_obj.LS_results)
                            if np.isnan(np.array(sect, dtype = np.float64)) == False: 
                                pdc_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                            else:
                                pdc_ac_per = temp_df['ac_period'].iloc[0]
                        if flux == 'SAP_flux': 
                            self.sap_rots_df = pd.concat([self.sap_rots_df,temp_df])
                            self.sap_rots_df.to_csv(self.sap_rot_fn, index = False)
                            sap_per,sap_power,sect = rot_obj.best_period(rot_obj.LS_results)
                            if np.isnan(np.array(sect, dtype = np.float64)) == False:
                                sap_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                            else:
                                sap_ac_per = temp_df['ac_period'].iloc[0]
                        print("Rotations added.")
            else:
                print("No SPOC data for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
                sap_per = sap_power = pdc_per = pdc_power = pdc_ac_per = sap_ac_per = np.nan            
                spoc_avail = False
        spoc_df = pd.DataFrame(data={'spoc_avail':[spoc_avail],'pdc_per':[pdc_per],'pdc_power':[pdc_power],
                                     'pdc_ac_per':[pdc_ac_per],
                                     'sap_per':[sap_per],'sap_power':[sap_power],
                                     'sap_ac_per':[sap_ac_per]}, columns = self.spoc_cols)
        return(spoc_df)
    
    def get_eleanor(self,row,i):
        print("Working on Eleanor LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
        eleanor_obj = eleanorObj(tic = row['tic'], ra = row['ra'], dec = row['dec'], download_dir = self.ele_dir, rotations = False)
        eleanor_obj.download(lc_only = True, pca = True)
        if eleanor_obj.no_eleanor == False:
            print("Found Eleanor lightcurve!")
            eleanor_avail = True
            for flux in self.ele_flux_type:
                rot_obj = myRotations(id_label = 'tic',id_num = row['tic'], lc_df = eleanor_obj.LC_df, 
                                      flux_type = flux, flux_err_available = False)
                rot_obj.LS_adv()
                rot_obj.get_exo_acf()
                # append LS_results to list
                if rot_obj.LS_results_available == True:
                    temp_df = rot_obj.LS_results
                    if rot_obj.exo_acf_avail == True:
                        temp_df2 = rot_obj.acf_result
                        #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                        temp_df = temp_df.merge(right = temp_df2, on = ['tic','sector'], how = 'left')
                    #tic_label = np.repeat(a = row['tic'], repeats = len(temp_df))
                    #temp_df['tic'] = tic_label
                    # self.ele_rots_df = pd.concat([self.ele_rots_df,temp_df])
                    # self.ele_rots_df.to_csv(self.ele_rot_fn, index = False)
                    # ele_per,ele_power = rot_obj.best_period(rot_obj.LS_results)
                    
                    if flux == 'corr': 
                        self.ele_corr_rots_df = pd.concat([self.ele_corr_rots_df,temp_df])
                        self.ele_corr_rots_df.to_csv(self.ele_corr_rot_fn, index = False)
                        ele_corr_per,ele_corr_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False:
                            ele_corr_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            ele_corr_ac_per = temp_df['ac_period'].iloc[0]
                    if flux == 'pca': 
                        self.ele_pca_rots_df = pd.concat([self.ele_pca_rots_df,temp_df])
                        self.ele_pca_rots_df.to_csv(self.ele_pca_rot_fn, index = False)
                        ele_pca_per,ele_pca_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False:
                            ele_pca_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            ele_pca_ac_per = temp_df['ac_period'].iloc[0]
                    print("Rotations added.")
        else:
            print("No Eleanor LC for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            ele_corr_per = ele_corr_power = ele_pca_per = ele_pca_power = ele_corr_ac_per = ele_pca_ac_per = np.nan
            eleanor_avail = False
        ele_df = pd.DataFrame(data={'ele_avail':[eleanor_avail],'ele_corr_per':[ele_corr_per],
                                    'ele_corr_power':[ele_corr_power],'ele_corr_ac_per':[ele_corr_ac_per],
                                    'ele_pca_per':[ele_pca_per],'ele_pca_power':[ele_pca_power],
                                    'ele_pca_ac_per':[ele_pca_ac_per]}, columns=self.ele_cols)
        return(ele_df)
        
    def get_cpm(self,row,i):
        print("Downlaoding TESS Cut and extracting CPM LC for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
        try:
            tesscut_obj = tesscutObj(tic = row['tic'], ra = row['ra'], dec = row['dec'], download_dir = self.cpm_dir)
            tesscut_obj.download(save_lc = True, keep_tesscut = False)
            if tesscut_obj.query_success == "success":
                tesscut_avail = True
                rotObj = myRotations(id_label = 'tic',id_num = row['tic'], lc_df = tesscut_obj.lc_df, 
                                     flux_type = 'cpm', flux_err_available = False)
                rotObj.LS_adv()
                rotObj.get_exo_acf()
                if rotObj.LS_results_available == True: 
                    temp_df = rotObj.LS_results
                    if rotObj.exo_acf_avail == True:
                        temp_df2 = rotObj.acf_result
                        #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                        temp_df = temp_df.merge(right = temp_df2, on = ['tic','sector'], how = 'left')
                    
                    #tic_label = np.repeat(a = row['tic'], repeats = len(temp_df))
                    #temp_df['tic'] = tic_label
                    self.cpm_rots_df = pd.concat([self.cpm_rots_df,temp_df])
                    self.cpm_rots_df.to_csv(self.cpm_rot_fn, index = False)
                    cpm_per,cpm_power,sect = rotObj.best_period(rotObj.LS_results)
                    if np.isnan(np.array(sect, dtype = np.float64)) == False:
                        cpm_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                    else:
                        cpm_ac_per = temp_df['ac_period']
                    print("Rotation added.")
            else:
                print("No TESS Cut/CPM LC for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
                cpm_per = cpm_power = cpm_ac_per = np.nan
                tesscut_avail = False
        except:
            print("No TESS Cut/CPM LC for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            cpm_per = cpm_power = cpm_ac_per = np.nan
            tesscut_avail = False
        cpm_df = pd.DataFrame(data={'tesscut_avail':[tesscut_avail],'cpm_per':[cpm_per],
                                    'cpm_power':[cpm_power],'cpm_ac_per':[cpm_ac_per]}, 
                              columns = self.cpm_cols)
        return(cpm_df)
    
    def get_k2(self,row,i,k2sff):
        if k2sff == True:
            print("Working on K2SFF LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2_dir)
        else:
            print("Working on K2 LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2_dir, products = ['LLC'])
        k2_obj.download(keep_fits = False, K2SFF = k2sff)
        if k2_obj.query_success == "success":
            k2_avail = True
            print("Found K2 lightcurve!")
            for flux in self.mast_flux_type:
                rot_obj = myRotations(id_label = 'epic',id_num = row['epic'], lc_df = k2_obj.lc_df, 
                                      flux_type = flux, flux_err_available = True)
                rot_obj.LS_adv()
                rot_obj.get_exo_acf()
                # append LS_results to list
                if rot_obj.LS_results_available == True:
                    temp_df = rot_obj.LS_results
                    if rot_obj.exo_acf_avail == True:
                        temp_df2 = rot_obj.acf_result
                        #temp_df2 = temp_df2[['epic','ac_period']]
                        temp_df = temp_df.merge(right = temp_df2, on = ['epic','sector'], how = 'left')
                    #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                    #temp_df['type'] = label_type
                    if flux == 'PDCSAP_flux': 
                        self.k2_pdc_rots_df = pd.concat([self.k2_pdc_rots_df,temp_df])
                        self.k2_pdc_rots_df.to_csv(self.k2_pdc_rot_fn, index = False)
                        pdc_per,pdc_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False:
                            pdc_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            pdc_ac_per = temp_df['ac_period'].iloc[0]
                    if flux == 'SAP_flux': 
                        self.k2_sap_rots_df = pd.concat([self.k2_sap_rots_df,temp_df])
                        self.k2_sap_rots_df.to_csv(self.k2_sap_rot_fn, index = False)
                        sap_per,sap_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        if np.isnan(np.array(sect, dtype = np.float64)) == False:
                            sap_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                        else:
                            sap_ac_per = temp_df['ac_period'].iloc[0]
                    print("Rotations added.")
        else:
            ### try to redownload, then try rotations. sometimes reorganization goes wrong/download gets interrupted
            print("Trying to redownload K2 data.")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2_dir, products = ['LLC'])
            k2_obj.download(keep_fits = False)
            if k2_obj.query_success == "success":
                k2_avail = True
                for flux in self.mast_flux_type:
                    rot_obj = myRotations(id_label = 'epic',id_num = row['epic'], lc_df = k2_obj.lc_df, 
                                          flux_type = flux, flux_err_available = True)
                    rot_obj.LS_adv()
                    rot_obj.get_exo_acf()
                    # append LS_results to list
                    if rot_obj.LS_results_available == True:
                        temp_df = rot_obj.LS_results
                        if rot_obj.exo_acf_avail == True:
                            temp_df2 = rot_obj.acf_result
                            #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                            temp_df = temp_df.merge(right = temp_df2, on = ['epic','sector'], how = 'left')
                        #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                        #temp_df['type'] = label_type
                        if flux == 'PDCSAP_flux': 
                            self.k2_pdc_rots_df = pd.concat([self.k2_pdc_rots_df,temp_df])
                            self.k2_pdc_rots_df.to_csv(self.k2_pdc_rot_fn, index = False)
                            pdc_per,pdc_power,sect = rot_obj.best_period(rot_obj.LS_results)
                            if np.isnan(np.array(sect, dtype = np.float64)) == False:
                                pdc_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                            else:
                                pdc_ac_per = temp_df['ac_period'].iloc[0]
                        if flux == 'SAP_flux': 
                            self.k2_sap_rots_df = pd.concat([self.k2_sap_rots_df,temp_df])
                            self.k2_sap_rots_df.to_csv(self.k2_sap_rot_fn, index = False)
                            sap_per,sap_power,sect = rot_obj.best_period(rot_obj.LS_results)
                            if np.isnan(np.array(sect, dtype = np.float64)) == False:
                                sap_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                            else:
                                sap_ac_per = temp_df['ac_period'].iloc[0]
                        print("Rotations added.")
            else:
                print("No K2 data for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
                sap_per = sap_power = pdc_per = pdc_power = pdc_ac_per = sap_ac_per = np.nan            
                k2_avail = False
        k2_df = pd.DataFrame(data={'k2_avail':[k2_avail],'k2_pdc_per':[pdc_per],'k2_pdc_power':[pdc_power],
                                   'k2_pdc_ac_per':[pdc_ac_per],
                                   'k2_sap_per':[sap_per],'k2_sap_power':[sap_power],
                                   'k2_sap_ac_per':[sap_ac_per]}, columns = self.k2_cols)
        return(k2_df)
    
    def get_k2sff(self,row,i,k2sff):
        if k2sff == True:
            print("Working on K2SFF LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2sff_dir)
        else:
            print("Working on K2 LCs for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2_dir, products = ['LLC'])
        k2_obj.download(keep_fits = False, K2SFF = k2sff)
        if k2_obj.query_success == "success":
            k2sff_avail = True
            print("Found K2SFF lightcurve!")
            for flux in self.k2sff_flux_type:
                rot_obj = myRotations(id_label = 'epic', id_num = row['epic'], lc_df = k2_obj.lc_df, 
                                      flux_type = flux, flux_err_available = False)
                rot_obj.LS_adv()
                rot_obj.get_exo_acf()
                # append LS_results to list
                if rot_obj.LS_results_available == True:
                    temp_df = rot_obj.LS_results
                    if rot_obj.exo_acf_avail == True:
                        temp_df2 = rot_obj.acf_result
                        #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                        temp_df = temp_df.merge(right = temp_df2, on = ['epic','sector'], how = 'left')
                    #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                    #temp_df['type'] = label_type
                    if flux == 'fcor': 
                        self.k2sff_cor_rots_df = pd.concat([self.k2sff_cor_rots_df,temp_df])
                        self.k2sff_cor_rots_df.to_csv(self.k2sff_cor_rot_fn, index = False)
                        fcor_per,fcor_power,sect = rot_obj.best_period(rot_obj.LS_results)
                        try:
                            test_sect = np.isnan(np.array(sect, dtype = np.float64))
                            if test_sect == False: 
                                fcor_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                            else:
                                fcor_ac_per = temp_df['ac_period'].iloc[0]
                        except:
                            fcor_ac_per = np.nan
                    # if flux == 'fraw': 
                    #     self.k2_sap_rots_df = pd.concat([self.k2_sap_rots_df,temp_df])
                    #     self.k2_sap_rots_df.to_csv(self.k2_sap_rot_fn, index = False)
                    #     sap_per,sap_power = rot_obj.best_period(rot_obj.LS_results)
                    print("Rotations added.")                    
        else:
            ### try to redownload, then try rotations. sometimes reorganization goes wrong/download gets interrupted
            print("Trying to redownload K2SFF data.")
            k2_obj = k2obj(epic = row['epic'], ra = row['ra'], dec = row['dec'], download_dir = self.k2sff_dir)
            k2_obj.download(keep_fits = False, K2SFF = k2sff)
            if k2_obj.query_success == "success":
                k2sff_avail = True
                for flux in self.k2sff_flux_type:
                    rot_obj = myRotations(id_label = 'epic', id_num = row['epic'], lc_df = k2_obj.lc_df, 
                                          flux_type = flux, flux_err_available = False)
                    rot_obj.LS_adv()
                    rot_obj.get_exo_acf()
                    # append LS_results to list
                    if rot_obj.LS_results_available == True:
                        temp_df = rot_obj.LS_results
                        if rot_obj.exo_acf_avail == True:
                            temp_df2 = rot_obj.acf_result
                            #temp_df2 = temp_df2[[self.key_col,'ac_period']]
                            temp_df = temp_df.merge(right = temp_df2, on = ['epic','sector'], how = 'left')
                        #label_type = np.repeat(a = flux.split("_")[0], repeats = len(temp_df))
                        #temp_df['type'] = label_type
                        if flux == 'fcor': 
                            self.k2sff_cor_rots_df = pd.concat([self.k2sff_cor_rots_df,temp_df])
                            self.k2sff_cor_rots_df.to_csv(self.k2sff_cor_rot_fn, index = False)
                            fcor_per,fcor_power,sect = rot_obj.best_period(rot_obj.LS_results)
                            try:
                                test_sect = np.isnan(np.array(sect, dtype = np.float64))
                                if test_sect == False: 
                                    fcor_ac_per = temp_df[temp_df['sector'] == sect]['ac_period'].to_numpy()[0]
                                else:
                                    fcor_ac_per = temp_df['ac_period'].iloc[0]
                            except:
                                fcor_ac_per = np.nan
                        # if flux == 'fraw': 
                        #     self.k2_sap_rots_df = pd.concat([self.k2_sap_rots_df,temp_df])
                        #     self.k2_sap_rots_df.to_csv(self.k2_sap_rot_fn, index = False)
                        #     sap_per,sap_power = rot_obj.best_period(rot_obj.LS_results)
                        print("Rotations added.")
            else:
                print("No K2SFF data for object " + str(i+1) + "/" + str(len(self.query_df)) + ".")
                fcor_per = fcor_power = fcor_ac_per = np.nan            
                k2sff_avail = False
        k2sff_df = pd.DataFrame(data={'k2sff_avail':[k2sff_avail],'k2sff_cor_per':[fcor_per],'k2sff_cor_power':[fcor_power],
                                     'k2sff_ac_per':[fcor_ac_per]}, columns = self.k2sff_cols)
        return(k2sff_df)
        
      
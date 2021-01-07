# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 22:21:39 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

#from openpyxl import load_workbook

import os
import shutil

from astroquery.mast import Observations
#from astroquery.mast import Catalogs

#from astroquery.gaia import Gaia

from astropy import units as u

#from astropy.coordinates import SkyCoord

from astropy.io import fits

from astropy.table import Table


class k2obj:    
    def __init__(self, epic=None, ra=None, dec=None, download_dir=None, products = "all"):
        self.epic = str(epic)
        self.ra = ra
        self.dec = dec
        self.download_dir = download_dir
        self.download_folder = os.path.join(download_dir,'epic' + str(epic))
        
    def download(self, keep_fits = False, K2SFF = False):
        """
        Dowload, reorganize, and simple rename desired data products.
        Updates:
            - get SPOC_df and product list ONLY option (maybe separate query() function)
            - simple rename function
        """
        def reorganize(K2SFF): #bring all files to top level
            print("Reorganizing...")
            ### EXTRACT ALL FILES
            #get down to where sector folders are, get list of folder names
            if K2SFF == True:
                get_down = os.path.join(self.download_folder,'mastDownload','HLSP')
            else:
                get_down = os.path.join(self.download_folder,'mastDownload','K2')
            
            sub_folders = os.listdir(get_down)
            
            for sub in sub_folders:
                
                #get file list of subfolder
                sub_path = os.path.join(get_down,sub)
                
                file_list = os.listdir(sub_path)
                
                #copy files into parent self.folder_name
                for file in file_list:
                    
                    file_path = os.path.join(sub_path,file)
                    
                    new_loc = os.path.join(self.download_folder,file)
                    
                    shutil.copyfile(src = file_path, dst = new_loc)
             
            ### DELETE 'mastDownload' FOLDER 
            delete_dir = os.path.join(self.download_folder,'mastDownload')
            shutil.rmtree(delete_dir)
        
        #def rename(): simple rename for all data products: tic#######_sector####_'type'.'type'
        
        def LC_extract(keep_fits,K2SFF): #extract LC data from all sectors, label by sector
            print("Extracting LC...")
            
            # IDENTIFY LC files
            
            file_list = os.listdir(self.download_folder)            
            LC_list = []            
            for file in file_list:
                file_type = os.path.splitext(file)[-2].split("_")[-1]
                if file_type == 'llc':
                    LC_list.append(file)            
            self.LC_list = LC_list            
            # COMBINE LC'S INTO SINGLE DATAFRAME            
            temp_df_list = []                    
            for lc in LC_list:                
                # append LC data
                lc_path = os.path.join(self.download_folder,lc)                
                temp_df = self.get_LC_data(lc_path,K2SFF = K2SFF)                
                # add sector label column
                sect = lc.split("-")[1].split("_")[0].split("c")[1]
                temp_sector = np.repeat(a=sect,repeats = len(temp_df['time']))
                temp_df['sector'] = temp_sector
                
                temp_df_list.append(temp_df)
                
            lc_df = pd.concat(temp_df_list, ignore_index = True)
            self.lc_df = lc_df
            if K2SFF == True:
                lc_file_name = "epic" + self.epic + "_K2SFF_LC.csv"
            else:
                lc_file_name = "epic" + self.epic + "_K2_LC.csv"            
            lc_file_path = os.path.join(self.download_folder,lc_file_name)            
            self.lc_path = lc_file_path
            
            #writer = pd.ExcelWriter(path = lc_file_path, engine = 'xlsxwriter')
            
            #lc_df.to_excel(writer, sheet_name = 'LC data', index = False)
            if keep_fits == False:
                save_path = os.path.join(self.download_dir,lc_file_name)
                lc_df.to_csv(save_path, index = False)
                shutil.rmtree(self.download_folder)
            else:
                lc_df.to_csv(lc_file_path, index = False)

        ### DO THE DOWNLOADING
        
        query_string = "epic " + self.epic
     
        try:
            obs_table = Observations.query_object(query_string, radius = 0.0005*u.deg)
            if K2SFF == True:
                K2_table = obs_table[obs_table['provenance_name'] == 'K2SFF']
            else:
                K2_table = obs_table[obs_table['provenance_name'] == 'K2']            
            self.K2_df = self.tab2df(K2_table)   
            self.timeseries_ids = self.K2_df[self.K2_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
            if K2SFF == True:
                Observations.download_products(self.timeseries_ids, download_dir = self.download_folder, description = ['FITS'])  
            else:
                Observations.download_products(self.timeseries_ids, download_dir = self.download_folder, productSubGroupDescription = ['LLC']) 
            reorganize(K2SFF = K2SFF)            
            LC_extract(keep_fits,K2SFF = K2SFF)           
                     
            self.query_success = "success"            
        except:
            self.query_success = "fail"
            
    def get_LC_data(self, lc_file, K2SFF, remove_nan = True):
        
        with fits.open(lc_file) as fits_file:
            #get light curve data time and flux arrays
            data = fits_file[1].data  
            if K2SFF == True:
                time = data['T']
                fraw = data['FRAW']
                fcor = data['FCOR']
                arclength = data['ARCLENGTH']
                moving = data['MOVING']
                cadenceno = data['CADENCENO']
                
                cols = ['time','fraw','fcor','arclength','moving','cadenceno']
                dat = Table([time,fraw,fcor,arclength,moving,cadenceno], names = cols)
            else:            
                time = data.TIME
                SAP_flux = data.SAP_FLUX
                SAP_flux_err = data.SAP_FLUX_ERR
                PDCSAP_flux = data.PDCSAP_FLUX
                PDCSAP_flux_err = data.PDCSAP_FLUX_ERR
            
                if remove_nan == True:
                    #get rid of all data points that don't have time values or don't have flux values
                    time_nan = ~np.isnan(time)        
                    SAP_flux_nan = ~np.isnan(SAP_flux)        
                    SAP_flux_err_nan = ~np.isnan(SAP_flux)            
                    PDCSAP_flux_nan = ~np.isnan(PDCSAP_flux)        
                    PDCSAP_flux_err_nan = ~np.isnan(PDCSAP_flux)                
                    no_nan = np.logical_and(np.logical_and(time_nan,SAP_flux_nan),SAP_flux_err_nan)            
                    no_nan = np.logical_and(np.logical_and(no_nan,PDCSAP_flux_nan),PDCSAP_flux_err_nan)
                
                    time = time[no_nan]
                    SAP_flux = SAP_flux[no_nan]
                    SAP_flux_err = SAP_flux_err[no_nan]
                    PDCSAP_flux = PDCSAP_flux[no_nan]
                    PDCSAP_flux_err = PDCSAP_flux_err[no_nan]
                    
                cols = ['time', 'SAP_flux', 'SAP_flux_err', 'PDCSAP_flux', 'PDCSAP_flux_err']
                dat = Table([time,SAP_flux,SAP_flux_err,PDCSAP_flux,PDCSAP_flux_err], names = cols)        
        df = self.tab2df(dat)
        df = pd.DataFrame(df, dtype = 'float')        
        return(df)
    
    def tab2df(self,table):
    
        df = pd.DataFrame()
        
        columns = table.colnames
        
        for col in columns:
            df[col] = table[col]
            
        return(df)
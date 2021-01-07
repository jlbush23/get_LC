# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:08:18 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

#from openpyxl import load_workbook

import os
import shutil

from astroquery.mast import Observations
from astroquery.mast import Catalogs

#from astroquery.gaia import Gaia

from astropy import units as u

#from astropy.coordinates import SkyCoord

from astropy.io import fits

from astropy.table import Table
#from astropy.table import Column


 

class mastObj:
    """
    Used for downloading TESS data products from MAST.
    """
    
    def __init__(self, tic=None, ra=None, dec=None, download_dir=None, products = "all", rotations = False):
        """
        Takes in TIC and/or RA/Dec, download directory, and product list.
        Updates: 
            - make tic,ra,dec flexible for float/str input
            - make sure download dir is proper format
            - make sure products is "all" or a list
            - specify ResolveError and No Data Products error exceptions
        """
        
        self.tic  = tic
        self.ra = ra
        self.dec = dec
        self.parent_folder = download_dir
        self.rotations = rotations
         
        if tic == None:
            
                 
            radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
            
            for rad in radii:
                if self.tic == None:
                    query_string = str(self.ra) + " " + str(self.dec) # make sure to have a space between the strings!
                    obs_table = Catalogs.query_object(query_string, radius = rad*u.deg, catalog = "TIC")
                    obs_df = self.tab2df(obs_table)
                    if len(obs_table['ID']) == 1:
                        self.tic = obs_table['ID'][0]
                        
                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]
                        
                        break
                    
                    if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                        temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                        
                        self.tic = temp_obs_df['ID'].iloc[0]
                        
                        self.bp_rp = (temp_obs_df['gaiabp'] - temp_obs_df['gaiarp']).iloc[0]
                                             
                        break
                    
                    if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                        self.tic = obs_table['ID'][0]
                        
                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]                        
                        
                        break
#                    
#                    if len(obs_table[obs_table['typeSrc'] == "tmgaia2"]) == 1:
#                        self.tic = obs_table['ID'][0]
#                        
#                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]
#                        
#                        break
            
            if self.tic == None:
                self.tic = "tic issue"
                self.bp_rp = 9999
        
        if ra == None:
            query_string = "tic " + self.tic # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = self.tab2df(obs_table)
            self.ra = obs_table['ra'][0]
            self.dec = obs_table['dec'][0]
            
        if products == "all":
            self.products = ['LC','DVM','TP','DVT','DVS','DVR']
        else:
            self.products = products
        
        if self.tic != "tic issue":
            if self.tic != None:
                if download_dir != None:
                    self.folder_name = os.path.join(download_dir,self.tic)
                    
#        ### GAIA QUERY
#        radii = np.linspace(start = 0.00027777, stop = 0.0069444, num = 25)
#        
#        for rad in radii:
#            coord = SkyCoord(ra = float(self.ra), dec = float(self.dec), unit = (u.degree,u.degree))
#            rad = u.Quantity(rad,u.deg)
#            
#            r = Gaia.query_object_async(coordinate = coord, radius = rad)
#            
#            if len(r) == 1:
#                self.gaia_info = self.tab2df(r)
#                break
        
        
    def download(self, keep_fits = False):
        """
        Dowload, reorganize, and simple rename desired data products.
        Updates:
            - get SPOC_df and product list ONLY option (maybe separate query() function)
            - simple rename function
        """
        def reorganize(): #bring all files to top level
            print("Reorganizing...")
            ### EXTRACT ALL FILES
            #get down to where sector folders are, get list of folder names
            get_down = os.path.join(self.folder_name,'mastDownload','TESS')
            
            sub_folders = os.listdir(get_down)
            
            for sub in sub_folders:
                
                #get file list of subfolder
                sub_path = os.path.join(get_down,sub)
                
                file_list = os.listdir(sub_path)
                
                #copy files into parent self.folder_name
                for file in file_list:
                    
                    file_path = os.path.join(sub_path,file)
                    
                    new_loc = os.path.join(self.folder_name,file)
                    
                    shutil.copyfile(src = file_path, dst = new_loc)
             
            ### DELETE 'mastDownload' FOLDER 
            delete_dir = os.path.join(self.folder_name,'mastDownload')
            shutil.rmtree(delete_dir)
        
        #def rename(): simple rename for all data products: tic#######_sector####_'type'.'type'
        
        def LC_extract(keep_fits): #extract LC data from all sectors, label by sector
            print("Extracting LC...")
            
            # IDENTIFY LC files
            
            file_list = os.listdir(self.folder_name)            
            LC_list = []            
            for file in file_list:
                file_type = os.path.splitext(file)[-2].split("_")[-1]
                if file_type == 'lc':
                    LC_list.append(file)            
            self.LC_list = LC_list            
            # COMBINE LC'S INTO SINGLE DATAFRAME            
            temp_df_list = []                    
            for lc in LC_list:                
                # append LC data
                lc_path = os.path.join(self.folder_name,lc)                
                temp_df = self.get_LC_data(lc_path)                
                # add sector label column
                sect = lc.split("-")[1].split("s")[-1][-2:]
                temp_sector = np.repeat(a=sect,repeats = len(temp_df['time']))
                temp_df['sector'] = temp_sector
                
                temp_df_list.append(temp_df)
                
            lc_df = pd.concat(temp_df_list, ignore_index = True)
            self.lc_df = lc_df            
            lc_file_name = "tic" + self.tic + "_mast_LC.csv"            
            lc_file_path = os.path.join(self.folder_name,lc_file_name)            
            self.lc_path = lc_file_path
            
            #writer = pd.ExcelWriter(path = lc_file_path, engine = 'xlsxwriter')
            
            #lc_df.to_excel(writer, sheet_name = 'LC data', index = False)
            if keep_fits == False:
                save_path = os.path.join(self.parent_folder,lc_file_name)
                lc_df.to_csv(save_path, index = False)
                shutil.rmtree(self.folder_name)
            else:
                lc_df.to_csv(lc_file_path, index = False)
#            if keep_fits == False:
#                for file in file_list:
#                    file_path = os.path.join(self.folder_name,file)
#                    os.remove(file_path)
                
            
#            if self.rotations == True:
#                flux_type_list = ['PDCSAP', 'SAP']
#                
#                for flux_type in flux_type_list:
#                    print("Working on " + str(flux_type) + " rotations.")
#                
#                    rotObj = myRotations(tic = self.tic, lc_df = self.lc_df, flux_type = flux_type + '_flux', flux_err_available = True)
#    
#                    rotObj.LS_adv()
#                    
#                    if rotObj.periodogram_LS_available == True:
#                        rotObj.periodogram_LS.to_excel(writer, sheet_name = flux_type + '_periodogram',index = False)
#                    if rotObj.LS_results_available == True:
#                        rotObj.LS_results.to_excel(writer, sheet_name = flux_type + '_LS_results', index = False)
#                
#                    del rotObj
#                    
#            writer.save()
#            writer.close()
            
        ### DO THE DOWNLOADING
        
        query_string = "tic " + self.tic
     
        try:
            obs_table = Observations.query_object(query_string, radius = 0.0005*u.deg)        
            SPOC_table = obs_table[obs_table['provenance_name'] == 'SPOC']            
            SPOC_df = self.tab2df(SPOC_table)            
            self.SPOC_df = SPOC_df            
            self.FFI_ids = SPOC_df[SPOC_df['dataproduct_type'] == 'image']['obsid'].to_numpy(dtype='str')            
            self.timeseries_ids = SPOC_df[SPOC_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
            Observations.download_products(self.timeseries_ids, download_dir = self.folder_name, productSubGroupDescription = self.products)            
            reorganize()            
            LC_extract(keep_fits)           
                     
            self.query_success = "success"            
        except:
            self.query_success = "fail"
            
        
    def check_spoc(self):
        query_string = "tic " + self.tic
     
        try:
            obs_table = Observations.query_object(query_string, radius = 0.0005*u.deg)        
            SPOC_df = obs_table[obs_table['provenance_name'] == 'SPOC'].to_pandas()
            self.timeseries_ids = SPOC_df[SPOC_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
            if len(self.timeseries_ids) > 0:
                self.spoc_avail = True
            else:
                self.spoc_avail = False
        except:
            self.spoc_avail = False
            
    
    ### EVENTUALLY
    
#    def info_summary(self): #create (a list ?) of query summary based on returned observation_table
#        """
#        Takes observation_table and create summary info attributes of the mastObj
#        to be fed into query_reports.
#        """
#    
#    
#    def query_report(self): #create PDF report on query, to be later appended with rotation analysis
#        """
#        Takes in self and creates PDF with 
#        sectors of obervations, number of FFI's and DVM reports, and availability of SPOC LCs.
#        """
#    
#    def check_local(self): #check to see what products have already been downloaded
#         """
#         Takes in self and checks download_dir for already downloaded products.
#         """
    
    def get_LC_data(self, lc_file, remove_nan = True):
        
        with fits.open(lc_file) as fits_file:
            #get light curve data time and flux arrays
            data = fits_file[1].data        
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
        
        
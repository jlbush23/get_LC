# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 13:09:29 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

import os
import shutil

from astroquery.mast import Catalogs
from astroquery.mast import Tesscut

from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.wcs import WCS

from astropy.io import fits

from astropy.table import Table
#from astropy.table import Column

import tess_cpm
#import lightkurve as lk


 

class tesscutObj:
    """
    Used for downloading TESS data products from MAST.
    """
    
    def __init__(self, tic=None, ra=None, dec=None, download_dir=None):#, products = "all"):
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
         
        if tic == None:
            
            query_string = str(self.ra) + " " + str(self.dec) # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = self.tab2df(obs_table)
            if len(obs_table['ID']) == 1:
                self.tic = obs_table['ID'][0]
                self.use_tic = True
            else:
                self.tic = "tic issue"
                self.use_tic = False
        else:
            self.use_tic = True
        
        if ra == None:
            query_string = "tic " + self.tic # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = self.tab2df(obs_table)
            self.ra = obs_table['ra'][0]
            self.dec = obs_table['dec'][0]
            
        download_path = os.path.join(self.parent_folder,str(self.tic))
        
        if os.path.exists(download_path):
            self.download_path = download_path
        else:
            self.download_path = download_path
            os.mkdir(download_path)
    
    def download(self, k=100, n=35, choose_pred_pix = "cosine_similarity", save_lc = False, keep_tesscut = False):   
        
        try:        
            if self.use_tic == True:
                #download_path = os.path.join(self.parent_folder,str(self.tic))
                #os.mkdir(download_path)
                #print(download_path)
                #cpm_lc_df_fn = os.path.join(self.download_path,"tic" + str(self.tic) + "_cpm_LC.xlsx")
                #print(cpm_lc_df_fn)
                manifest = Tesscut.download_cutouts(size=31, sector=None, path=self.download_path, inflate=True, objectname="TIC " + str(self.tic))
                self.manifest = manifest
            else:
                #download_path = os.path.join(self.parent_folder,"ra" + str(self.ra) + "_dec" + str(self.dec))
                #os.mkdir(download_path)
                #cpm_lc_df_fn = os.path.join(download_path,"ra" + str(self.ra) + "_dec" + str(self.dec) + "_cpm_LC.xlsx")
                
                cutout_coord = SkyCoord(float(self.ra), float(self.dec), unit="deg")
                manifest = Tesscut.download_products(coordinates=cutout_coord, size = 31, sector = None, path=self.download_path, inflate=True)
                self.manifest = manifest
            
            if len(manifest) > 0:
                self.query_success = "success"
            if len(manifest) == 0:
                self.query_success = "fail"
                shutil.rmtree(self.download_path)
        except:
            self.query_success = "fail"
            
        
        if self.query_success == "success":
            self.lc_extract(k=k, n=n, choose_pred_pix = choose_pred_pix, save_lc = save_lc, keep_tesscut = keep_tesscut)
            
    
    def lc_extract(self, k=100, n=35, choose_pred_pix = "cosine_similarity", save_lc = False, keep_tesscut = False):
        if self.use_tic == True:
            cpm_lc_df_fn = "tic" + str(self.tic) + "_cpm_LC.csv"
            self.save_fn = cpm_lc_df_fn
        if self.use_tic == False:
            cpm_lc_df_fn = "ra" + str(self.ra) + "_dec" + str(self.dec) + "_cpm_LC.csv"
            self.save_fn = cpm_lc_df_fn
        
        TESS_cuts = os.listdir(self.download_path)
        cpm_lc_df_list = []
        
        for cut in TESS_cuts:
            sector = cut.split("-")[1][-2:]
            
            temp_cut_fn = os.path.join(self.download_path,cut)
            
            with fits.open(temp_cut_fn, mode="readonly") as hdu:
                x_cen = int(round(hdu[1].header["1CRPX4"]))
                y_cen = int(round(hdu[1].header["2CRPX4"]))
                         
            temp_source = tess_cpm.Source(temp_cut_fn, remove_bad=True)            
            temp_source.set_aperture(rowlims=[y_cen-1,y_cen+1], collims=[x_cen-1, y_cen+1])            
            temp_source.add_cpm_model(n=n, predictor_method = choose_pred_pix);        
            #temp_source.add_poly_model();        
            temp_source.set_regs([0.01])#, 0.1])            
            temp_source.holdout_fit_predict(k=k)            
            time = temp_source.time            
            flux = temp_source.get_aperture_lc(data_type="cpm_subtracted_flux")            
            sector = np.repeat(a=sector, repeats = len(time))            
            lc_table = Table([time,flux,sector], names = ['time','cpm','sector'])            
            lc_df = self.tab2df(lc_table)            
            cpm_lc_df_list.append(lc_df) 
        
        del temp_source
        
        cpm_lc_df = pd.concat(cpm_lc_df_list)
        if save_lc == True:
            if keep_tesscut == False:
                save_path = os.path.join(self.parent_folder,cpm_lc_df_fn)
                cpm_lc_df.to_csv(save_path, index = False)
                shutil.rmtree(self.download_path)
            else:
                save_path = os.path.join(self.download_path,cpm_lc_df_fn)
                cpm_lc_df.to_csv(save_path, index = False)
        self.lc_df = cpm_lc_df
#        if keep_tesscut == False:
#            for cut in TESS_cuts:
#                file_path = os.path.join(self.download_path,cut)
#                os.remove(file_path)        
        
    def tab2df(self,table):
    
        df = pd.DataFrame()
        
        columns = table.colnames
        
        for col in columns:
            df[col] = table[col]
            
        return(df)      
        

        
        
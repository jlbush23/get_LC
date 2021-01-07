# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:08:35 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

import os

import eleanor
from lightkurve.lightcurve import LightCurve as LC

from astroquery.mast import Catalogs

from astropy import units as u

from astropy.coordinates import SkyCoord

from astropy.io import fits

from astropy.table import Table
from astropy.table import Column


 

class eleanorObj:
    """
    Used for downloading TESS data products from MAST.
    """
    
    def __init__(self, tic=None, ra=None, dec=None, download_dir=None, use_tic = True, use_coord = False, rotations = False):#, products = "all"):
        """
        Takes in TIC and/or RA/Dec, download directory, and product list.
        Updates: 
            - make tic,ra,dec flexible for float/str input
            - make sure download dir is proper format
            - make sure products is "all" or a list
            - specify ResolveError and No Data Products error exceptions
        """
        
        self.tic  = tic
        self.use_tic = use_tic
        self.use_coord = use_coord
        self.ra = ra
        self.dec = dec
        self.parent_folder = download_dir
        self.rotations = rotations
         
        if tic == None:
            
            query_string = str(self.ra) + " " + str(self.dec) # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = self.tab2df(obs_table)
            if len(obs_table['ID']) == 1:
                self.tic = obs_table['ID'][0]
            else:
                self.tic = "tic issue"
        
        if ra == None:
            query_string = "tic " + self.tic # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = self.tab2df(obs_table)
            self.ra = obs_table['ra'][0]
            self.dec = obs_table['dec'][0]
        
        ### GET SOURCE LIST AND SECTOR LIST
        try: 
            if self.use_tic == True:
                self.source_list = eleanor.multi_sectors(sectors = 'all', tic = int(self.tic), tc = True)            
                sector_list = []    
                for source in self.source_list:
                    sector_list.append(source.sector)                
                self.sector_list = sector_list
            if self.use_coord == True:
                self.source_list = eleanor.multi_sectors(sectors = 'all', coords = (self.ra,self.dec), tc = False)            
                sector_list = []    
                for source in self.source_list:
                    sector_list.append(source.sector)                
                self.sector_list = sector_list         
            
            self.observed_yet = True
        except:
            self.observed_yet = False
        
        
#        if products == "all":
#            self.products = ['LC','DVM','TP','DVT','DVS','DVR']
#        else:
#            self.products = products
#        
#        if self.tic != "tic issue":
#            if self.tic != None:
#                self.folder_name = os.path.join(download_dir,self.tic)
                
    def download(self, save_hlsp = False, keep_tesscut = False, lc_only = True, pca = False):
        
        #initialize storage lists for later concat        
        dat_list = []
        tpf_list = []
        aperture_list = []
        LC_df_list = []
        
        #fill lists with source info from each available sector        
        try:
        
            for source in self.source_list:
                
                sector = source.sector
                #print("Working on sector " + str(sector))                
                dat,tpf,aperture,LC_df = self.get_my_eleanor_data(source = source, sector = sector, pca = pca)#, psf = False)                
                dat_list.append(dat)
                tpf_list.append(tpf)
                aperture_list.append(aperture)
                LC_df_list.append(LC_df)
                
            #store dat_list,tpf_list,aperture_list and concat LC_df and store            
            LC_df_all = pd.concat(LC_df_list, ignore_index = True)
            
            self.dat_list = dat_list
            self.tpf_list = tpf_list
            self.aperture_list = aperture_list
            self.LC_df = LC_df_all
            
            ### WRITE ALL INFO TO EXCEL FILE TITLED tic.xlsx             
            #create storage folder for tic if not already created            
            parent_dir_folder_list = os.listdir(self.parent_folder)
            storage_folder = os.path.join(self.parent_folder,str(self.tic))            
            is_in = np.in1d(self.tic,parent_dir_folder_list)[0]
            
            if lc_only == False: 
                if is_in == False:                                            
                    os.mkdir(path=storage_folder)                       
                temp_file_name = os.path.join(storage_folder,"tic" + str(self.tic) + "_eleanor_LCinfo.xlsx")
                writer = pd.ExcelWriter(path = temp_file_name, engine = 'xlsxwriter')            
                LC_df_all.to_excel(writer, sheet_name = 'LC data', index = False)
                
    #            if self.rotations == True:
    #                print("Working on corrected flux rotations.")
    #            
    #                rotObj = myRotations(tic = self.tic, lc_df = self.LC_df, flux_type = 'corr', flux_err_available = False)
    #    
    #                rotObj.LS_adv()
    #                
    #                if rotObj.periodogram_LS_available == True:
    #                    rotObj.periodogram_LS.to_excel(writer, sheet_name = 'corr' + '_periodogram',index = False)
    #                if rotObj.LS_results_available == True:
    #                    rotObj.LS_results.to_excel(writer, sheet_name = 'corr' + '_LS_results',index = False)
    #                
    #                del rotObj
    
                for i,sector in enumerate(self.sector_list):                
                    tpf_list[i].to_excel(writer, sheet_name = 'tpf_sector' + str(sector), index = False, header = False)                
                    aperture_list[i].to_excel(writer, sheet_name = 'aperture_sector' + str(sector), index = False, header = False)      
                writer.save()
                writer.close()
            
                temp_fn = os.path.join(storage_folder,"tic" + str(self.tic) + "_eleanor_LC.csv")
                LC_df_all.to_csv(temp_fn,index = False)    
                if save_hlsp == True:
                    for data in self.dat_list:
                        data.save(directory = storage_folder)    
            else:
                temp_fn = os.path.join(self.parent_folder,"tic" + str(self.tic) + "_eleanor_LC.csv")
                LC_df_all.to_csv(temp_fn,index = False) 
                
            
            self.no_eleanor = False        
        except:
            self.no_eleanor = True
                  
            
    def get_my_eleanor_data(self,source,sector,psf = True,pca = False):
    
        ### initialize return values        
        LC_df = pd.DataFrame()            
        ### get TargetData object
        #print("pca is " + str(pca))
        dat = eleanor.TargetData(source = source, height = 15, width = 15, bkg_size = 31, do_psf=psf, do_pca=pca, try_load = False, save_postcard = False)
        q = dat.quality == 0
        #print(dat.time[q])
        
        tpf = pd.DataFrame(np.nanmedian(dat.tpf, axis = 0))        
        aperture = pd.DataFrame(dat.aperture)
        
        LC_df['time'] = dat.time[q]
        LC_df['raw'] = dat.raw_flux[q]
        LC_df['corr'] = dat.corr_flux[q]
        if psf == True:
            LC_df['psf'] = dat.psf_flux[q]
        if pca == True:
            LC_df['pca'] = dat.pca_flux[q]
        
        sector4df = np.repeat(a=sector,repeats = len(LC_df['time']))
        LC_df['sector'] = sector4df
            
        return(dat,tpf,aperture,LC_df)
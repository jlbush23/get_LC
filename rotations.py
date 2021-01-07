# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:59:02 2020

@author: jlbus
"""
import pandas as pd
import numpy as np

import os
import shutil

#from astroquery.mast import Observations
#from astroquery.mast import Catalogs

from astropy import units as u

from astropy.coordinates import SkyCoord

#from astropy.io import fits

from astropy.timeseries import LombScargle

from astropy.table import Table
from astropy.table import Column

from scipy.stats import binned_statistic as bin_stat
from scipy.interpolate import interp1d
from scipy.signal import medfilt

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

## import df and flux type

class myRotations(object):
    """
    Used for measuring rotation periods.
    """
    
    def __init__(self, id_label, id_num, lc_df, flux_type, flux_err_available = True):
        """
        Takes in light curve identifiers, time, flux, 
        and flux errors.
        """
        self.id_label = id_label
        self.id_num = id_num
        if self.id_label == 'tic':
            self.name = 'TIC ' + str(id_num)
        if self.id_label == 'epic':
            self.name = 'EPIC ' + str(id_num)
        self.lc_df = lc_df
        self.sector  = lc_df['sector']
        self.time = lc_df['time'].to_numpy(dtype='float')
        self.flux = lc_df[flux_type].to_numpy(dtype='float')
        self.flux_type = flux_type
        self.flux_err_available = flux_err_available
        self.sector_list = np.unique(self.sector.dropna())
        
        if flux_err_available == True:
            self.flux_err = lc_df[flux_type + '_err'].to_numpy(dtype='float')
        else:
            self.flux_err = []

        
        #normalize, median bin flux and smooth
        
        #self.flux_norm
        
        #self.flux_norm_smooth
        
        #self.bin_per
        
        #self.bin_flux
        
        #self.bin_flux_smooth
        
    ### GET ROTATION RESULTS
    
    
    ### flag for peak1 > .96*peak2
    
    ### masking for aliases
    
    ### Gaussian shape flag
        
    ### LS_method to reutrn: period,LS_power,LS_Per/Power - 123    

    def LS_simple(self):
        if self.name[:4] == 'EPIC':
            min_freq = 1/50.0
        else:
            min_freq = 1/30.0
        
        time = self.lc_df['time']
        flux = []
        #for sector in self.sector_list        
    
        if self.flux_err_available == False:                
            freq,power = LombScargle(time,flux).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
        if self.flux_err_available == True:
            flux_err = self.lc_df[self.flux_type + '_err']
            freq,power = LombScargle(time,flux,flux_err).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
        
        threshhold = 0.005 
        peak_locs = []
        pow_peaks = []
        per_peaks = [] 
        
        for i in range(len(power)-2):
            if (power[i] > power[i+1]) and (power[i]>power[i-1]) and (power[i]>threshhold):
                if (power[i] > power[i-2]) and (power[i] > power[i+2]):
                    peak_locs.append(i)
                    pow_peaks.append(power[i])
                    per_peaks.append((1/freq)[i])
                
        #take top 3 periods/powers
        sorted_pow_peaks = np.sort(pow_peaks)[::-1]
        top3 = sorted_pow_peaks[0:3]
        per_top3 = []
        for pow in top3:
            per_top3.append((per_peaks)[np.where(pow_peaks == pow)[0][0]])
            
        period = 1/freq
                    
        if len(top3) == 0:
            Per1 = np.nan
            Power1 = np.nan
            Per2 = np.nan
            Power2 = np.nan
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 1:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = np.nan
            Power2 = np.nan
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 2:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 3:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = per_top3[2]
            Power3 = top3[2]
            
        self.simple_LS_periodogram = (pd.DataFrame(data = np.transpose(np.array([period,power])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
        temp_LS_results_tab = Table([[self.id_num],[Per1],[Per2],[Per3],[Power1],[Power2],[Power3]], 
                                    names = [self.id_label,'LS_Per1','LS_Per2','LS_Per3','LS_Power1','LS_Power2','LS_Power3'])
        self.simple_LS_results = self.tab2df(temp_LS_results_tab)
    
    def LS_adv(self):
        if self.name[:4] == 'EPIC':
            min_freq = 1/50.0
        else:
            min_freq = 1/30.0
               
        periodogram_list = []
        LS_results_list = []
        
        for sector in self.sector_list:
            
            time = self.lc_df[self.lc_df['sector'] == sector]['time']
            flux = self.lc_df[self.lc_df['sector'] == sector][self.flux_type]
        
            if self.flux_err_available == False:                
                freq,power = LombScargle(time,flux).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
            if self.flux_err_available == True:
                flux_err = self.lc_df[self.lc_df['sector'] == sector][self.flux_type + '_err']
                freq,power = LombScargle(time,flux,flux_err).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
            
            threshhold = 0.005 
            peak_locs = []
            pow_peaks = []
            per_peaks = [] 
            
            for i in range(len(power)-2):
                if (power[i] > power[i+1]) and (power[i]>power[i-1]) and (power[i]>threshhold):
                    if (power[i] > power[i-2]) and (power[i] > power[i+2]):
                        peak_locs.append(i)
                        pow_peaks.append(power[i])
                        per_peaks.append((1/freq)[i])
                    
            #take top 3 periods/powers
            sorted_pow_peaks = np.sort(pow_peaks)[::-1]
            top3 = sorted_pow_peaks[0:3]
            per_top3 = []
            for pow in top3:
                per_top3.append((per_peaks)[np.where(pow_peaks == pow)[0][0]])
                
            period = 1/freq
            
            temp_sector = np.repeat(a=str(sector),repeats = len(period))
                        
            if len(top3) == 0:
                Per1 = np.nan
                Power1 = np.nan
                Per2 = np.nan
                Power2 = np.nan
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 1:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = np.nan
                Power2 = np.nan
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 2:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = per_top3[1]
                Power2 = top3[1]
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 3:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = per_top3[1]
                Power2 = top3[1]
                Per3 = per_top3[2]
                Power3 = top3[2]
                
            temp_periodogram_df = (pd.DataFrame(data = np.transpose(np.array([period,power])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
            temp_periodogram_df['sector'] = temp_sector
            temp_LS_results_tab = Table([[self.id_num],[sector],[Per1],[Per2],[Per3],[Power1],[Power2],[Power3]], 
                                        names = [self.id_label, 'sector','LS_Per1','LS_Per2','LS_Per3','LS_Power1','LS_Power2','LS_Power3'])
            temp_LS_results_df = self.tab2df(temp_LS_results_tab)
            temp_LS_results_df['sector'] = temp_LS_results_df['sector'].to_numpy(dtype = 'str')
            
            periodogram_list.append(temp_periodogram_df)
            LS_results_list.append(temp_LS_results_df)
        
        try:
            self.periodogram_LS = pd.concat(objs = periodogram_list, ignore_index = True)
            self.LS_results = pd.concat(objs = LS_results_list, ignore_index = True)
            self.periodogram_LS_available = True
            self.LS_results_available = True
        except: 
            self.periodogram_LS_available = False
            self.LS_results_available = False
            
    def get_exo_acf(self):        
        flux_type = self.flux_type
        flux_err_avail = self.flux_err_available
        import exoplanet as xo
        
        #folder,file = os.path.split(lc_fn)
        #read lc file name
        #tic = file.split("_")[0].split("tic")[1]
        #lc_df = pd.read_csv(lc_fn)
        
        #sector_list = np.unique(self.lc_df['sector'].to_numpy())
        
        result_list = []
        self.acf_periodograms = {}
        
        for sector in self.sector_list:
            temp_lc_df = self.lc_df[self.lc_df['sector'] == sector]
            
            time=temp_lc_df['time'].to_numpy()
            flux = temp_lc_df[self.flux_type].to_numpy()
            if self.flux_err_available: flux_err = temp_lc_df[self.flux_type + '_err'].to_numpy()
            
            mu = np.mean(flux)
            flux = (flux / mu - 1)
            if self.flux_err_available: flux_err = flux_err / mu
            
            try:
                #AC Results
                if self.flux_err_available:
                    ac_results = xo.autocorr_estimator(x=time, y=flux, yerr=flux_err, min_period=0.1, max_period=45)
                else:
                    ac_results = xo.autocorr_estimator(x=time, y=flux, min_period=0.1, max_period=45)
                
                ac_peak1 = ac_results['peaks'][0]
                ac_period = ac_peak1['period']
                
                temp_periodogram = pd.DataFrame(data = np.transpose(ac_results['autocorr']), columns = ['Period','AC Power'])
                #ac_error = ac_peak1['period_uncert']
                self.acf_periodograms[sector] = temp_periodogram
            except:
                ac_period = np.nan
            
            # temp_result = pd.DataFrame(data = {'sector':str(sector),'ls_period1':[ls_period],'ls_logpower':[ls_power],
            #                               'ls_error':[ls_error],'ac_period':[ac_period]}, 
            #                       columns = ['sector','ls_period1','ls_logpower','ls_error','ac_period'])
            
            temp_result = pd.DataFrame(data = {'sector':[str(sector)],'ac_period':[ac_period]}, 
                                  columns = ['sector','ac_period'])            
            result_list.append(temp_result)
            
            
            
        self.acf_result = pd.concat(result_list)
        id_repeats = np.repeat(self.id_num, len(self.acf_result))
        self.acf_result.insert(loc = 0, column = self.id_label, value = id_repeats)
        self.exo_acf_avail = True
            
    def AC_simple(self):
        time = self.lc_df['time']
        flux = self.lc_df[self.flux_type]
  
        #interpolate and get new flux
        num_points = int(round((np.max(time)-np.min(time))*24*60))
        #num_points = int(round(len(time)*10))
        
        time_new = np.linspace(start=np.min(time),stop=np.max(time),num=num_points)
        
        npts = len(time_new)
        lags = np.arange(-npts + 1, npts)
        lags_pos = lags[lags > .05]
        lags_neg = lags[lags < -0.05]
        
        f = interp1d(time, flux, kind = 'linear')
        
        f_interp = f(time_new)
        
        autocov = np.correlate(f_interp - np.mean(f_interp),f_interp - np.mean(f_interp), mode = 'full')
        autocorrelation = autocov/(len(f_interp)*(f_interp.std()**2))
        
        auto_pos = autocorrelation[lags > 0.05]
        
        #ac_roll = medfilt(auto_pos, kernel_size = 20)
        
        bin_auto,bin_lag,bin_AC_idx = bin_stat(x=lags_pos,values=auto_pos,statistic = 'median', bins = 0.05*len(auto_pos))
        
        bin_lag = bin_lag[0:len(bin_lag)-1]
        
        width = np.median(np.gradient(time_new)) ## this is the offset between x points
        
        bin_per = bin_lag*width
        
        #period = lags_pos*width
        
#            def fwhm2sigma(fwhm):
#                return fwhm / np.sqrt(8 * np.log(2))
#        
#        
#            FWHM = .5
#            sigma = fwhm2sigma(FWHM)
#            
#            ac_smooth = np.zeros(auto_pos.shape)
#            for i,x_position in enumerate(period):
#                kernel = np.exp(-(period - x_position) ** 2 / (2 * sigma ** 2))
#                kernel = kernel/sum(kernel)
#                ac_smooth[i] = sum(auto_pos * kernel)
  

        #find first minimum
        threshhold = 0.1
        min_locs = []
        pow_min = []
        per_min = []
        
        for i in range(len(bin_auto)-2):
            if (bin_auto[i] < bin_auto[i+1]) and (bin_auto[i]<bin_auto[i-1]) and (bin_auto[i]<threshhold):
                #if (bin_auto[i] < bin_auto[i-2]) and (bin_auto[i] < bin_auto[i+2]):
                min_locs.append(i)
                pow_min.append(bin_auto[i])
                per_min.append(bin_per[i])
        
#            for i in range(len(ac_smooth)-2):
#                if (ac_smooth[i] < ac_smooth[i+1]) and (ac_smooth[i]<ac_smooth[i-1]) and (ac_smooth[i]<threshhold):
#                    #if (bin_auto[i] < bin_auto[i-2]) and (bin_auto[i] < bin_auto[i+2]):
#                    min_locs.append(i)
#                    pow_min.append(ac_smooth[i])
#                    per_min.append(period[i])
#                    
        #take top 3 periods/powers
        sorted_pow_min = np.sort(pow_min)
        if len(sorted_pow_min) > 0:
            bottom = sorted_pow_min[0]
            per_bottom = (per_min)[np.where(pow_min == bottom)[0][0]]
        if len(sorted_pow_min) == 0:
            bottom = np.nan
            per_bottom = np.nan
            
        #find peak after first minimum
        threshhold = 0.018 
        peak_locs = []
        pow_peaks = []
        per_peaks = []
        
        if np.isnan(per_bottom) == False:
            check_auto_here = bin_auto[bin_per>per_bottom]
            check_bin_per = bin_per[bin_per>per_bottom]
        if np.isnan(per_bottom) == True:
            check_auto_here = bin_auto[bin_per>0.19]
            check_bin_per = bin_per[bin_per>0.19]
        
#            if np.isnan(per_bottom) == False:
#                check_auto_here = ac_smooth[period>per_bottom]
#                check_bin_per = period[period>per_bottom]
#            if np.isnan(per_bottom) == True:
#                check_auto_here = ac_smooth[period>0.19]
#                check_bin_per = period[period>0.19]
        
        for i in range(len(check_auto_here)-1):
            if (check_auto_here[i] > check_auto_here[i+1]) and (check_auto_here[i]>check_auto_here[i-1]) and (check_auto_here[i]>threshhold):
                #if (check_auto_here[i] > check_auto_here[i-2]) and (check_auto_here[i] > check_auto_here[i+2]):
                peak_locs.append(i)
                pow_peaks.append(check_auto_here[i])
                per_peaks.append(check_bin_per[i])
                
        #take top 3 periods/powers
        sorted_pow_peaks = np.sort(pow_peaks)[::-1]
        top3 = sorted_pow_peaks[0:3]
        per_top3 = []
        for pow in top3:
            per_top3.append((per_peaks)[np.where(pow_peaks == pow)[0][0]])
        
        if len(top3) == 0:
            Per1 = np.nan
            Power1 = np.nan
            Per2 = np.nan
            Power2 = np.nan
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 1:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = np.nan
            Power2 = np.nan
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 2:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = np.nan
            Power3 = np.nan
        
        if len(top3) == 3:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = per_top3[2]
            Power3 = top3[2]
            
        self.simple_AC_periodogram = (pd.DataFrame(data = np.transpose(np.array([bin_per,bin_auto])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
        temp_AC_results_tab = Table([[self.id_num],[Per1],[Per2],[Per3],[Power1],[Power2],[Power3]], 
                                    names = [self.id_label,'AC_Per1','AC_Per2','AC_Per3','AC_Power1','AC_Power2','AC_Power3'])
        self.simple_AC_results = self.tab2df(temp_AC_results_tab)
            
    def AC_method(self):
        
        periodogram_list = []
        
        AC_results_list = []
        
        for sector in self.sector_list:
            
            time = self.lc_df[self.lc_df['sector'] == sector]['time']
            flux = self.lc_df[self.lc_df['sector'] == sector][self.flux_type]
      
            #interpolate and get new flux
            num_points = int(round((np.max(time)-np.min(time))*24*60))
            #num_points = int(round(len(time)*10))
            
            time_new = np.linspace(start=np.min(time),stop=np.max(time),num=num_points)
            
            npts = len(time_new)
            lags = np.arange(-npts + 1, npts)
            lags_pos = lags[lags > .05]
            lags_neg = lags[lags < -0.05]
            
            f = interp1d(time, flux, kind = 'linear')
            
            f_interp = f(time_new)
            
            autocov = np.correlate(f_interp - np.mean(f_interp),f_interp - np.mean(f_interp), mode = 'full')
            autocorrelation = autocov/(len(f_interp)*(f_interp.std()**2))
            
            auto_pos = autocorrelation[lags > 0.05]
            
            #ac_roll = medfilt(auto_pos, kernel_size = 20)
            
            bin_auto,bin_lag,bin_AC_idx = bin_stat(x=lags_pos,values=auto_pos,statistic = 'median', bins = 0.05*len(auto_pos))
            
            bin_lag = bin_lag[0:len(bin_lag)-1]
            
            width = np.median(np.gradient(time_new)) ## this is the offset between x points
            
            bin_per = bin_lag*width
            
            #period = lags_pos*width
            
#            def fwhm2sigma(fwhm):
#                return fwhm / np.sqrt(8 * np.log(2))
#        
#        
#            FWHM = .5
#            sigma = fwhm2sigma(FWHM)
#            
#            ac_smooth = np.zeros(auto_pos.shape)
#            for i,x_position in enumerate(period):
#                kernel = np.exp(-(period - x_position) ** 2 / (2 * sigma ** 2))
#                kernel = kernel/sum(kernel)
#                ac_smooth[i] = sum(auto_pos * kernel)
      
            temp_sector = np.repeat(a=sector,repeats = len(bin_per))
            
            #find first minimum
            threshhold = 0.1
            min_locs = []
            pow_min = []
            per_min = []
            
            for i in range(len(bin_auto)-2):
                if (bin_auto[i] < bin_auto[i+1]) and (bin_auto[i]<bin_auto[i-1]) and (bin_auto[i]<threshhold):
                    #if (bin_auto[i] < bin_auto[i-2]) and (bin_auto[i] < bin_auto[i+2]):
                    min_locs.append(i)
                    pow_min.append(bin_auto[i])
                    per_min.append(bin_per[i])
            
#            for i in range(len(ac_smooth)-2):
#                if (ac_smooth[i] < ac_smooth[i+1]) and (ac_smooth[i]<ac_smooth[i-1]) and (ac_smooth[i]<threshhold):
#                    #if (bin_auto[i] < bin_auto[i-2]) and (bin_auto[i] < bin_auto[i+2]):
#                    min_locs.append(i)
#                    pow_min.append(ac_smooth[i])
#                    per_min.append(period[i])
#                    
            #take top 3 periods/powers
            sorted_pow_min = np.sort(pow_min)
            if len(sorted_pow_min) > 0:
                bottom = sorted_pow_min[0]
                per_bottom = (per_min)[np.where(pow_min == bottom)[0][0]]
            if len(sorted_pow_min) == 0:
                bottom = np.nan
                per_bottom = np.nan
                
            #find peak after first minimum
            threshhold = 0.018 
            peak_locs = []
            pow_peaks = []
            per_peaks = []
            
            if np.isnan(per_bottom) == False:
                check_auto_here = bin_auto[bin_per>per_bottom]
                check_bin_per = bin_per[bin_per>per_bottom]
            if np.isnan(per_bottom) == True:
                check_auto_here = bin_auto[bin_per>0.19]
                check_bin_per = bin_per[bin_per>0.19]
            
#            if np.isnan(per_bottom) == False:
#                check_auto_here = ac_smooth[period>per_bottom]
#                check_bin_per = period[period>per_bottom]
#            if np.isnan(per_bottom) == True:
#                check_auto_here = ac_smooth[period>0.19]
#                check_bin_per = period[period>0.19]
            
            for i in range(len(check_auto_here)-1):
                if (check_auto_here[i] > check_auto_here[i+1]) and (check_auto_here[i]>check_auto_here[i-1]) and (check_auto_here[i]>threshhold):
                    #if (check_auto_here[i] > check_auto_here[i-2]) and (check_auto_here[i] > check_auto_here[i+2]):
                    peak_locs.append(i)
                    pow_peaks.append(check_auto_here[i])
                    per_peaks.append(check_bin_per[i])
                    
            #take top 3 periods/powers
            sorted_pow_peaks = np.sort(pow_peaks)[::-1]
            top3 = sorted_pow_peaks[0:3]
            per_top3 = []
            for pow in top3:
                per_top3.append((per_peaks)[np.where(pow_peaks == pow)[0][0]])
            
            if len(top3) == 0:
                Per1 = np.nan
                Power1 = np.nan
                Per2 = np.nan
                Power2 = np.nan
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 1:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = np.nan
                Power2 = np.nan
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 2:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = per_top3[1]
                Power2 = top3[1]
                Per3 = np.nan
                Power3 = np.nan
            
            if len(top3) == 3:
                Per1 = per_top3[0]
                Power1 = top3[0]
                Per2 = per_top3[1]
                Power2 = top3[1]
                Per3 = per_top3[2]
                Power3 = top3[2]
                
            temp_periodogram_df = (pd.DataFrame(data = np.transpose(np.array([bin_per,bin_auto])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
            temp_periodogram_df['sector'] = temp_sector
            temp_AC_results_tab = Table([[self.id_num],[sector],[Per1],[Per2],[Per3],[Power1],[Power2],[Power3]], 
                                        names = [self.id_label, 'sector','AC_Per1','AC_Per2','AC_Per3','AC_Power1','AC_Power2','AC_Power3'])
            temp_AC_results_df = self.tab2df(temp_AC_results_tab)
            
            periodogram_list.append(temp_periodogram_df)
            AC_results_list.append(temp_AC_results_df)
        
        try:
            self.periodogram_AC = pd.concat(objs = periodogram_list, ignore_index = True)
            self.AC_results = pd.concat(objs = AC_results_list, ignore_index = True)
            self.periodogram_AC_available = True
            self.AC_results_available = True
        except: 
            self.periodogram_AC_available = False
            self.AC_results_available = False
            
        #return(bin_per,bin_auto,Per1,Per2,Per3,Power1,Power2,Power3)
        
    def best_period(self,rot_df):
        temp_Power1 = (rot_df['LS_Power1'].dropna()).to_numpy(dtype='float')        
        if len(temp_Power1) > 0:        
            row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
            best_df = rot_df.iloc[row_idx,:]
            best_per = best_df['LS_Per1']
            best_pow = best_df['LS_Power1']
            best_sect = str(best_df['sector'])
            return((best_per,best_pow,best_sect))
        else:
            return((np.nan,np.nan,np.nan))
            
        
    def simple_period_graph(self, MAST = False, eleanor = False, flux_type = None):
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        #This stuff is everything, use it for any python plot to make it nicer.
        mpl.rcParams['lines.linewidth'] =3
        mpl.rcParams['axes.linewidth'] = 2
        mpl.rcParams['xtick.major.width'] =2
        mpl.rcParams['ytick.major.width'] =2
        mpl.rcParams['xtick.minor.width'] =1.5
        mpl.rcParams['ytick.minor.width'] =1.5
        mpl.rcParams['ytick.labelsize'] = 17
        mpl.rcParams['xtick.labelsize'] = 17
        mpl.rcParams['axes.labelsize'] = 17
        #mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['axes.labelweight']='semibold'
        mpl.rcParams['mathtext.fontset']='stix'
        mpl.rcParams['font.weight'] = 'semibold'
        mpl.rcParams['axes.titleweight']='semibold'
        mpl.rcParams['axes.titlesize']=17
        
        AC_power = self.simple_AC_periodogram['power']
        
        AC_height = np.max(AC_power) - np.min(AC_power)
        
        #make plots
    
        fig = plt.figure(figsize=(15,12))
        gridspec.GridSpec(3,2)
        
        ### PLOT FULL LC
        plt.subplot2grid((3,2), (0,0), colspan = 2)
        
        for i,sector in enumerate(self.sector_list):
            
            time = self.lc_df[self.lc_df['sector'] == sector]['time']
            flux = self.lc_df[self.lc_df['sector'] == sector][self.flux_type]
            
            bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.05*len(time)))
            bin_time = bin_time[0:len(bin_time)-1]
                  
            if self.flux_type != 'cpm':
                bin_flux = bin_flux/np.nanmedian(bin_flux)
                flux = flux/np.median(flux) 
                plt.scatter(time,flux, s = 0.75, label = sector)
                plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
            if self.flux_type == 'cpm':
                plt.scatter(time,flux, s = 1, label = sector)#, c = c_sector[i])
                #plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)
        
        if self.flux_type == 'cpm':
            plt.ylim((-0.15,0.15))
        plt.xlabel("Time (days)")
        plt.ylabel("Normalized Flux")
        plt.legend(loc = 'upper right', fontsize = 'xx-large')
        plt.title("Light Curve of " + self.name)
        
        ### PLOT LS periodograms
    
        plt.subplot2grid((3,2), (1,0),colspan=1)
        
        for i,sector in enumerate(self.sector_list):
            LS_results = self.LS_results[self.LS_results['sector'] == sector]
            
            period = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['period']
            power = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['power']
            
#            ymax1 = LS_results['LS_Power1'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax2 = LS_results['LS_Power2'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax3 = LS_results['LS_Power3'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            
            plt.plot(period,power)#,c = c_sector[i])
            plt.scatter(LS_results['LS_Per1'],LS_results['LS_Power1'], c = 'r')
            plt.scatter(LS_results['LS_Per2'],LS_results['LS_Power2'], c = 'r')
            plt.scatter(LS_results['LS_Per3'],LS_results['LS_Power3'], c = 'r')
#            plt.axvline(x=LS_results['LS_Per1'].to_numpy()[0], ymin=0.01,ymax = ymax1, c= 'r', linewidth = 4)
#            plt.axvline(x=LS_results['LS_Per2'].to_numpy()[0], ymin=0.01,ymax = ymax2, c= 'r', linewidth = 2)
#            plt.axvline(x=LS_results['LS_Per3'].to_numpy()[0], ymin=0.01,ymax = ymax3, c= 'r', linewidth = 1)
          
        plt.plot(self.simple_LS_periodogram['period'],self.simple_LS_periodogram['power'], c = 'black', 
                 label = str(self.simple_LS_results['LS_Per1'][0]) + ' days')
        plt.scatter(self.simple_LS_results['LS_Per1'][0],self.simple_LS_results['LS_Power1'][0], c = 'r')
        plt.scatter(self.simple_LS_results['LS_Per2'][0],self.simple_LS_results['LS_Power2'][0], c = 'r')
        plt.scatter(self.simple_LS_results['LS_Per3'][0],self.simple_LS_results['LS_Power3'][0], c = 'r')
       
        temp_Power1 = (self.LS_results['LS_Power1'].dropna()).to_numpy(dtype='float')
    
        if len(temp_Power1) > 0:        
            row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
            best_per = self.LS_results['LS_Per1'].iloc[row_idx]            
            best_sector = self.LS_results['sector'].iloc[row_idx]
        else:
            best_period = best_sector = row_idx = np.nan
        
        plt.xlim((0,30))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("LS Power")
        plt.title("LS Period = " + str(round(best_per,4)) + " - Found in Sector" + str(best_sector))
        plt.legend(loc = 'upper left')
    
        ### PLOT AC PERIODOGRAMS
        
        plt.subplot2grid((3,2), (1,1),colspan=1)
        
        for i,sector in enumerate(self.sector_list):
            AC_results = self.acf_result[self.acf_result['sector'] == sector]
            period = self.acf_periodograms[sector]['Period']
            power = self.acf_periodograms[sector]['AC Power']
            
            plt.plot(period,power)#, c = c_sector[i])
            #plt.scatter(AC_results['ac_period'],AC_results['AC_Power1'], c = 'r')
            # plt.scatter(AC_results['AC_Per2'],AC_results['AC_Power2'], c = 'r')
            # plt.scatter(AC_results['AC_Per3'],AC_results['AC_Power3'], c = 'r')
            plt.axvline(x=AC_results['ac_period'],color = 'red', ymin=0.01, ymax = 0.99,linewidth = 4)
#            plt.axvline(x=AC_results['AC_Per2'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power2'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 2)
#            plt.axvline(x=AC_results['AC_Per3'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power3'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 1)
#       
        # plt.plot(self.simple_AC_periodogram['period'],self.simple_AC_periodogram['power'], label = str(self.simple_AC_results['AC_Per1'][0]) + ' days')
        # plt.scatter(self.simple_AC_results['AC_Per1'][0],self.simple_AC_results['AC_Power1'][0], c = 'r')
        # plt.scatter(self.simple_AC_results['AC_Per2'][0],self.simple_AC_results['AC_Power2'][0], c = 'r')
        # plt.scatter(self.simple_AC_results['AC_Per3'][0],self.simple_AC_results['AC_Power3'][0], c = 'r')
        
        ## title for ac periodogram
        # if np.isnan(row_idx) == False:
        #     assoc_AC_per = self.AC_results['AC_Per1'].iloc[row_idx]
        # else:
        #     assoc_AC_per = np.nan
        
        ac_title = ''
        for i,sector in enumerate(self.sector_list):
            AC_results = self.acf_result[self.acf_result['sector'] == sector]
            if i>0: ac_title = ac_title + ', '
            ac_title = ac_title + 'Sector ' + str(sector) + ' : ' + str(AC_results['ac_period']) + ' d'
            
        
        plt.xlim((0,30))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("ACF Power")
        plt.title(ac_title)
        
        ### PLOT PHASE PHOLDED CURVES FOR BEST SECTOR
        
        time_best_sector = self.lc_df['time'].to_numpy(dtype = 'float')
        flux_best_sector = self.lc_df[self.flux_type].to_numpy(dtype = 'float')
        
        bin_flux,bin_time,bin_idx = bin_stat(x=time_best_sector,values=flux_best_sector,statistic = 'median', bins = round(0.05*len(time)))
        bin_time = bin_time[0:len(bin_time)-1]
        bin_flux = bin_flux/np.nanmedian(bin_flux)
        flux = flux/np.median(flux)    
        
        #normalized phases arrays
        
        # simple_LS_per = self.simple_LS_results['LS_Per1'][0]
        # simple_AC_per = self.simple_AC_results['AC_Per1'][0]
        
        if self.flux_type != 'cpm':
            phase_norm_LS = (bin_time % simple_LS_per)/simple_LS_per #from LS method        
            phase_norm_AC = (bin_time % assoc_AC_per)/assoc_AC_per #from AC
        if self.flux_type == 'cpm':
            phase_norm_LS = (time_best_sector % simple_LS_per)/simple_LS_per #from LS method 
            phase_norm_AC = (time_best_sector % assoc_AC_per)/assoc_AC_per #from AC
        #detect changes in phase
    
        c1 = []
        for i in range(len(phase_norm_LS)-1):
            if phase_norm_LS[i+1] - phase_norm_LS[i] < 0:
                c1.append(i+1)
    
        c2 = []
        for i in range(len(phase_norm_AC)-1):
            if phase_norm_AC[i+1] - phase_norm_AC[i] < 0:
                c2.append(i+1)
        
        # LS phase fold
    
        plt.subplot2grid((3,2), (2,0),colspan=1)
        first=0
        if self.flux_type != 'cpm':
            for change in c1:
                plt.scatter(phase_norm_LS[first:change],bin_flux[first:change])
                first=change+1
        if self.flux_type == 'cpm':
            for change in c1:
                plt.scatter(phase_norm_LS[first:change],flux_best_sector[first:change])
                first=change+1
        plt.xlabel("Fraction of Period")
        plt.ylabel("Normalized Flux")
        plt.title("LS Phase-Folded - Full LC")
        
        # AC phase fold
    
        plt.subplot2grid((3,2), (2,1),colspan=1)
        first = 0
        if self.flux_type != 'cpm':
            for change in c2:
                plt.scatter(phase_norm_AC[first:change],bin_flux[first:change])
                first=change+1
        if self.flux_type == 'cpm':
            for change in c2:
                plt.scatter(phase_norm_AC[first:change],flux_best_sector[first:change])
                first=change+1
                
        plt.xlabel("Fraction of Period")
        plt.ylabel("Normalized Flux")
        plt.title("AC Phase-Folded")#" - Full LC")
        
        fig.tight_layout()
        plt.close(fig=fig)
        
        #decide on best period based on highest power.... not sure if this is great still. In fact, I know its not the end of the story...
    #    best_power_LS = ((power)[np.where(power == np.max(power))])[0]
    #    best_power_AC = ((AC_power)[np.where(AC_per == bestper_AC)])[0]
    #    
        #period_results = np.array([bestper_LS,best_power_LS,bestper_AC,best_power_AC])
        
    #     choose = np.array([best_power1,best_power2])
    #     between = np.array([bestper1,bestper2])
    #     #choice = ((between)[np.where(choose == np.argmax(choose))])[0]
        
        self.simple_rotFig = fig
    
    def period_graph(self, MAST = False, eleanor = False, flux_type = None):
       
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        #This stuff is everything, use it for any python plot to make it nicer.
        mpl.rcParams['lines.linewidth'] =3
        mpl.rcParams['axes.linewidth'] = 2
        mpl.rcParams['xtick.major.width'] =2
        mpl.rcParams['ytick.major.width'] =2
        mpl.rcParams['xtick.minor.width'] =1.5
        mpl.rcParams['ytick.minor.width'] =1.5
        mpl.rcParams['ytick.labelsize'] = 17
        mpl.rcParams['xtick.labelsize'] = 17
        mpl.rcParams['axes.labelsize'] = 17
        #mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['axes.labelweight']='semibold'
        mpl.rcParams['mathtext.fontset']='stix'
        mpl.rcParams['font.weight'] = 'semibold'
        mpl.rcParams['axes.titleweight']='semibold'
        mpl.rcParams['axes.titlesize']=17
        
        # AC_power = self.periodogram_AC['power']
        
        # AC_height = np.max(AC_power) - np.min(AC_power)
        
        #make plots
    
        fig = plt.figure(figsize=(15,12))
        gridspec.GridSpec(3,2)
        
        ### PLOT FULL LC
        plt.subplot2grid((3,2), (0,0), colspan = 2)
        
        for i,sector in enumerate(self.sector_list):
            
            time = self.lc_df[self.lc_df['sector'] == sector]['time']
            flux = self.lc_df[self.lc_df['sector'] == sector][self.flux_type]
            
            bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.05*len(time)))
            bin_time = bin_time[0:len(bin_time)-1]
                  
            if self.flux_type != 'cpm':
                bin_flux = bin_flux/np.nanmedian(bin_flux)
                flux = flux/np.median(flux) 
                plt.scatter(time,flux, s = 0.75, label = sector)
                plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
            if self.flux_type == 'cpm':
                plt.scatter(time,flux, s = 1, label = sector)#, c = c_sector[i])
                #plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)
        
        if self.flux_type == 'cpm':
            plt.ylim((-0.15,0.15))
        plt.xlabel("Time (days)")
        plt.ylabel("Normalized Flux")
        plt.legend(loc = 'upper right', fontsize = 'xx-large')
        plt.title("Light Curve of " + self.name)
        
        ### PLOT LS periodograms
    
        plt.subplot2grid((3,2), (1,0),colspan=1)
        
        for i,sector in enumerate(self.sector_list):
            LS_results = self.LS_results[self.LS_results['sector'] == sector]
            
            period = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['period']
            power = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['power']
            
#            ymax1 = LS_results['LS_Power1'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax2 = LS_results['LS_Power2'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax3 = LS_results['LS_Power3'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            
            plt.plot(period,power)#,c = c_sector[i])
            plt.scatter(LS_results['LS_Per1'],LS_results['LS_Power1'], c = 'r')
            plt.scatter(LS_results['LS_Per2'],LS_results['LS_Power2'], c = 'r')
            plt.scatter(LS_results['LS_Per3'],LS_results['LS_Power3'], c = 'r')
#            plt.axvline(x=LS_results['LS_Per1'].to_numpy()[0], ymin=0.01,ymax = ymax1, c= 'r', linewidth = 4)
#            plt.axvline(x=LS_results['LS_Per2'].to_numpy()[0], ymin=0.01,ymax = ymax2, c= 'r', linewidth = 2)
#            plt.axvline(x=LS_results['LS_Per3'].to_numpy()[0], ymin=0.01,ymax = ymax3, c= 'r', linewidth = 1)
            
            
        temp_Power1 = (self.LS_results['LS_Power1'].dropna()).to_numpy(dtype='float')
    
        if len(temp_Power1) > 0:        
            row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
            best_per = self.LS_results['LS_Per1'].iloc[row_idx]            
            best_sector = self.LS_results['sector'].iloc[row_idx]
        else:
            best_per = best_sector = row_idx = np.nan
        
        plt.xlim((0,15))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("LS Power")
        plt.title("LS Period = " + str(round(best_per,4)) + " - Found in Sector" + str(best_sector))
    
        ### PLOT AC PERIODOGRAMS
        
        plt.subplot2grid((3,2), (1,1),colspan=1)
        
        for key in self.acf_periodograms:
            # nan_test = sum(np.isnan(self.lc_df[self.lc_df['sector'] == sector][self.flux_type].to_numpy(dtype = 'float')))
            # flux_len = len(self.lc_df[self.lc_df['sector'] == sector])
            # if nan_test != flux_len:
            AC_results = self.acf_result[self.acf_result['sector'] == str(key)]
            period = self.acf_periodograms[key]['Period']
            power = self.acf_periodograms[key]['AC Power']
            #print(AC_results)
            plt.plot(period,power)#, c = c_sector[i])
#             plt.scatter(AC_results['AC_Per1'],AC_results['AC_Power1'], c = 'r')
#             plt.scatter(AC_results['AC_Per2'],AC_results['AC_Power2'], c = 'r')
#             plt.scatter(AC_results['AC_Per3'],AC_results['AC_Power3'], c = 'r')
            plt.axvline(x=AC_results['ac_period'][0],color = 'red', ymin=0.01, ymax = 0.99, linewidth = 4)
#            plt.axvline(x=AC_results['AC_Per2'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power2'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 2)
#            plt.axvline(x=AC_results['AC_Per3'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power3'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 1)
        
        # ac periodogram title
        ac_title = ''
        for i,key in enumerate(self.acf_periodograms):
            AC_results = self.acf_result[self.acf_result['sector'] == str(key)]
            if i>0: ac_title = ac_title + ', '
            ac_title = ac_title + 'Sector ' + str(key) + ' : ' + str(round(AC_results['ac_period'][0],2)) + ' d'
        
        if len(self.acf_result['ac_period']) >= 1: 
            assoc_AC_per = self.acf_result[self.acf_result['sector'] == best_sector]['ac_period'][0]
        else:
            assoc_AC_per = np.nan
        # if np.isnan(row_idx) == False:
        #     assoc_AC_per = self.AC_results['AC_Per1'].iloc[row_idx]
        # else:
        #     assoc_AC_per = np.nan
        
        plt.xlim((0,15))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("AC Power")
        plt.title(ac_title)
        
        ### PLOT PHASE PHOLDED CURVES FOR BEST SECTOR
        
        time_best_sector = self.lc_df[self.lc_df['sector'] == best_sector]['time'].to_numpy(dtype = 'float')
        flux_best_sector = self.lc_df[self.lc_df['sector'] == best_sector][self.flux_type].to_numpy(dtype = 'float')
        
        bin_flux,bin_time,bin_idx = bin_stat(x=time_best_sector,values=flux_best_sector,statistic = 'median', bins = round(0.05*len(time)))
        bin_time = bin_time[0:len(bin_time)-1]
        bin_flux = bin_flux/np.nanmedian(bin_flux)
        flux = flux/np.median(flux)    
        
        #normalized phases arrays
       
        if (self.flux_type == 'cpm') | (self.flux_type == 'fcor'):
            phase_norm_LS = (time_best_sector % best_per)/best_per #from LS method 
            phase_norm_AC = (time_best_sector % assoc_AC_per)/assoc_AC_per #from AC
        else:
            phase_norm_LS = (bin_time % best_per)/best_per #from LS method        
            phase_norm_AC = (bin_time % assoc_AC_per)/assoc_AC_per #from AC
        #detect changes in phase
    
        c1 = []
        for i in range(len(phase_norm_LS)-1):
            if phase_norm_LS[i+1] - phase_norm_LS[i] < 0:
                c1.append(i+1)
    
        c2 = []
        for i in range(len(phase_norm_AC)-1):
            if phase_norm_AC[i+1] - phase_norm_AC[i] < 0:
                c2.append(i+1)
        
        # LS phase fold
    
        plt.subplot2grid((3,2), (2,0),colspan=1)
        first=0
        if (self.flux_type == 'cpm') | (self.flux_type == 'fcor'):
            for change in c1:
                plt.scatter(phase_norm_LS[first:change],flux_best_sector[first:change])
                first=change+1
        else:
            for change in c1:
                plt.scatter(phase_norm_LS[first:change],bin_flux[first:change])
                first=change+1
        plt.xlabel("Fraction of Period")
        plt.ylabel("Normalized Flux")
        plt.title("LS Phase-Folded - Sector " + str(best_sector))
        
        # AC phase fold
    
        plt.subplot2grid((3,2), (2,1),colspan=1)
        first = 0
        if (self.flux_type == 'cpm') | (self.flux_type == 'fcor'):
            for change in c2:
                plt.scatter(phase_norm_AC[first:change],flux_best_sector[first:change])
                first=change+1
        else:
            for change in c2:
                plt.scatter(phase_norm_AC[first:change],bin_flux[first:change])
                first=change+1        
        plt.xlabel("Fraction of Period")
        plt.ylabel("Normalized Flux")
        plt.title("AC Phase-Folded")
        
        fig.tight_layout()
        plt.close(fig=fig)
        
        #decide on best period based on highest power.... not sure if this is great still. In fact, I know its not the end of the story...
    #    best_power_LS = ((power)[np.where(power == np.max(power))])[0]
    #    best_power_AC = ((AC_power)[np.where(AC_per == bestper_AC)])[0]
    #    
        #period_results = np.array([bestper_LS,best_power_LS,bestper_AC,best_power_AC])
        
    #     choose = np.array([best_power1,best_power2])
    #     between = np.array([bestper1,bestper2])
    #     #choice = ((between)[np.where(choose == np.argmax(choose))])[0]
        
        self.rotFig = fig
        
    def lc_graph(self, MAST = False, eleanor = False, flux_type = None):
       
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        #This stuff is everything, use it for any python plot to make it nicer.
        mpl.rcParams['lines.linewidth'] =3
        mpl.rcParams['axes.linewidth'] = 2
        mpl.rcParams['xtick.major.width'] =2
        mpl.rcParams['ytick.major.width'] =2
        mpl.rcParams['xtick.minor.width'] =1.5
        mpl.rcParams['ytick.minor.width'] =1.5
        mpl.rcParams['ytick.labelsize'] = 17
        mpl.rcParams['xtick.labelsize'] = 17
        mpl.rcParams['axes.labelsize'] = 17
        #mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['axes.labelweight']='semibold'
        mpl.rcParams['mathtext.fontset']='stix'
        mpl.rcParams['font.weight'] = 'semibold'
        mpl.rcParams['axes.titleweight']='semibold'
        mpl.rcParams['axes.titlesize']=17
        
        AC_power = self.periodogram_AC['power']
        
        AC_height = np.max(AC_power) - np.min(AC_power)
        
        #make plots
    
        fig = plt.figure(figsize=(15,12))
        gridspec.GridSpec(2,2)
        
        ### PLOT FULL LC
        plt.subplot2grid((2,2), (0,0), colspan = 2)
        
        for i,sector in enumerate(self.sector_list):
            
            time = self.lc_df[self.lc_df['sector'] == sector]['time']
            flux = self.lc_df[self.lc_df['sector'] == sector][self.flux_type]
            
            bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.05*len(time)))
            bin_time = bin_time[0:len(bin_time)-1]
                  
            if self.flux_type != 'cpm':
                bin_flux = bin_flux/np.nanmedian(bin_flux)
                flux = flux/np.median(flux) 
                plt.scatter(time,flux, s = 0.75, label = sector)
                plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
            if self.flux_type == 'cpm':
                plt.scatter(time,flux, s = 1, label = sector)#, c = c_sector[i])
                #plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)
        
        if self.flux_type == 'cpm':
            plt.ylim((-0.15,0.15))
        plt.xlabel("Time (days)")
        plt.ylabel("Normalized Flux")
        plt.legend(loc = 'upper right', fontsize = 'xx-large')
        plt.title("Light Curve of " + self.name)
        
        ### PLOT LS periodograms
    
        plt.subplot2grid((2,2), (1,0),colspan=1)
        
        for i,sector in enumerate(self.sector_list):
            LS_results = self.LS_results[self.LS_results['sector'] == sector]
            
            period = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['period']
            power = self.periodogram_LS[self.periodogram_LS['sector'] == sector]['power']
            
#            ymax1 = LS_results['LS_Power1'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax2 = LS_results['LS_Power2'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax3 = LS_results['LS_Power3'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            
            plt.plot(period,power)#,c = c_sector[i])
            plt.scatter(LS_results['LS_Per1'],LS_results['LS_Power1'], c = 'r')
            plt.scatter(LS_results['LS_Per2'],LS_results['LS_Power2'], c = 'r')
            plt.scatter(LS_results['LS_Per3'],LS_results['LS_Power3'], c = 'r')
#            plt.axvline(x=LS_results['LS_Per1'].to_numpy()[0], ymin=0.01,ymax = ymax1, c= 'r', linewidth = 4)
#            plt.axvline(x=LS_results['LS_Per2'].to_numpy()[0], ymin=0.01,ymax = ymax2, c= 'r', linewidth = 2)
#            plt.axvline(x=LS_results['LS_Per3'].to_numpy()[0], ymin=0.01,ymax = ymax3, c= 'r', linewidth = 1)
            
            
        temp_Power1 = (self.LS_results['LS_Power1'].dropna()).to_numpy(dtype='float')
    
        if len(temp_Power1) > 0:        
            row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
            best_period = self.LS_results['LS_Per1'].iloc[row_idx]            
            best_sector = self.LS_results['sector'].iloc[row_idx]
        else:
            best_period = best_sector = row_idx = np.nan
        
        plt.xlim((0,15))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("LS Power")
        plt.title("LS Period = " + str(round(best_period,4)) + " - Found in Sector" + str(best_sector))
    
        ### PLOT AC PERIODOGRAMS
        
        plt.subplot2grid((2,2), (1,1),colspan=1)
        
        for i,sector in enumerate(self.sector_list):
            AC_results = self.AC_results[self.AC_results['sector'] == sector]
            period = self.periodogram_AC[self.periodogram_AC['sector'] == sector]['period']
            power = self.periodogram_AC[self.periodogram_AC['sector'] == sector]['power']
        
            plt.plot(period,power)#, c = c_sector[i])
            plt.scatter(AC_results['AC_Per1'],AC_results['AC_Power1'], c = 'r')
            plt.scatter(AC_results['AC_Per2'],AC_results['AC_Power2'], c = 'r')
            plt.scatter(AC_results['AC_Per3'],AC_results['AC_Power3'], c = 'r')
#            plt.axvline(x=AC_results['AC_Per1'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power1'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 4)
#            plt.axvline(x=AC_results['AC_Per2'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power2'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 2)
#            plt.axvline(x=AC_results['AC_Per3'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power3'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 1)
#   
        if np.isnan(row_idx) == False:
            assoc_AC_per = self.AC_results['AC_Per1'].iloc[row_idx]
        else:
            assoc_AC_per = np.nan
        
        plt.xlim((0,15))
        #plt.ylim((0,1))
        plt.xlabel("Period (days)")
        plt.ylabel("AC Power")
        plt.title("Associated AC Period = " + str(round(assoc_AC_per,4)))
        
        self.lc_fig = fig
    
    def tab2df(self,table):
    
        df = pd.DataFrame()
        
        columns = table.colnames
        
        for col in columns:
            df[col] = table[col]
            
        return(df)
    
    
    
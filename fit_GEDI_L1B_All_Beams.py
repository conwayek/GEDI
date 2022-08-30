# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 09:34:48 2022

@author: e.conway@kri.neu.edu

Code to analyze the statistics of the GEDI LIDAR L1B .h5 files for waveform maxima/minima feature extraction

Inputs are a .h5 file

Output is a NetCDF4 file with the detections

Requirements: 
The EGM2008 dataset interpreter interp_2p5min.f is located from:/scratch/e.conway/GEDI/
Once compiled, it reads the EGM2008 parameters Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE located at /scratch/e.conway/GEDI/
It wil compile in your cwd. 
This code, below, creates the inputs.dat and the fortran code outputs the GEOID elevation above the wgs84 ellipsoid in outputs.dat
We pass the fortran code the location of inputs.dat and outputs.dat respectively

These allow us to get height above MSL.
GEOID stored in NetCDF4 file.

"""
import numba as nb
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import netCDF4 as nc4
from sys import exit
import os
import glob
import faulthandler

faulthandler.enable()

@nb.jit(nopython=True)
def func_maxima(dem_srtm_full,ystretchers_full,wfLat_full_zero,wfLon_full_zero,wfLat_full_one,wfLon_full_one,noise_full,noise_stdev_full,wv_form_full,start_idx_full,count,shot_number_full,zstart_full,zend_full):

    maxima = []
    minima = []
    elev_maxima = []
    elev_minima = []
    wv_maxima = []
    wv_minima = []
    wv_maxima_lon = []
    wv_minima_lon = []
    wv_maxima_lat = []
    wv_minima_lat = []
    integrated_wv = []
    dem_srtm = []
    
    n = len(shot_number_full)

    for i in range(n):
        nmaxima = 0
        nminima = 0
        
        shot = shot_number_full[i]
        
        index = np.where(shot==shot_number_full)[0][0]
        wfCount = count[index]
        wfStart = int(start_idx_full[index] - 1)
        
        ystretcher = ystretchers_full[index]
    
         
        # start end points of wave elevation
        zstart = zstart_full[index]
        zend = zend_full[index]
        # start end points of wave lon
        zstart_lon = wfLon_full_zero[index]
        zend_lon = wfLon_full_one[index]
        # start end points of wave lat
        zstart_lat = wfLat_full_zero[index]
        zend_lat = wfLat_full_one[index]        
        #try: 
        # subset of waveform
        wv_form = wv_form_full[wfStart: wfStart + wfCount]

        
        
        wv = wv_form    
        noise = noise_full[i]
        noise_stdev = noise_stdev_full[i]
        # noise corrected waveform
        wv = wv - noise
        if(len(wv)> 10):
            if (np.nanmax(wv) > 30):
                    
                #Integrated Energy in the noise-corrected wave
                integrated_signal = np.nansum(wv_form)
                integrated_wv.append(integrated_signal)
                
                #srtm dem
                dem_srtm.append(dem_srtm_full[int(wfCount )])
                
                # Elevation in the waveform
                stretch = np.zeros(wfCount)
                for k in range(1,wfCount+1):
                    stretch[k-1] = np.add(zend,np.multiply(k,((zstart - zend) / int(wfCount))))
                stretch = stretch[::-1]
                
                # Latitude Points in the waveform
                stretch_lon = np.zeros(wfCount)
                for k in range(1,wfCount+1):
                    stretch_lon[k-1] = np.add(zend_lon,np.multiply(k,((zstart_lon - zend_lon) / int(wfCount))))
                stretch_lon = stretch_lon[::-1]        
                
                # Longitude in the waveform
                stretch_lat = np.zeros(wfCount)
                for k in range(1,wfCount+1):
                    stretch_lat[k-1] = np.add(zend_lat,np.multiply(k,((zstart_lat - zend_lat) / int(wfCount))))
                stretch_lat = stretch_lat[::-1]      
                

            
                idx = np.where(wv<noise_stdev)
                wv[idx] = np.nan
            
                CompleteSearch = False
                start_idx = 0
                detected = False
                
                while CompleteSearch == False: 

                    if(np.isfinite(wv[start_idx]) and detected == False):
                        if(wv[start_idx] > 5*noise_stdev):
                            nminima = nminima + 1
                            elev_minima.append(stretch[start_idx])
                            wv_minima.append(wv[start_idx])
                            wv_minima_lon.append(stretch_lon[start_idx])
                            wv_minima_lat.append(stretch_lat[start_idx])
                        for j in range(start_idx,len(wv)-1):
                            if(j >= (len(wv) - 2) and detected == False):
                                CompleteSearch = True
                            if(wv[j+1]<wv[j] and detected == False and wv[j] > 5*noise_stdev):
                                nmaxima = nmaxima + 1
                                elev_maxima.append(stretch[int(j)])
                                wv_maxima.append(wv[int(j)])
                                wv_maxima_lon.append(stretch_lon[int(j)])
                                wv_maxima_lat.append(stretch_lat[int(j)])          
                                start_idx = j 
                                detected = True
                             
                    elif(detected == True and np.isfinite(wv[start_idx])):  
                        for j in range(start_idx,len(wv)-1):
                            if(j >= (len(wv) - 2) and detected == True):
                                CompleteSearch = True
                            if(wv[j+1]>wv[j] and detected == True):
                                start_idx = j
                                detected = False
                    elif(np.nansum(wv[start_idx:]) == np.nan):
                        CompleteSearch = True    
                    elif(np.isfinite(wv[start_idx]!=True)):
                        start_idx = start_idx + 1
                        if(start_idx >= (len(wv) - 2)):
                            CompleteSearch = True
                            
                minima.append(nminima)
                maxima.append(nmaxima)
            else:
                integrated_wv.append(-1)
                minima.append(0)
                maxima.append(0)
                dem_srtm.append(np.nan)
        else:
            integrated_wv.append(-1)
            minima.append(0)
            maxima.append(0)
            dem_srtm.append(np.nan)
    return(dem_srtm,wv_minima_lon,wv_minima_lat,wv_maxima_lon,wv_maxima_lat,elev_maxima,elev_minima,minima,maxima,wv_maxima,wv_minima,integrated_wv)


l1bdir = '/scratch/e.conway/GEDI/L1B/'
outdir = '/scratch/e.conway/GEDI/MyResults/'

files = glob.glob(l1bdir+'*.h5')
beams = ['BEAM0000','BEAM0001','BEAM0010','BEAM0011','BEAM0101','BEAM0110','BEAM1000','BEAM1011']
for file in files:
    fnew = file.split('.h5')[0] + '.nc'
    fnew = os.path.basename(fnew)
    fnew=os.path.join(outdir,fnew)
    if(os.path.exists(fnew)==False):
        for beam in beams:
        
        
            df = h5.File(file,'r')
            count = df[beam]['rx_sample_count'][()]
            start_idx_full = df[beam]['rx_sample_start_index'][()]
            shot_number_full = df[beam]['shot_number'][()]
            wfLat_full_zero = df[beam]['geolocation/latitude_bin0'][()]
            wfLon_full_zero = df[beam]['geolocation/longitude_bin0'][()]
            wfLat_full_one = df[beam]['geolocation/latitude_lastbin'][()]
            wfLon_full_one = df[beam]['geolocation/longitude_lastbin'][()]
            zstart_full = df[beam]['geolocation/elevation_bin0'][()]
            zend_full = df[beam]['geolocation/elevation_lastbin'][()]
            wv_form_full = df[beam]['rxwaveform'][()]
            noise_stdev_full=df[beam]['noise_stddev_corrected'][()]
            noise_full=df[beam]['noise_mean_corrected'][()]
            master_int=df[beam]['master_int'][()]
            master_frac=df[beam]['master_frac'][()]
            dem_full=df[beam]['geolocation/digital_elevation_model'][()]
            dem_srtm_full = df[beam]['geolocation/digital_elevation_model_srtm'][()]
            ystretchers_full = df[beam]['selection_stretchers_y'][()]
            df.close()
            
            time = master_int + master_frac
        
            # Subsetting based on lon/lat sizes
            
                
            dem_srtm,lon_min,lat_min,lon_max,lat_max,elev_maxima,elev_minima,minima,maxima,wv_maxima,wv_minima,integrated_signal = func_maxima(dem_srtm_full,ystretchers_full,wfLat_full_zero,wfLon_full_zero,wfLat_full_one,wfLon_full_one,noise_full,noise_stdev_full,wv_form_full,start_idx_full,count,shot_number_full,zstart_full,zend_full)
            elev_maxima = np.array(elev_maxima)           
            elev_minima = np.array(elev_minima)           
            wv_maxima = np.array(wv_maxima)           
            wv_minima = np.array(wv_minima)   
            time = np.array(time)      
            integrated_signal = np.array(integrated_signal)
            lon_min=np.array(lon_min)
            lat_min=np.array(lat_min)
            lon_max=np.array(lon_max)
            lat_max=np.array(lat_max)
            
            # Let us get the GEOID data       


            inputdat = 'INPUT.DAT'
            with open(inputdat,'w') as f:
                for i in range(len(wfLon_full_zero)):
                    if(wfLat_full_zero[i] == np.nan or wfLon_full_zero[i] == np.nan):
                        f.write(str(0.0)+' '+str(0.0)+'\n')
                    else:
                        f.write(str(wfLat_full_zero[i])+' '+str(wfLon_full_zero[i])+'\n')
            
            os.system("gfortran /scratch/e.conway/GEDI/interp_2p5min.f -o interp_2p5min.x")
            os.system("./interp_2p5min.x "+str(os.getcwd()+'/ ')+str(os.getcwd()+'/ '))
            
            outputdat = 'OUTPUT.DAT'
            with open(outputdat,'r') as f:
                data = f.readlines()
            
            add = 0
            geoid = np.zeros(len( wfLon_full_zero))
            for i in range(len(data)):
                if(wfLat_full_zero[i+add] == np.nan or wfLon_full_zero[i+add] == np.nan):
                    geoid[i+add] = np.nan
                    add = add + 1
                    geoid[i+add] = data[i][24:]
                else:
                    geoid[i+add] = data[i][24:]

  
            
            
            npts = len(master_int)
            n_elev_minima= len(elev_minima)
            n_elev_maxima= len(elev_maxima)
            
            

            if(beam=='BEAM0000'):
                mode = 'w'
            else:
                mode = 'a'
            with nc4.Dataset(fnew,mode,format='NETCDF4') as f:
        
                newgrp = f.createGroup(str(beam))
                
                newgrp.createDimension('Npts',npts)
                newgrp.createDimension('Nminima',n_elev_minima)
                newgrp.createDimension('Nmaxima',n_elev_maxima)
                
                time_nc = newgrp.createVariable('Time','f8','Npts',zlib=True)
                #lon_max_nc = newgrp.createVariable('LonMaxima','f8','Nmaxima',zlib=True)
                #lat_max_nc = newgrp.createVariable('LatMaxima','f8','Nmaxima',zlib=True)
                #lon_min_nc = newgrp.createVariable('LonMinima','f8','Nminima',zlib=True)
                #lat_min_nc = newgrp.createVariable('LatMinima','f8','Nminima',zlib=True)
                lon_max_nc = newgrp.createVariable('Lon','f8','Npts',zlib=True)
                lat_max_nc = newgrp.createVariable('Lat','f8','Npts',zlib=True)
                geoid_nc = newgrp.createVariable('GEOID','f8','Npts',zlib=True)
                dem_srtm_nc = newgrp.createVariable('DEM_Srtm','f4','Npts',zlib=True)
                dem_x90_nc = newgrp.createVariable('DEM_x90','f4','Npts',zlib=True)
                elev_maxima_nc = newgrp.createVariable('MaximaElevation','f4','Nmaxima',zlib=True)
                elev_minima_nc = newgrp.createVariable('MinimaElevation','f4','Nminima',zlib=True)
                wv_maxima_nc = newgrp.createVariable('WaveformSignalMax','f4','Nmaxima',zlib=True)
                wv_minima_nc = newgrp.createVariable('WaveformSignalMin','f4','Nminima',zlib=True)
                integrated = newgrp.createVariable('IntegratedSignal','f4','Npts',zlib=True)
                maxima_nc = newgrp.createVariable('MaximaCounter','f4','Npts',zlib=True)
                minima_nc = newgrp.createVariable('MinimaCounter','f4','Npts',zlib=True)
                
                #lon_max_nc[:] = lon_max[:]
                #lon_min_nc[:] = lon_min[:]
                #lat_max_nc[:] = lat_max[:]
                #lat_min_nc[:] = lat_min[:]
                lon_max_nc[:] = wfLon_full_zero[:]
                lat_max_nc[:] = wfLat_full_zero[:]
                geoid_nc[:] = geoid[:]          
 
                dem_srtm_nc[:] = dem_srtm[:]
                dem_x90_nc[:] = dem_full[:]
                elev_maxima_nc[:] = elev_maxima[:]
                elev_minima_nc[:] = elev_minima[:]
                wv_maxima_nc[:] = wv_maxima[:]
                wv_minima_nc[:] = wv_minima[:]
                integrated[:] = integrated_signal[:]
                maxima_nc[:] = maxima[:]
                minima_nc[:] = minima[:]
                time_nc[:] = time[:]
            
        

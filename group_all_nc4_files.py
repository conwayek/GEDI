# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 15:32:28 2022

@author: e.conway@kri.neu.edu

Code to analyze the data of the GEDI LIDAR L1B NetCDF4 files for waveform maxima/minima feature extraction
produced by the fit_gedi_v3.py code

Input is a folder location to group all .nc files together than have a lon/lat point within our defined/input polygon

Output is a .nc file with data collated

"""

import netCDF4 as nc4
import numpy as np
import numba as nb
from sys import exit
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import make_kmz_v2
import os
from osgeo import osr
import glob

parameters = {'axes.labelsize': 40,
          'axes.titlesize': 40,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'legend.fontsize': 30}
plt.rcParams.update(parameters)

# Make a kml of all results?
# Warning:: It may be very large
make_kml = False

#@nb.jit(fastmath=True)#nopython=True)
def listlists(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11):
    n1=[]
    n2=[]
    n3=[]
    n4=[]
    n5=[]
    n6=[]
    n7=[]
    n8=[]
    n9=[]
    n10=[]
    n11=[]
    n1 = [t1 for sublist in list1 for t1 in sublist ]
    n2 = [t2 for sublist in list2 for t2 in sublist ]
    n3 = [t3 for sublist in list3 for t3 in sublist ]
    n4 = [t4 for sublist in list4 for t4 in sublist ]
    n5 = [t5 for sublist in list5 for t5 in sublist ]
    n6 = [t6 for sublist in list6 for t6 in sublist ]
    n7 = [t7 for sublist in list7 for t7 in sublist ]
    n8 = [t8 for sublist in list8 for t8 in sublist ]
    n9 = [t9 for sublist in list9 for t9 in sublist ]
    n10 =[t10 for sublist in list10 for t10 in sublist ]
    n11 =[t11 for sublist in list11 for t11 in sublist ]
    return(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11)

myresults = '/scratch/e.conway/GEDI/MyResults/'
outdir = '/scratch/e.conway/GEDI/MyResults/'

#files=glob.glob(myresults+'*O02458_04_T00568_02_005_01_V002.nc')
files = glob.glob(myresults+'GEDI01_B_2020*'+'*.nc')


beams = ['BEAM0000','BEAM0010','BEAM0001','BEAM0011','BEAM0101','BEAM0110','BEAM1000','BEAM1011']

dic = {}
count=0
for beam in beams:
    count=count+1
    dic[beam] = count
time=[]
geoid = []
maxi = []
maxi_lat = []
maxi_lon = []
dem = []
beam_number = []

# for all the strong points
nmax_elev = []
nmax_elev_beam = []
nmax_elev_points = []
nmax_elev_lon = [] 
nmax_elev_lat = [] 

for file in files:
    for beam in beams:
        maxi_zero = []
        geoid_zero = []
        maxi_lat_zero = []
        maxi_lon_zero = []
        dem_zero = []
        tzero = []
        bzero = []
        nmax_elev_zero = []
        nmax_elev_beam_zero = []
        nmax_elev_points_zero = []
        nmax_elev_lon_zero = [] 
        nmax_elev_lat_zero = [] 
        with nc4.Dataset(file,'r') as f:
    
    
            ftime = f.groups[beam].variables['Time'][:]     
            max_elev = f.groups[beam].variables['MaximaElevation'][:]
            min_elev = f.groups[beam].variables['MinimaElevation'][:]
            max_wf = f.groups[beam].variables['WaveformSignalMax'][:]
            min_wf = f.groups[beam].variables['WaveformSignalMin'][:]
            int_sig = f.groups[beam].variables['IntegratedSignal'][:]
            max_count = f.groups[beam].variables['MaximaCounter'][:]
            min_count = f.groups[beam].variables['MinimaCounter'][:]
            lon = f.groups[beam].variables['Lon'][:]
            lat = f.groups[beam].variables['Lat'][:]
            dem_srtm = f.groups[beam].variables['DEM_Srtm'][:]
            dem_x90=f.groups[beam].variables['DEM_x90'][:]
            geoid_raw=f.groups[beam].variables['GEOID'][:]
      
        max_counter = 0
        npts = len(time)
        
    
        dem_choice = dem_x90
        
        # routine for pulling the maximum 'max elevation' point from our previous waveform analysis.
        # we want the lon/lat too.
        # elevation is w.r.t WGS84 Ellipsoid
        nmax = len(lon)  
        for i in range(int(nmax)):
            if(int_sig[i] != -1 and max_count[i] != 0):
                nmaxpnt = max_count[i]
                maxi_zero.append( np.nanmax(max_elev[int(max_counter):int(max_counter+nmaxpnt)]))
                geoid_zero.append(geoid_raw[i])
                # get the lon + lat's
                maxi_lat_zero.append(lat[i])
                maxi_lon_zero.append(lon[i])
                dem_zero.append(dem_choice[i])
                tzero.append(ftime[i])
                bzero.append(dic[beam])

                if(nmaxpnt > 1):
                    # here, let us pull out all the strongest features including the minima from the data
                    elev_maxima = max_elev[int(max_counter):int(max_counter+nmaxpnt)]
                    wv_maxima = max_wf[int(max_counter):int(max_counter+nmaxpnt)]
                    # standard deviation in the noise is consistently approximately 4.0
                    noise_stdev = 4.0
                    mx_signal = np.nanmax(wv_maxima)
                    thresh = 0.10 * mx_signal
                    # kick out any weak signals
                    idx = np.where((wv_maxima>thresh) & (wv_maxima>3.0*noise_stdev))
                    elev_maxima=elev_maxima[idx] 
                    wv_maxima=wv_maxima[idx] 
                    for k in range(len(wv_maxima)-1):
                        if(np.isfinite(elev_maxima[k]) and np.isfinite(elev_maxima[k+1])):
                            if(abs(elev_maxima[k] - elev_maxima[k+1]) <= 2.0):
                                idx_max = max(wv_maxima[k],wv_maxima[k+1])
                                wx_max = np.where(wv_maxima==idx_max)
                                elev_maxima[wx_max] = np.nan
                                wv_maxima[wx_max] = np.nan
                    ct = np.nansum(np.isfinite(elev_maxima))
                    # split lon/lat into ct points
                    tlat = np.linspace(lat[i]*0.99985,lat[i]*1.00015,ct)
                    tlon = np.linspace(lon[i]*0.99985,lon[i]*1.00015,ct)
                    nmax_elev_points_zero.append(ct)
                    iterable=0
                    for k in range(len(wv_maxima)):
                        if(np.isfinite(elev_maxima[k])):
                            nmax_elev_beam_zero.append(dic[beam])
                            nmax_elev_zero.append(elev_maxima[k])
                            nmax_elev_lon_zero.append(tlon[iterable])
                            nmax_elev_lat_zero.append(tlat[iterable])
                            iterable=iterable+1
                else:
                    nmax_elev_points_zero.append(0)                   
                    nmax_elev_beam_zero.append(np.nan)
                    nmax_elev_zero.append(np.nan)
                    nmax_elev_lon_zero.append(np.nan)
                    nmax_elev_lat_zero.append(np.nan)
                max_counter = max_counter + nmaxpnt
            else:
                maxi_zero.append(np.nan)
                maxi_lat_zero.append(np.nan)
                maxi_lon_zero.append(np.nan)
                dem_zero.append(np.nan)
                geoid_zero.append(np.nan)
                tzero.append(ftime[i])
                bzero.append(dic[beam])
     
        maxi.append(maxi_zero)
        maxi_lat.append(maxi_lat_zero)
        maxi_lon.append(maxi_lon_zero)
        dem.append(dem_zero)    
        geoid.append(geoid_zero)    
        beam_number.append(bzero) 
        time.append(tzero)

        nmax_elev.append(nmax_elev_zero)
        nmax_elev_beam.append(nmax_elev_beam_zero)
        nmax_elev_points.append(nmax_elev_points_zero)
        nmax_elev_lat.append(nmax_elev_lat_zero)
        nmax_elev_lon.append(nmax_elev_lon_zero)
    


maxi,maxi_lat,maxi_lon,geoid,beam_number,time,nmax_elev,nmax_elev_beam,nmax_elev_points,nmax_elev_lon,nmax_elev_lat = listlists(maxi,maxi_lat,maxi_lon,geoid,beam_number,time,nmax_elev,nmax_elev_beam,nmax_elev_points,nmax_elev_lon,nmax_elev_lat)



time = np.array(time)
maxi = np.array(maxi)
maxi_lat=np.array(maxi_lat)
maxi_lon=np.array(maxi_lon)
geoid = np.array(geoid)
beam_number = np.array(beam_number)
nmax_elev=np.array(nmax_elev)
nmax_elev_beam=np.array(nmax_elev_beam)
nmax_elev_points=np.array(nmax_elev_points)
nmax_elev_lon=np.array(nmax_elev_lon)
nmax_elev_lat=np.array(nmax_elev_lat)

msl = maxi - geoid

idx = np.where(np.isfinite(maxi)&np.isfinite(maxi_lon))
time = time[idx]
maxi = maxi[idx]
maxi_lat=maxi_lat[idx]
maxi_lon=maxi_lon[idx]
geoid = geoid[idx]
beam_number = beam_number[idx]
msl = msl[idx]



fname = 'Combined_Congo_2020_South.nc'
fname = os.path.join(outdir,fname)

npts = len(maxi)

nmaxima = len(nmax_elev)

with nc4.Dataset(fname,'w',format='NETCDF4') as f:
    f.createDimension('Npts',npts)
    time_nc = f.createVariable('Time','f8','Npts',zlib=True)
    lon_nc = f.createVariable('Lon','f4','Npts',zlib=True)
    lat_nc = f.createVariable('Lat','f4','Npts',zlib=True)
    elev_msl = f.createVariable('ElevMSL','f4','Npts',zlib=True)
    elev_wgs = f.createVariable('ElevWGS84','f4','Npts',zlib=True)
    beam_nc = f.createVariable('BeamNumber','i2','Npts',zlib=True)

    f.createDimension('Nmaxima',nmaxima)
    nmax_elev_lat_nc = f.createVariable('MaximaLat','f4','Nmaxima',zlib=True) 
    nmax_elev_lon_nc = f.createVariable('MaximaLon','f4','Nmaxima',zlib=True) 
    nmax_elev_nc = f.createVariable('MaximaElev','f4','Nmaxima',zlib=True) 
    nmax_elev_beam_nc = f.createVariable('MaximaBeam','i2','Nmaxima',zlib=True)
    nmax_elev_points_nc = f.createVariable('NpointsPerWv','i2','Npts',zlib=True)

    time_nc[:] = time
    lon_nc[:] = maxi_lon
    lat_nc[:] = maxi_lat
    elev_wgs[:] = maxi
    elev_msl[:] = msl
    beam_nc[:] = beam_number
    nmax_elev_nc[:] =  nmax_elev
    nmax_elev_beam_nc[:] = nmax_elev_beam
    nmax_elev_points_nc[:] = nmax_elev_points
    nmax_elev_lat_nc[:] =  nmax_elev_lat
    nmax_elev_lon_nc[:] =  nmax_elev_lon

if(make_kml == True):
    color = 255*maxi.flatten() / vmax
    idx = np.where(np.isfinite(maxi) & np.isfinite(maxi_lon))[0]
    fname = os.path.basename(file).split('nc')[0]
    fname=fname+'kml'
    fname = os.path.join(outdir,fname)
    
    msl = msl.flatten()[idx]
    maxi = maxi.flatten()[idx]
    maxi_lat = maxi_lat.flatten()[idx]
    maxi_lon=maxi_lon.flatten()[idx]
    color = color[idx]
    
    make_kmz_v2.make_kml(maxi_lon,maxi_lat,msl,color,fname)
    #make_kmz_v2.make_kml(maxi_lon,maxi_lat,maxi,color,fname)



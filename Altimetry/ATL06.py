import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd
import os
import h5py
from shapely.geometry import point,polygon
import icepyx as ipd
from datetime import date
from dateutil.relativedelta import relativedelta
from mpl_toolkits.axes_grid1 import make_axes_locatable
from playsound import playsound

def get_data(bbox,date_range,path) :
    try:
        os.mkdir(path+'/'+date_range[0]+'--'+date_range[1])
    except:
        None
    path = path+'/'+date_range[0]+'--'+date_range[1]
    #creating the icepyx object
   
    region = ipd.Query('ATL06',bbox,date_range)
    print(region.avail_granules())
    region.granules.avail

    #logging into earthdata

    earthdata_uid = input("Enter your Earthdata username:")
    email = input("Enter your Eathdata email:")
    region.earthdata_login(earthdata_uid,email)

    #creating a default variable list

    region.order_vars.append(defaults=True)
    #print(region.order_vars.wanted,sep='/n')
    region.order_vars.remove(all=True)

    #modifying the default variable list

    #print(region.order_vars.wanted)
    region.order_vars.append(var_list=['latitude'])
    region.order_vars.append(var_list=['longitude'])
    region.order_vars.append(var_list=['h_li'])
    region.order_vars.append(var_list=['x_atc'])
    region.order_vars.append(var_list=['atl06_quality_summary'])
    print("The requested data is:")
    print(region.order_vars.wanted)
   
   
    
    region.subsetparams(Coverage=region.order_vars.wanted)
    region.reqparams['page_size']=int(input("Enter desired number of granules per order:"))

    #ordering data 
    email=input("Do you want an email containing information of your order requests(y/n)")
    email=True if email=='y' else False
    region.order_granules(email=email)

    #downloading data

    region.download_granules(path)

def data_to_csv(path_in):
    group = ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']
    try:
        os.mkdir(path_in+'/CSV') 
    except:
        None
    path_out = path_in+'/CSV'
    a=os.listdir(path_in)
    try:
        a.remove('.ipynb_checkpoints')
    except:
        None
    for g in group:
        beam = pd.DataFrame()
        beam['lat']=[]
        beam['lon']=[]
        beam['h_li']=[]
        beam['x_atc']=[]
        beam['q_flag']=[]

        for fname in a:
            df = pd.DataFrame()
            fname = path_in+'/'+fname
            try:
                with h5py.File(fname,'r') as f:
                    try:
                        df['lat'] = f['/'+g+'/land_ice_segments/latitude'][:]
                        df['lon'] = f['/'+g+'/land_ice_segments/longitude'][:]
                        df['h_li'] = f['/'+g+'/land_ice_segments/h_li'][:]
                        df['x_atc'] = f['/'+g+'/land_ice_segments/ground_track/x_atc'][:]
                        df['q_flag'] = f['/'+g+'/land_ice_segments/atl06_quality_summary'][:]
                        beam=beam.append(df,ignore_index=True)
                    except:
                        print(fname+" has no relevant data")
                        continue
            except:
                print(fname+" is not a hdf5 file.")
                continue
        beam=beam[beam['h_li']< 8611]
        beam.to_csv(path_out+'/'+g+'.csv')


def h_li_plot(region,end_time):
    year, month, day = map(int, end_time.split('-'))
    start_time = date(year, month, day)+relativedelta(months=-3)
    start_time = str(start_time)
    print(start_time)
    date_range=[start_time,end_time]
    if region in ['Karakoram','West Himalaya','East Himalaya','Central Himalaya']:
        #Data download
        try:
            os.mkdir(os.getcwd().rstrip('\package') + '/ISDAT/'+region)
        except:
            None
        try:
            os.mkdir(os.getcwd().rstrip('\package')+'/'+region+'/data')
        except:
            None
        basemap = gpd.read_file(os.getcwd().rstrip('\package')+'/'+region+'/shapefile/'+region+'.shp')
        if (os.path.isdir(os.getcwd().rstrip('\package')+'/'+region+'/data/'+date_range[0]+'--'+date_range[1])==False or len(os.listdir(os.getcwd().rstrip('\package')+'/'+region+'/data/'+date_range[0]+'--'+date_range[1]))<=1):
            path = os.getcwd().rstrip('\package')+'/'+region+'/data'
            fname = os.getcwd().rstrip('\package')+'/'+region+'/shpfile/'+region+'.shp'
            print("Downloding data")
            get_data(fname,date_range,path)
            data_to_csv(path+'/'+date_range[0]+'--'+date_range[1])
        else:
            path = os.getcwd().rstrip('\package')+'/'+region+'/data/'+date_range[0]+'--'+date_range[1]
            fname = os.getcwd().rstrip('\package')+'/'+region+'/shpfile/'+region+'.shp'
            print("Data aldready exists")
        df1=pd.read_csv(path+'/CSV/gt1l.csv')
        df2=pd.read_csv(path+'/CSV/gt2l.csv')
        df3=pd.read_csv(path+'/CSV/gt3l.csv')    
        beam_gdf1=gpd.GeoDataFrame(df1,geometry=gpd.points_from_xy(df1['lon'],df1['lat']),crs='EPSG:4326')
        beam_gdf2=gpd.GeoDataFrame(df2,geometry=gpd.points_from_xy(df2['lon'],df2['lat']),crs='EPSG:4326')
        beam_gdf3=gpd.GeoDataFrame(df3,geometry=gpd.points_from_xy(df3['lon'],df3['lat']),crs='EPSG:4326')
        fig, ax = plt.subplots(1, 1)
        ax.set_facecolor('black')
        divider = make_axes_locatable(ax)
        cax_ = divider.append_axes("right", size="5%", pad=0.05)
        cat_= divider.append_axes("bottom", size="5%", pad=0.5)
        cat_.plot(date_range,np.zeros_like(date_range))
        basemap.plot(ax=ax,edgecolor='y',color='gainsboro')
        beam_gdf1.plot('h_li',cmap='BuPu',markersize=0.1,ax=ax,legend=True,cax=cax_,label="land",vmin=3000)
        beam_gdf2.plot('h_li',ax=ax,cmap='BuPu',markersize=0.1,vmin=3000)
        beam_gdf3.plot('h_li',ax=ax,cmap='BuPu',markersize=0.1,vmin=3000)
        ax.set_title("Land Ice Height (CRS:WGS 84 EPSG:4326)")
        left, width = .15,.5
        bottom, height = .1, .5
        right = left + width
        top = bottom + height
        ax.text(left,bottom,region,horizontalalignment='left', verticalalignment='center',transform=ax.transAxes,bbox=dict(facecolor='gainsboro',alpha=0.3))
        playsound(os.getcwd().rstrip('\package')+'/mixkit-musical-reveal-961.wav')
    else:
        print('Enter correct region') 

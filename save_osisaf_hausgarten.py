#!/usr/bin/env python
# coding: utf-8

# ### Get OSI SAF SIC at Hausgarten stations

# #### Modules and functions

# In[3]:


import numpy as np
import pylab as plt
import xarray as xr
import os,sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
import pandas as pd
import datetime
import glob


# In[4]:


def convert_to_decimal_degrees(coord_str):
    # Remove the degree symbol and split the string
    degree_str, minute_str_with_direction = coord_str.split('°')
    minute_str, direction = minute_str_with_direction.split("'")

    # Remove any commas from the minute string and convert to float
    minutes = float(minute_str.replace(',', '.'))

    # Convert degree string to float
    degrees = float(degree_str)

    # Calculate decimal degrees
    decimal_degrees = degrees + minutes / 60

    # If the direction is 'S' or 'W', make the decimal degrees negative
    if direction.strip() in ['S', 'W']:
        decimal_degrees = -decimal_degrees

    return decimal_degrees


# In[5]:


def xydist(x,y,x_2d,y_2d):
    return np.sqrt((x-x_2d)**2 + (y-y_2d)**2)


# In[6]:


def XY_To_LatLon(x,y,projstring):
    p1=pyproj.Proj(projstring,preserve_units=True)
    (lat,lon)=p1(x,y,inverse=True,)
    return(lat,lon)


# In[7]:


def LatLon_To_XY(lat,lon,projstring):
    p1=pyproj.Proj(projstring,preserve_units=True)
    (x,y)=p1(lat,lon,inverse=False,)
    return(x,y)


# In[8]:


def get_coorddicts(df):
    londict = {stationname:np.round(lon,3) for stationname,lon in zip(df["Station ID"],df["lon_decimal"])}
    latdict = {stationname:np.round(lat,3) for stationname,lat in zip(df["Station ID"],df["lat_decimal"])}
    return londict,latdict


# In[9]:


def design_map(ax,title,hemisphere):
    if hemisphere == "nh":
        ax.set_extent([-30, 30, 75,85], ccrs.PlateCarree())
    elif hemisphere == "sh":
        ax.set_extent([-180, 180, -60,-90], ccrs.PlateCarree())
    else:
        sys.exit("Hemisphere needs to be either 'nh' or 'sh'!")
    ax.set_title(title,fontsize = 14) # set title
    ax.coastlines()
    ax.gridlines()
    #ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='#F5E5A0')  # Adjust facecolor as needed
    #ax.set_facecolor('#F5F5F5')
    #ax.gridlines().label_style = {'fontsize': 14}
    #gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, alpha=0.5, linestyle='--', draw_labels=True)
    #gl.xlabel_style = {'size': 14, 'color': 'black'}
    #gl.ylabel_style = {'size': 14, 'color': 'black'}


# #### Load data

# In[10]:


datadir_osisaf = "/albedo/work/user/vludwig/05_SINXS/02_DATA/00_MISC/osisaf/reduced/"
datadir_hausgarten = "/albedo/work/user/vludwig/06_HAUSGARTEN/02_DATA"
plotdir = "/albedo/home/vludwig/06_HAUSGARTEN/03_PLOTS/"


# ##### Filenames of OSI SAF data

# In[12]:


fns = list(np.sort([os.path.basename(fn) for fn in glob.glob(os.path.join(datadir_osisaf,"ice_conc_nh*nc"))]))


# #### Get coordinates

# ##### OSI SAF coordinates

# In[13]:


coordfile = xr.open_dataset(os.path.join(datadir_osisaf,fns[0]))
x_2d,y_2d = np.meshgrid(coordfile.coords["xc"],coordfile.coords["yc"])


# In[14]:


lon_2d,lat_2d = XY_To_LatLon(x_2d,y_2d,projstring = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")


# ##### Hausgarten coordinates

# In[15]:


df_hausgarten = pd.read_excel(os.path.join(datadir_hausgarten,"average HG-Stations.xlsx"))


# In[16]:


# Define the list of Station IDs you want to filter by
station_ids = ['EG-I', 'EG-II', 'EG-III', 'EG-IV', 'N3', 'N4', 'N5']

# Select rows where "Station ID" is in the specified list
hausgarten_ice = df_hausgarten[df_hausgarten['Station ID'].isin(station_ids)]
inds = hausgarten_ice.index.to_list()


# In[17]:


hausgarten_ice["lat_decimal"] = [convert_to_decimal_degrees(coordstr) for coordstr in hausgarten_ice["average lat"]]
hausgarten_ice["lon_decimal"] = [convert_to_decimal_degrees(coordstr) for coordstr in hausgarten_ice["average long"]]
hausgarten_ice["color"] = plt.rcParams['axes.prop_cycle'].by_key()['color'][0:hausgarten_ice.shape[0]]


# In[18]:


londict,latdict = get_coorddicts(hausgarten_ice)


# #### Plot map

# In[20]:


if False:
    fig, ax = plt.subplots(nrows=1,ncols=1,subplot_kw={'projection': ccrs.LambertAzimuthalEqualArea(central_longitude=0, central_latitude=90)},figsize=(8,8))

    design_map(ax,title = f"Hausgarten stations", hemisphere="nh")
    
    
    #im2 = ax.scatter(hausgarten_ice["lon_decimal"], hausgarten_ice["lat_decimal"],c = hausgarten_ice["color"].to_list(),label = hausgarten_ice["Station ID"].to_list(), transform = ccrs.PlateCarree())
    for i in [16,17,18,19,21,22,23]:
        ax.scatter(hausgarten_ice['lon_decimal'][i], hausgarten_ice['lat_decimal'][i], label=hausgarten_ice['Station ID'][i], transform = ccrs.PlateCarree())
    
    # Add legend (one entry per unique label)
    ax.legend(ncol = 3)
    
    
    fig.savefig(os.path.join(plotdir,f"stations_hausgarten.pdf"),bbox_inches = "tight", dpi = 150)
    
    fig.show()
    plt.close(fig)
    #sys.exit()


# #### Get data

# In[30]:


# Define years and months
years = list(range(1999, 2024))
months = [f"{d:02d}" for d in range(1,13)]

# Generate station information
station_ids = hausgarten_ice["Station ID"]

projstring = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

# Create base columns: Station ID, lon, lat, sic, year, and month
for station_id in station_ids:
    print(f"### {station_id} ###")
    lon_station = londict[station_id]
    lat_station = latdict[station_id]
    x_station,y_station = LatLon_To_XY(lon_station,lat_station,projstring)
    dist_station = xydist(x_station, y_station, x_2d,y_2d)
    minind_station = np.where(dist_station == dist_station.min())
    data_station = {
        "station_id":np.repeat(station_id,len(years)*12),
        'lon': londict[station_id],
        'lat': latdict[station_id],
        'year': np.repeat(years,len(months)),
        'month': months*len(years),

    }
    
    # Create the DataFrame
    df_station = pd.DataFrame(data_station)
    sic_mean_station,sic_std_station = [],[]
    for year in years:
        print(f"### {year} ###")
        print(f"Time: {datetime.datetime.strftime(datetime.datetime.now(),'%H:%M:%S')}")
        if year <=2020:
            prefix = "ice_conc_nh_ease2-250_cdr-v3p0_"
        else:
            prefix = "ice_conc_nh_ease2-250_icdr-v3p0_"
        print("Loading data...")
        ds_year = xr.open_mfdataset(f"{datadir_osisaf}{prefix}{year}*nc")
        print("Extracting station...")
        ds_point = ds_year.isel(xc =minind_station[1], yc = minind_station[0])
        print("Getting mean and standard deviation...")
        monthly_mean = ds_point['sea_ice_concentration'].resample(time='1ME').mean()
        monthly_std = ds_point['sea_ice_concentration'].resample(time='1ME').std()
        print("Adding to dataframe...")
        df_station.loc[df_station["year"] == year,"sic_mean"] = np.round(monthly_mean.to_numpy().flatten(),3)
        df_station.loc[df_station["year"] == year,"sic_std"] = np.round(monthly_std.to_numpy().flatten(),3)
    print(f"Saving dataframe for {station_id}...")
    df_station.to_csv(os.path.join(datadir_hausgarten,f"hausgarten_sic_mean_std_{station_id}_{years[0]}-{years[-1]}.csv"),index = False)

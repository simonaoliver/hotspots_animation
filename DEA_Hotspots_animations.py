#!/usr/bin/env python
# coding: utf-8

# Loads DEA Hotspots data for a given time and location, and animates the data over an image underlay with fading colours that represent the age of the hotspot.


##################
# Import modules #
##################

import os
import logging as logger
logger.basicConfig(format='%(levelname)s:%(message)s', level=logger.INFO)
import pandas as pd
import geopandas as gpd
import datetime as dt
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser()
from jinja2 import Environment, FileSystemLoader
import yaml
import xarray as xr
from owslib.wms import WebMapService
from matplotlib.cm import get_cmap

# Create custom cmap with dark grey at end 
parser.add_argument('--configuration', dest='configuration', default='config.yaml',help='animation configuration')
args = parser.parse_args()

# Set parameters used to load and visualise DEA Hotspots data

################################
# Load and clean hotspots data #
################################
def filter(sensors):

    for sensordict in sensors:
        
        filter_string = ''
        count = 0
        
        for sensor in sensordict.keys():
            filter_string = filter_string+'(sensor=%27'+sensor+'%27%20AND%20(product=%27'
            product_count = 0
            for product in sensordict[sensor]:
                filter_string = filter_string+product+'%27'
                if product_count < (len(sensordict[sensor])-1):
                    filter_string = filter_string+'%20OR%20product=%27'
                else:
                    filter_string = filter_string+'))' 
                product_count = product_count + 1
            if count < (len(sensordict.keys())-1):        
                filter_string = filter_string+'%20OR%20'
            count = count+1

    return(filter_string)
    
# Load WFS query data

    
def load_hotspots(filter_string, time_period, y_min, x_min, y_max, x_max, max_features, min_confidence):
    
    to_date = dt.datetime.today().strftime('%Y-%m-%d')  
    from_date = (dt.datetime.today() - dt.timedelta(days=time_period)).strftime('%Y-%m-%d')
    #url = f"https://hotspots.dea.ga.gov.au/geoserver/public/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=public:hotspots&outputFormat=application/json&CQL_FILTER=((sensor=%27AVHRR%27%20AND%20(product=%27SRSS%27))%20OR%20(sensor=%27MODIS%27%20AND%20(product=%27MOD14%27))%20OR%20(sensor=%27VIIRS%27%20AND%20(product=%27AFMOD%27%20OR%20product=%27EDR%27)))%20AND%20datetime%20%3E%20%27{from_date}%27%20AND%20datetime%20%3C%20%27{to_date}%27%20AND%20INTERSECTS(location,%20POLYGON(({y_min}%20{x_min},%20{y_min}%20{x_max},%20{y_max}%20{x_max},%20{y_max}%20{x_min},%20{y_min}%20{x_min})))&maxFeatures={max_features}&startIndex=0&sortBy=sensor%20A"
    url = f"https://hotspots.dea.ga.gov.au/geoserver/public/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=public:hotspots&outputFormat=application/json&CQL_FILTER=({filter_string})%20AND%20datetime%20%3E%20%27{from_date}%27%20AND%20datetime%20%3C%20%27{to_date}%27%20AND%20INTERSECTS(location,%20POLYGON(({y_min}%20{x_min},%20{y_min}%20{x_max},%20{y_max}%20{x_max},%20{y_max}%20{x_min},%20{y_min}%20{x_min})))&maxFeatures={max_features}&startIndex=0&sortBy=sensor%20A"
    
    hotspots_gdf = gpd.read_file(url)

    # Filter by confidence
    hotspots_gdf = hotspots_gdf.loc[hotspots_gdf.confidence >= min_confidence]

    # Fix datetime
    hotspots_gdf['datetime'] = pd.to_datetime(hotspots_gdf['start_dt'])

    # Extract required columns
    hotspots_gdf = hotspots_gdf.loc[:, [
            'datetime', 'latitude', 'longitude', 'confidence', 'geometry'
            ]]
    hotspots_gdf.sort_values('datetime', ascending=True, inplace=True)
    logger.info('Hotspots loaded successfully '+str(hotspots_gdf.geometry.total_bounds))
    return(hotspots_gdf)

################################
# Load WMS to xarray #
################################
    
# Create a query object

def wms_xarray(name, url,layer, y_min, x_min, y_max, x_max):

    #TODO - check image exists before recreating
    infile = name+'.tif'
    outfile = name+'_georef.tif'
    
    wms = WebMapService(url, version='1.3.0')
    crs = sorted(wms[layer].crsOptions)
    time = wms[layer].timepositions
    output = wms.getmap(layers=[layer],
                            # TODO Styles true colour = tc doesn't work - figure out why
                            Styles='tc',
                            srs=crs[2],
                            bbox=(x_min, y_min, x_max, y_max),
                            size=(512, 512),
                            format='image/geotiff',
                            # TODO remove specific time reference
                            time=time[5]
                            )

    with open(infile, 'wb') as out:
        out.write(output.read())
    
    epsg ='EPSG:4326'
    logger.info('gdal_translate -a_srs '+epsg+' -a_ullr '+str(x_min)+' '+str(y_max)+' '+str(x_max)+' '+str(y_min)+' '+infile+' '+outfile)
    os.system('gdal_translate -a_srs '+epsg+' -a_ullr '+str(x_min)+' '+str(y_max)+' '+str(x_max)+' '+str(y_min)+' '+infile+' '+outfile)
    logger.info("Background image georeferencing complete") 
    ds = xr.open_rasterio(outfile)
    logger.info("Background image loaded to xarray for plotting")
    return(ds)

#############################
# Generate animation frames #
#############################

# If output folder doesn't exist, create it
def create_outdir(name):
    output_dir = f'frames_{name}'
    os.makedirs(output_dir, exist_ok=True)
    logger.info("Output directory created")
    return(output_dir)
    
    
# Get date/times to iterate through (1 per frame)
def get_dates(hotspots_gdf, frame_freq):
    comp_dates = pd.date_range(hotspots_gdf.datetime.min(), 
                           hotspots_gdf.datetime.max(), 
                           freq=frame_freq)
    logger.info("Dataframe of dates fitting frame frequency created")
    return(comp_dates)
    
    
def run_animation(comp_dates, hotspots_gdf, ds, hotspots_markersize, hotspots_alpha,fade_hours,fade_cmap,hotspots_cmap, x_min, y_min, x_max, y_max,output_dir, timezone):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    for i, comp_date in enumerate(comp_dates):
    
        # Extract only hotspots occuring prior to frame date/time
        hotspots_prev = hotspots_gdf.loc[
            hotspots_gdf['datetime'] < comp_date].copy()
        
        # Calculate hours between frame date and hotspot and sort
        hotspots_prev['hours_before'] = ((
            (comp_date - hotspots_prev['datetime'])).astype('timedelta64[m]') / 60)
        hotspots_prev.sort_values('hours_before', inplace=True, ascending=False)
    
        # Plot Geomedian as underlay
        ds[[0,1,2]].plot.imshow(ax=ax, vmax=500)
        
        # Plot hotspots
        hotspots_prev.plot(ax=ax,
                           column='hours_before',
                           cmap=hotspots_cmap,
                           markersize=hotspots_markersize,
                           alpha=hotspots_alpha,
                           vmin=0,
                           vmax=fade_hours)
        
        # Customise plot and add title
        ax.set_facecolor('black')
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_max, y_min])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('');
        ax.text(0.03, 0.95,
                f"{comp_date.tz_localize(tz='UTC').tz_convert(timezone):%Y-%m-%d}",
                ha='left', 
                va='center', 
                transform=ax.transAxes,
               fontdict={'fontsize': 20, 
                         'color': 'white', 
                         'fontname':'Liberation Sans'})
        
        # Export frame to file
        fig.savefig(f'{output_dir}/hotspots_{i}.jpeg', 
                    bbox_inches='tight',
                    dpi=100,
                    pad_inches=0)
        plt.cla()

    ###########################################
    # Combine into MP4 animation using FFMPEG #
    ###########################################
    #TODD replace with subprocess
    #os.system('ffmpeg -y -r 12 -i $output_dir/hotspots_%d.png -c:v libx264 -vf crop=in_w-15:in_h-15 -pix_fmt yuv420p $output_dir/hotspots_animation.mp4')
    os.system('ffmpeg -y -i $output_dir/hotspots_%d.jpeg -vf crop=in_w-15:in_h-15,minterpolate=fps=24 $output_dir/hotspots_animation.gif')

    
if __name__ == '__main__':
    file_loader = FileSystemLoader("templates")
    env = Environment(loader=file_loader)

    # Get configurations
    satellites = []
    with open(args.configuration, 'r') as config:
        cfg = yaml.load(config, Loader=yaml.Loader)
    
    for configuration in cfg['configurations']:
        logger.info(str(configuration['time_period'])+' day animation for '+configuration['description'])
        
        # Load Hotspots
        hotspots_gdf = load_hotspots(filter(configuration['sensors']),
                                     configuration['time_period'],
                                     configuration['y_min'],
                                     configuration['x_min'],
                                     configuration['y_max'],
                                     configuration['x_max'],
                                     configuration['max_features'], 
                                     configuration['min_confidence'])

        # Load background image as xarray
        ds = wms_xarray(configuration['name'],
                        configuration['url'],
                        configuration['layer'],
                        configuration['y_min'],
                        configuration['x_min'],
                        configuration['y_max'],
                        configuration['x_max']
                        )
        
        ds[[0,1,2]].plot.imshow(vmax=500)
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        comp_dates = get_dates(hotspots_gdf, configuration['frame_freq'])
        output_dir = create_outdir(configuration['name'])
        timezone = configuration['timezone']
        hotspots_markersize = configuration['hotspots_markersize']
        hotspots_alpha = configuration['hotspots_alpha']
        fade_hours = configuration['fade_hours']
        fade_cmap = configuration['fade_cmap']
        
        hotspots_cmap = configuration['hotspots_cmap']
        y_min = configuration['y_min']
        x_min = configuration['x_min']
        y_max = configuration['y_max']
        x_max = configuration['x_max']

        cmap = get_cmap(hotspots_cmap)
        cmap.set_over(fade_cmap)
        
        for i, comp_date in enumerate(comp_dates):
        
            # Extract only hotspots occuring prior to frame date/time
            hotspots_prev = hotspots_gdf.loc[
                hotspots_gdf['datetime'] < comp_date].copy()
            
            # Calculate hours between frame date and hotspot and sort
            hotspots_prev['hours_before'] = ((
                (comp_date - hotspots_prev['datetime'])).astype('timedelta64[m]') / 60)
            hotspots_prev.sort_values('hours_before', inplace=True, ascending=False)
        
            # Plot Geomedian as underlay
            ds[[0,1,2]].plot.imshow(ax=ax, vmax=500)
            
            
            
            # Plot hotspots
            hotspots_prev.plot(ax=ax,
                               column='hours_before',
                               cmap=get_cmap(cmap),
                               markersize=hotspots_markersize,
                               alpha=hotspots_alpha,
                               vmin=0,
                               vmax=fade_hours)
            
            # Customise plot and add title
            ax.set_facecolor('black')
            ax.set_xlim([x_min, x_max])
            ax.set_ylim([y_min, y_max])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('');
            ax.text(0.03, 0.95,
                    f"{comp_date.tz_localize(tz='UTC').tz_convert(timezone):%Y-%m-%d}",
                    ha='left', 
                    va='center', 
                    transform=ax.transAxes,
                   fontdict={'fontsize': 20, 
                             'color': 'white', 
                             'fontname':'Liberation Sans'})
            
            # Export frame to file
            fig.savefig(f'{output_dir}/hotspots_{i}.png', 
                        bbox_inches='tight',
                        dpi=100,
                        pad_inches=0)
            plt.cla()
        os.system('ffmpeg -y -i $output_dir/hotspots_%d.png -vf crop=in_w-15:in_h-15,minterpolate=fps=24 $output_dir/hotspots_animation.gif')
        #exit()
        '''
        run_animation(get_dates(hotspots_gdf, configuration['frame_freq']), hotspots_gdf, ds,
                      configuration['hotspots_markersize'],
                      configuration['hotspots_alpha'],
                      configuration['fade_hours'],
                      configuration['fade_cmap'],
                      configuration['hotspots_cmap'],
                      configuration['y_min'],
                      configuration['x_min'],
                      configuration['y_max'],
                      configuration['x_max'],
                      create_outdir(configuration['name']),
                      configuration['timezone'],
                      )
        '''
    #shape = gpd.read_file('Australian_States_Shapefile/States Map.shp')
    #print(shape.loc[shape['NAME'] == 'Tasmania'].geometry.bounds)
    #TODO - use the ows library to retrieve hotspots
    #wfs11 = WebFeatureService(url='https://hotspots.dea.ga.gov.au/geoserver/public/wfs', version='1.1.0')
    #wfs11.identification.title
    #print(wfs11.get_schema('public:hotspots'))
    #print([operation.name for operation in wfs11.operations])
    #print(list(wfs11.contents))
    #response = wfs11.getfeature(typename='public:hotspots',
    #                            bbox=(132,-26,133,-25),
    #                            srsname='urn:x-ogc:def:crs:EPSG:4326',
    #                            filter='application/json')
    #print(response.read())
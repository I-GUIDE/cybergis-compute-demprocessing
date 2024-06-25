## Import the modules/packages/libraries required
import math
import matplotlib.pyplot as plt
import os
import threading

from pynhd import NLDI
import urllib.request
import rasterio
import rasterio.plot

import geopandas # for reading the shapefile

## Merging/ Mosaicking
from rasterio.merge import merge

## Reprojecting
from rasterio.warp import calculate_default_transform, reproject, Resampling
import pyproj

## Clipping
from rasterio.mask import mask
from shapely.geometry import mapping
from shapely.geometry import Polygon

## Define a function for making a directory depending on whether is exists or not.
## We are creating a function so that it can be used later for creating three folders in the later modules
def check_create_path_func(path):
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print(f"The new directory {path} is created!")
    else:
        print(f"The new directory {path} is not created as it already exists!")

def download_dem_file_func(usgs_filename, local_fileloc_filename):
    try:
        #urllib.request.urlretrieve(url, local_fileloc_filename,MyProgressBar())
        print(f'Beginning file download for {usgs_filename}...')
        urllib.request.urlretrieve(f'https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/{resolution}/TIFF'
                                   f'/current/{usgs_filename}/USGS_{resolution}_{usgs_filename}.tif',
        local_fileloc_filename)
        print(f"Completed file download and saved to {local_fileloc_filename}")
    except Exception as e_value:
        print(f"Error downloading {usgs_filename}.tif: {e_value}")

def merge_dem_raster_func(input_foldername, output_filename):
    ## Get a list of all DEM files in the input folder
    dem_files = [f for f in os.listdir(input_foldername) if f.startswith(f"USGS_{resolution}_")]

    ## Create a list to store the raster datasets
    datasets = []

    ## Open each DEM file and append it to the datasets list
    for dem_file in dem_files:
        file_path = os.path.join(input_foldername, dem_file)
        src = rasterio.open(file_path)
        datasets.append(src)

    ## Merge the raster datasets into a single mosaic
    mosaic, out_trans = merge(datasets)

    ## Copy the metadata from one of the datasets (assuming they all have the same metadata)
    out_meta = datasets[0].meta.copy()
    out_meta.update({
        'height': mosaic.shape[1],
        'width': mosaic.shape[2],
        'transform': out_trans
    })

    # Write the mosaic to the output file
    with rasterio.open(output_filename, 'w', **out_meta) as dest:
        dest.write(mosaic)

    ## Close all the opened datasets
    for dataset in datasets:
        dataset.close()

    print(f"Merging completed for {num_tiles_download} DEM raster files")

def reproject_raster_func(input_raster, output_raster, target_crs):
    # Load the input DEM
    with rasterio.open(input_raster) as src:
        
        # Calculate the transformation and new dimensions
        transform, width, height = calculate_default_transform(src.crs, target_crs, src.width, src.height, *src.bounds)

        # Update metadata for the new dataset
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': target_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        # Create the output dataset and perform the reprojection
        with rasterio.open(output_raster, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest
                )

    print("Reprojection complete.")

def clip_raster_with_shapefile_func(input_raster, output_raster, shapefile):
    ## Open the shapefile using geopandas
    shapefile_gdf = geopandas.read_file(shapefile).buffer(buffer_value_in_meters)
    
    ## Open the input raster
    with rasterio.open(input_raster) as src:
        ## Convert the shapefile geometry to the same CRS as the raster
        shapefile_gdf = shapefile_gdf.to_crs(src.crs)
        
        ## Convert the shapefile geometry to GeoJSON-like format
        geoms = [mapping(geom) for geom in shapefile_gdf.geometry]
        
        ## Clip the raster using the shapefile geometry
        clipped, out_transform = mask(src, geoms, crop=True)
        
        ## Update metadata for the new dataset
        out_meta = src.meta.copy()
        out_meta.update({
            'height': clipped.shape[1],
            'width': clipped.shape[2],
            'transform': out_transform
        })
        
        ## Write the clipped raster to the output file
        with rasterio.open(output_raster, 'w', **out_meta) as dst:
            dst.write(clipped)

    print("Raster clipping complete.")

### main workflow starts here ###

# read input variables from CyberGIS Compute
site_id = os.environ['param_site_id']
resolution=os.environ['param_resolution']

# hardcoded buffer value for now for the clip step
buffer_value_in_meters=150

folder_input=f'/job/executable/input_{site_id}'
folder_results=f'/job/result'
dem_files_store=f'{folder_input}/raw_{site_id}'
check_create_path_func(folder_input)
check_create_path_func(dem_files_store)

watershed=NLDI().get_basins(site_id,fsource='nwissite')
shapefile_fileloc_filename=f'{folder_input}/shape_{site_id}.shp'
watershed.to_file(filename=shapefile_fileloc_filename,
                  driver= 'ESRI Shapefile',
                  mode='w')

extents_basin=watershed.total_bounds
extent_left=abs(math.floor(extents_basin[0]))
extent_right=abs(math.floor(extents_basin[2]))
extent_bottom=abs(math.ceil(extents_basin[1]))
extent_top=abs(math.ceil(extents_basin[3]))

extent = (-extent_left,extent_bottom-1, -extent_right+1, extent_top)
print("Extents of the watershed:",extent)
polygons = []
for lon in range(extent[0],extent[2]):
    for lat in range(extent[1],extent[3]):
        polygon = Polygon([(lon, lat), (lon+1, lat), (lon+1, lat+1), (lon, lat+1)])
        polygons.append({'geometry': polygon, 'lat': lat+1, 'lon': abs(lon)})
        #lat+1 is used as the top left corner is used as the identifier of DEM tile to be downloaded.

gdf=geopandas.GeoDataFrame(polygons,crs="EPSG:4326")

# Perform the intersection
intersection = geopandas.overlay(gdf,watershed, how='intersection')

## Select all columns except geometry
#del attributes_intersect
attributes_intersect = intersection.drop(columns=intersection.geometry.name)

# Function to remove negative sign
overlap_lonlat=list(zip(attributes_intersect['lon'],attributes_intersect['lat']))

### Step 1: Download DEM files based on overlap_lonlat list
usgs_file_list=[]
local_fileloc_filename_list=[]
for location in overlap_lonlat:
    usgs_filename=f'n{location[1]:02d}w{location[0]:03d}'
    usgs_file_list.append(f'n{location[1]:02d}w{location[0]:03d}')
    local_fileloc_filename_list.append(fr'{dem_files_store}/USGS_{resolution}_{usgs_filename}.tif')

## Threading for parallel download to reduce time
threads = []
for usgs_file, filename in zip(usgs_file_list, local_fileloc_filename_list):
    thread = threading.Thread(target=download_dem_file_func,
                              args=(usgs_file, filename))
    thread.start()
    threads.append(thread)

## Wait for all threads to finish to avoid unexpected behavior or incorrect output
for thread in threads:
    thread.join()

## Step 2: Merge the downloaded rasters
num_tiles_download = len(overlap_lonlat)
if num_tiles_download==1:
    merged_raster_filename=f'{dem_files_store}/USGS_{resolution}_{usgs_filename}.tif'
    print("No need to merge as only one tiff file")
else:
    merged_raster_filename=f'{folder_input}/merged_{resolution}_{site_id}.tif'
    merge_dem_raster_func(dem_files_store, merged_raster_filename)

## Step 3: Reproject merged raster to watershed CRS
target_crs = watershed.estimate_utm_crs(datum_name='WGS 84')

## reproject watershed
watershed_file=geopandas.read_file(shapefile_fileloc_filename)
watershed_proj=watershed_file.to_crs(target_crs)
proj_shapefile_filename=f'{folder_results}/proj_shp_{site_id}.shp'
watershed_proj.to_file(proj_shapefile_filename, driver='ESRI Shapefile',mode='w')

## reproject raster
reprojected_raster_filename=f'{folder_input}/reprojected_{resolution}_{site_id}.tif'
reproject_raster_func(merged_raster_filename, reprojected_raster_filename, target_crs)

## Step 4: Clip to watershed extents
clipped_raster_filename=f'{folder_results}/clipped_raster_{site_id}.tif'
clip_raster_with_shapefile_func(reprojected_raster_filename,clipped_raster_filename,proj_shapefile_filename)

## Step 5: create a preview
with rasterio.open(clipped_raster_filename) as src:
    fig, ax = plt.subplots(figsize = (8,8))
    mapped=rasterio.plot.show(src,
                              ax=ax, 
                              cmap='viridis')
    im = mapped.get_images()[0]
    fig.colorbar(im,
                 ax=ax,
                 shrink=0.5,
                 label='Elevation (in meters)')

## Uee the reprojected watershed shapefile
watershed_proj.plot(ax=ax, 
           facecolor='none', 
           edgecolor='red')
src.close()

ax.set_xlabel("Easting (meters)")
ax.set_ylabel("Northing (meters)")
ax.set_title("Clipped Raster DEM "+site_id)
ax.set_title(f"Clipped Raster DEM with a buffer of {buffer_value_in_meters} m and with {target_crs} Projected CRS\n(for USGS {site_id})")
check_create_path_func(f'{folder_results}/preview') 
clip_fig_filepath = f'{folder_results}/preview/clipped_raster_{site_id}.png'
fig.savefig(clip_fig_filepath)
plt.close(fig)

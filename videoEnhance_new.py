import cv2
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
import geopandas as gpd
import rasterio as rio
from rasterio.plot import plotting_extent
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
import shapefile
from os.path import isfile, join
from pathlib import Path


# Set figure size and title size of plots
plt.rcParams['figure.figsize'] = (14, 14)
plt.rcParams['axes.titlesize'] = 20


def get_im_pairs(mission_csv, enhanceFolder, saveFolder):
    df = pd.read_csv(mission_csv)
    waypoints = df.waypoint.unique()
    new_waypoints = ['WP7']+list(waypoints[0:-1])
    waypoints = new_waypoints
    
    i=0
    pairs = []
    for wp in waypoints:
        df_filt = df[df['waypoint']==wp]
        df_filt = df_filt.reset_index()
        for i in range(len(df_filt)):
            im_path = df_filt['local_path'][i]
            new_im_path = os.path.splitext(os.path.basename(im_path))[0]
            if wp == 'WP7':
                wp = 'WP07'
            enhance = enhanceFolder
            for image in glob.glob(enhance + '/*.jpg'):
                enhanceName = os.path.basename(image)
                if enhanceName[0:8]==new_im_path:
                    enhancePair = (im_path, image)
                    pairs.append(enhancePair)
    return pairs


# =============================================================================
# Working but file does not have a spatial coordinate system, have to project in Arc
# =============================================================================
def pyshp_geobox_to_shapefiles(inFile, outFile):
    # funtion to generate a .prj file
    box = shapefile.Writer(outFile, shapeType=shapefile.POLYGON)
    box.field('label', 'C', size=20)
    df = pd.read_csv(inFile)
    for i in range(len(df)):
        llx = df['llx'][i]
        lly = df['lly'][i]
        lrx = df['lrx'][i]
        lry = df['lry'][i]
        urx = df['urx'][i]
        ury = df['ury'][i]
        ulx = df['ulx'][i]
        uly = df['uly'][i]
        image = df['image'][i]
        box.poly([[[llx,lly],[lrx,lry],[urx,ury],[ulx,uly]]])
        box.record(image)
    # create the PRJ file
    prj = open("%s.prj" % outFile, "w")
    epsg = 'GEOGCS["WGS 84",'
    epsg += 'DATUM["WGS_1984",'
    epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
    epsg += ',PRIMEM["Greenwich",0],'
    epsg += 'UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()
    
def plotMap(image_path, shapefile_path, raster_path, ax):
    image = os.path.basename(image_path)

    # Open fire boundary data with geopandas
    shape_boundary = gpd.read_file(shapefile_path)
    shape_boundary = shape_boundary[shape_boundary['label']==image].reset_index()
    

    # Open NAIP data in read ('r') mode
    with rio.open(raster_path) as raster_src:
        raster_data = raster_src.read()

        # Project fire boundary to match NAIP data
        shape_bound = shape_boundary.to_crs(raster_src.crs)

        # Create plotting extent from DatasetReader object
        raster_plot_extent = plotting_extent(raster_src)


    ep.plot_rgb(raster_data,
                rgb=[0, 1, 2],
                ax=ax,
                extent=raster_plot_extent)  # Use plotting extent from DatasetReader object

    shape_bound.plot(ax=ax, color='k',alpha=1)
    plt.title('Image Footprint')
    plt.tight_layout()

def makePlots(pairs, saveFolder, shapefile_path, raster_path):
    i=0
    for pair in pairs:
        print(i/len(pairs))
        if i<441:
            i=i+1
            continue
        # Plot uncropped array
        plt.suptitle(os.path.splitext(os.path.basename(pair[0]))[0])
        plt.subplot(2,2,1)
        plt.title('Original')
        plt.imshow(cv2.imread(pair[0]))
        plt.xticks([],[])
        plt.yticks([],[])
        plt.subplot(2,2,2)
        plt.imshow(cv2.imread(pair[1]))
        plt.xticks([],[])
        plt.yticks([],[])
        plt.title('Enhanced')
        plt.subplot(2,1,2)
        ax = plt.gca()
        plotMap(pair[0],
                shapefile_path,
                raster_path, ax)
        #plt.show()
        plt.savefig(os.path.join(saveFolder,str(i) + '.png'), dpi=300)
        plt.close()
        i=i+1
def makeVideo(frameFolder, vidPath, rate):
    pathIn= frameFolder
    pathOut = vidPath
    fps = rate
    files = sorted(Path(frameFolder).iterdir(), key=os.path.getmtime)
    im = cv2.imread(os.path.join(frameFolder, files[0]))
    height,width,layers = im.shape
    size = (width,height)
    out = cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
    for i in range(len(files)):
        print(i/len(files))
        filename=os.path.join(pathIn,files[i])
        #reading each frame
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        out.write(img)
        #inserting the frames into an image array
        img=None
    out.release()

#pyshp_geobox_to_shapefiles(r'C:\MarkLundineSurface\iver\Camera\uss_nina_footprint_seafloor_enhance.csv', r'C:\MarkLundineSurface\iver\Camera\uss_nina_footprint_seafloor_enhance.shp')
##plotMap(r'C:\MarkLundineSurface\iver\Camera\WP27\VC1\WP27_130.jpg',
##        r'C:\MarkLundineSurface\iver\Camera\uss_nina_footprint_seafloor_enhance.shp',
##        r'C:\MarkLundineSurface\iver\nina_sidescan.tif')
#pairs = get_im_pairs(r'C:\MarkLundineSurface\iver\Camera\uss_nina_footprint_seafloor_enhance.csv', r'C:\MarkLundineSurface\iver\compare_enhance')
#makePlots(pairs)
#makeVideo(r'C:\MarkLundineSurface\iver\compare_enhance', r'C:\MarkLundineSurface\iver\nina_track.mp4',10)


        

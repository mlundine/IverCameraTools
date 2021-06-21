"""
author: Mark Lundine, CSHEL, Univ. of Delaware
"""


import numpy as np
import simplekml 
import cv2
import os
import glob
import pandas as pd
from PIL import Image, ExifTags
import geopy
from geopy.distance import geodesic

def get_missions(camera_path):
    """
    Gets all the subdirectories for the
    images taken on the Iver Mission
    Inputs: camera_path (str), full filepath to the Camera data taken by the Iver
    Outputs: subfolders (list), list of all the waypoint folders in the camera folder
    """
    subfolders = [ f.path for f in os.scandir(camera_path) if f.is_dir() ]
    for i in range(len(subfolders)):
        elem = subfolders[i]
        new_elem = os.path.join(elem, 'VC1')
        subfolders[i] = new_elem
    return subfolders


def get_footprint(center_lat, center_long,
                  altitude, pitch, roll, heading,
                  sensor_x=10.67, sensor_y=8.00, focal_length=12):
    """
    Calculates footprint of Iver image
    inputs:
    center_lat is latitude in degrees
    center_long is longitude in degrees
    altitude is altitude in m
    pitch is pitch in degrees
    roll is roll in degrees
    heading is heading in degrees
    
    optional inputs:
    sensor_x is sensor width in mm 10.67 for Iver which has a Type 1/1.2 sensor
    sensor_y is sensor height in mm 8.00 for Iver which has a Type 1/1.2 sensor
    focal_length is focal length of lens in mm

    outputs:
    footprint_dict containing the following keys
    ulx:  upper left lon coordinate of image in degrees
    uly:  upper left lat coodrinate of image in degrees
    urx:  upper right lon coordinate of image in degrees
    ury:  upper right lat coordinate of image in degrees
    llx:  lower left lon coordinate of image in degrees
    lly:  lower left lat coordinate of image in degrees
    lrx:  lower right lon coordinate of image in degrees
    lry:  lower right lat coodrinate of image in degrees
    geo_width:  width of image in meters
    geo_height: height of image in meters
    """

    ##get footprint in meters based on lens focal length,
    ##sensor dimensions,
    ##pitch, roll, and altitude of Iver
    fov_wide = 2*np.arctan(sensor_x/(2*focal_length)) #angle radians
    fov_tall = 2*np.arctan(sensor_y/(2*focal_length)) #angle radians
    pitch = pitch*np.pi/180
    roll = roll*np.pi/180
    
    sensor_to_bottom = altitude*np.tan(pitch-0.5*fov_tall) #m
    sensor_to_top = altitude*np.tan(pitch+0.5*fov_tall) #m
    sensor_to_left = altitude*np.tan(roll-0.5*fov_wide) #m
    sensor_to_right = altitude*np.tan(roll+0.5*fov_wide) #m
    
    footprint_height = sensor_to_top-sensor_to_bottom #m
    footprint_width = sensor_to_right-sensor_to_left #m
    footprint_diagonal = np.sqrt((footprint_height)**2+(footprint_width)**2)
    center_to_corner = footprint_diagonal/2

    ##use image GPS data as center point
    ##use Iver heading
    ##find corner coordinates

    ##define center and various bearings
    origin = geopy.Point(center_lat, center_long)
    ##center to upperleft, center to upperright, center to lowerright, center to lowerleft
    bearings = [heading-45, heading+45, heading+135, heading+225]

    ##calculate terminal points, center to each corner of image
    ##upperleft
    destination = geodesic(kilometers=center_to_corner/1000).destination(origin, bearings[0])
    uly, ulx = destination.latitude, destination.longitude
    ##upperright
    destination = geodesic(kilometers=center_to_corner/1000).destination(origin, bearings[1])
    ury, urx = destination.latitude, destination.longitude    
    ##lowerright
    destination = geodesic(kilometers=center_to_corner/1000).destination(origin, bearings[2])
    lry, lrx = destination.latitude, destination.longitude
    ##lowerleft
    destination = geodesic(kilometers=center_to_corner/1000).destination(origin, bearings[3])
    lly, llx = destination.latitude, destination.longitude

    footprint_dict = {
                      'uly':uly,
                      'ulx':ulx,
                      'ury':ury,
                      'urx':urx,
                      'lry':lry,
                      'lrx':lrx,
                      'lly':lly,
                      'llx':llx,
                      'geo_width':footprint_width,
                      'geo_height':footprint_height,
                      }
    return footprint_dict

def get_meta_data_indv(image_path):
    """
    reads Iver image metadata and outputs as a dictionary
    input:
    image_path (str), filepath to image
    output:
    meta_data_dict (dictionary), contains all of the image metadata
    """

    ### heading, depth, altitude, roll, pitch and speed
    img = Image.open(image_path)
    width,height = img.size
    exifData = {}    
    exifDataRaw = img._getexif()
    for tag, value in exifDataRaw.items():
        decodedTag = ExifTags.TAGS.get(tag, tag)
        exifData[decodedTag] = value
    time = exifData['DateTimeOriginal']
    userComment = exifData['UserComment'].decode('utf-8')

    lat_index = userComment.find('La')+3
    substr = userComment[lat_index:]
    lat = substr[0:substr.find(',')]

    long_index = userComment.find('Ln')+3
    substr = userComment[long_index:]
    long = substr[0:substr.index(',')]
    
    heading_index = userComment.find('H')+2
    substr = userComment[heading_index:]
    heading = substr[0:substr.index(',')]
    
    depth_index = userComment.find('D')+2
    substr = userComment[depth_index:]
    depth = substr[0:substr.index(',')]
    
    altitude_index = userComment.find('A')+2
    substr = userComment[altitude_index:]
    altitude = substr[0:substr.index(',')]
    
    roll_index = userComment.find('R')+2
    substr = userComment[roll_index:]
    roll = substr[0:substr.index(',')]
    
    pitch_index = userComment.find('P')+2
    substr = userComment[pitch_index:]
    pitch = substr[0:substr.index(',')]
    
    speed_index = userComment.find('S')+2
    substr = userComment[speed_index:]
    speed = substr[0:substr.index(',')]

    meta_data_dict = {'im':os.path.basename(image_path),
                      'path':image_path,
                      'latitude':lat,
                      'longitude':long,
                      'heading':heading,
                      'depth':depth,
                      'altitude':altitude,
                      'roll':roll,
                      'pitch':pitch,
                      'speed':speed,
                      'time':time,
                      'width':width,
                      'height':height}
    return meta_data_dict

def get_meta_data_mission(mission_id,camera_path):
    """
    saves metadata for each image to a csv
    inputs:
    mission_id (str), just a name for the mission
    camera_path (str), full filepath to the Camera directory of the Iver mission
    output:
    save_path, the filepath to the metadata csv this function saves
    """
    waypoints = get_missions(camera_path)
    mission_meta_data = [['image','local_path','latitude','longitude',
                          'heading','depth','altitude','roll','pitch',
                          'speed', 'time','width','height','waypoint',
                          'uly','ulx',
                          'ury','urx',
                          'lry','lrx',
                          'lly','llx',
                          'geo_width','geo_height']]
    for wp in waypoints:
        for im in glob.glob(wp + '/*.jpg'):
            meta = get_meta_data_indv(im)
            element = list(meta.values())
            element.append(os.path.basename(os.path.dirname(wp)))
            footprint = get_footprint(float(meta['latitude']),-float(meta['longitude']),
                                      float(meta['altitude']),float(meta['pitch']),
                                      float(meta['roll']),float(meta['heading']))
            element = element + list(footprint.values())
            mission_meta_data.append(element)
    save_path = os.path.join(camera_path,mission_id+'.csv')
    np.savetxt(save_path, mission_meta_data, delimiter=",", fmt='%s')
    return save_path

def kmz(mission_csv, photo_or_ground):
    """
    makes kmzs of Iver images
    inputs:
    mission_csv, the filepath to the csv with all of the image metadata
    photo_or_ground, 'photo' for photo overlay kmz, 'ground' for ground overlay kmz
    'ground' makes the images displayed as their actual footprint size
    'photo' just places them at the lat, lon coordinate, not accurate footprint
    """
    df = pd.read_csv(mission_csv)
    waypoints = df.waypoint.unique()
    for wp in waypoints:
        df_filt = df[df['waypoint']==wp]
        df_filt = df_filt.reset_index()
        kml = simplekml.Kml()
        for i in range(len(df_filt)):
            im_path = df_filt['local_path'][i]
            new_im_path = os.path.splitext(os.path.basename(im_path))[0]
            enhance = r'C:\imagecor\Image_Correction_App\20210609_USS_Nina_Iver_Photos_Enhanced'
            if wp == 'WP7':
                wp = 'WP07'
            enhance = os.path.join(enhance,wp)
            for image in glob.glob(enhance + '/*.jpg'):
                enhanceName = os.path.basename(image)
                if enhanceName[0:8]==new_im_path:
                    im_path = df_filt['local_path'][i]
                    kml.addfile(image)
                    im = df_filt['image'][i]
                    lat = df_filt['latitude'][i]
                    long = df_filt['longitude'][i]
                    alt = df_filt['altitude'][i]
                    head = df_filt['heading'][i]
                    llx = df_filt['llx'][i]
                    lly = df_filt['lly'][i]
                    ll = (llx,lly)
                    lrx = df_filt['lrx'][i]
                    lry = df_filt['lry'][i]
                    lr = (lrx,lry)
                    urx = df_filt['urx'][i]
                    ury = df_filt['ury'][i]
                    ur = (urx,ury)
                    ulx = df_filt['ulx'][i]
                    uly = df_filt['uly'][i]
                    ul = (ulx,uly)
                    if photo_or_ground == 'ground':
                        ground = kml.newgroundoverlay(name=im,
                                                      altitude=alt,
                                                      altitudemode = simplekml.GxAltitudeMode.relativetoseafloor)
                        ground.icon.href = image
                        ground.gxlatlonquad.coords = [ll,lr,
                                                      ur,ul]
                    if photo_or_ground =='photo':
                        photo = kml.newphotooverlay(name=im)
                        photo.camera = simplekml.Camera(longitude=-long, latitude=lat, altitude=alt,heading=head,
                                                        altitudemode=simplekml.AltitudeMode.clamptoground)
                        photo.point.coords = [(-long,lat)]
                        photo.style.iconstyle.icon.href = im_path
                        photo.icon.href = im_path
                        photo.viewvolume = simplekml.ViewVolume(-25,25,-15,15,1)
        kmlpath = os.path.splitext(mission_csv)[0]+wp+'.kmz'
        kml.savekmz(kmlpath)
        kml = None
    return kmlpath
    
def main(mission_id,camera_folder, photo_or_ground):
    """
    Will run all functions, reads the metadata of Iver mission photos
    Saves metadata to csv
    Makes kmzs for each line
    csvs and kmzs get saved to the Camera folder you specify
    inputs:
    mission_id (str): feed a string that defines the name of the mission, whatever you want
    camera_folder (str): full filepath to the camera folder of the iver mission
    photo_or_ground (str): 'ground' for ground overlay kmz, 'photo' for photo overlays
    
    """
    mission_csv = get_meta_data_mission(mission_id,camera_folder)
    print('Metadata read and saved to csv')
    kmz_path = kmz(mission_csv,photo_or_ground)
    print('Kmzs made')
main('uss_nina_footprint_seafloor_enhance',r'C:\MarkLundineSurface\iver\Camera','ground')

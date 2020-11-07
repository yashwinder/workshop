import os
import gpxpy
import gpxpy.gpx
import numpy as np
import pandas as pd
from geopy import distance
import matplotlib.pyplot as plt
from datetime import datetime

data = []

def main(path):

    global data

    files = os.listdir(path)
    
    for file in files:
        gpx_file = open(os.path.join(path, file), 'r')
        gpx = gpxpy.parse(gpx_file)

        segment = gpx.tracks[0].segments[0]
        coords = pd.DataFrame([{
            'lat': p.latitude,
            'long': p.longitude,
            'ele': p.elevation,
            'time': p.time} for p in segment.points])
        
        pair_of_coords = {}
        for i, p in enumerate(segment.points):
            pair_of_coords[(p.latitude, p.longitude)] = i
        
        data.append((coords, pair_of_coords))

    return data

def get_distance_elevation(route):

    lat, long, ele = route

    distance_covered = 0
    elevation_gain = 0

    for i in range(len(lat) - 1):
        distance_covered += distance.geodesic((lat[i], long[i]), (lat[i + 1], long[i + 1])).km
        elevation_gain += max(0, ele[i + 1] - ele[i])
    
    return (distance_covered, elevation_gain)
def get_all_stats(routes):


    distance_covered, elevation_gain = get_distance_elevation(routes[0])
    
    indexes = [i for i in range(len(routes))]
    speeds = [distance_covered / (len(lat) / 3600) for lat,_,_ in routes]

    ret_info = {
        'distance_covered': distance_covered,
        'elevation_gain': elevation_gain,
        'speed_plot': (indexes, speeds)
    }

    return ret_info

def check_uniqueness(routes):

    dist_arr = np.array([get_distance_elevation(route)[0] for route in routes])
    ele_arr = np.array([get_distance_elevation(route)[1] for route in routes])

    dist_mean = sum(dist_arr) / len(dist_arr)
    ele_mean = sum(ele_arr) / len(ele_arr)

    dist_arr = dist_arr - [dist_mean]
    ele_arr = ele_arr - [ele_mean]

    for i in range(len(dist_arr)):
        if (abs(dist_arr[i]) > 0.1 or abs(ele_arr[i]) > 1):
            return False
    
    return True

def get_coordinates_info(start, end, mid = (0, 0)):

    global data

    routes = []
    for i in range(len(data)):
        coords, pair_of_coordinates = data[i]
        if start in pair_of_coordinates and end in pair_of_coordinates and (mid == (0, 0) or mid in pair_of_coordinates):
            
            idx_start = pair_of_coordinates[start]
            idx_mid = -1 if mid == (0, 0) else pair_of_coordinates[mid]
            idx_end = pair_of_coordinates[end]
            
            if (idx_start < idx_end and (idx_mid == -1 or (idx_start < idx_mid and idx_mid < idx_end))):
                routes.append((coords['lat'].tolist()[idx_start : idx_end], coords['long'].tolist()[idx_start : idx_end], coords['ele'].tolist()[idx_start : idx_end]))
        else:
            continue
    
    if not(check_uniqueness(routes)):
        return False

    return get_all_stats(routes)


def get_attr_per_day():
    '''
    
    Returns a tuple of dictionaries for distance vs day, elevation-gain vs day and speed vs day plots
    
    ''' 

    dist_map, ele_map, speed_map = dict(), dict(), dict()
    for i in range(len(data)):
        time_series, ele_series, lat_series, long_series = data[i][0]['time'], data[i][0]['ele'], data[i][0]['lat'], data[i][0]['long'] 
        for j in range(len(data[i][0])):
            day = time_series[j].strftime("%x")
            if day not in dist_map:
                dist_map[day] = []
            if day not in ele_map:
                ele_map[day] = []
            ele_map[day].append(ele_series[j])
            dist_map[day].append((lat_series[j], long_series[j], time_series[j]))

    for day, arr in ele_map.items():
        ele_gain = 0
        for i in range(len(arr) - 1):
            ele_gain += max(0, arr[i + 1] - arr[i])
        ele_map[day] = ele_gain
    
    FMT = "%H:%M:%S"
    
    for day, arr in dist_map.items():
        total_dist = total_time = 0
        for i in range(len(arr) - 1):
            tmp_dist = distance.geodesic((arr[i][0], arr[i][1]), (arr[i + 1][0], arr[i + 1][1])).km
            time_elapsed = datetime.strptime(arr[i + 1][2].strftime("%X"),FMT) - datetime.strptime(arr[i][2].strftime("%X"),FMT)
            days, seconds = time_elapsed.days, time_elapsed.seconds
            hours = days * 24 + seconds // 3600
            minutes = (seconds % 3600) // 60
            seconds = seconds % 60
            total_dist += tmp_dist
            total_time += (hours + minutes/60 + seconds/3600)
        dist_map[day] = total_dist
        speed_map[day] = total_dist/total_time
    
    return (dist_map, ele_map, speed_map)

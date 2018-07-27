#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:42:06 2018

@author: kcooper
"""
import pandas as pd
import numpy as np

data = pd.read_csv("index_obsidNames.csv")
uniquenames = pd.read_csv("HST_obsids.csv")
data = pd.DataFrame(data)


#data = np.genfromtxt("index_obsidNames.csv", delimiter=',', dtype='|U20')
#uniquenames = np.genfromtxt("HST_obsids.csv", delimiter=',', dtype='|U20')[1:,1:]
realnames = list(data.columns.values)[1:]
#print(realnames)
datanew = data.values[:,1:]
uniquenames = uniquenames.values[:,1:]
uniquenames = pd.DataFrame(uniquenames)

#realnames = data[0, 1:]
#datanew = data[1:, 1:]

y_shape, x_shape = datanew.shape

t = [(u, []) for u in uniquenames[0]]
dic = dict(t)

for u in uniquenames[0]:
    for i in range(x_shape):
        for j in range(y_shape):
            if datanew[j, i] == u:
                dic[u].append(realnames[i])
            else:
                pass


from astroquery.mast import Observations
from astropy import units as u
from astropy.coordinates import SkyCoord
import operator

def find_sources():
    chandra = pd.read_csv("chandraSources.csv")
    sourceNames = pd.DataFrame(chandra)
    obsTable = Observations.query_criteria(obs_id=uniquenames[0], calib_level=3, dataRights="public", obs_collection="HST", instrument_name=["ACS/WFC", "WFC3/UVIS"], em_min=[3e-07, 9.7e-07], em_max=[3e-07, 9.7e-07], filters=["%W"])
    obs = list(obsTable['obs_id'])
    overlap = pd.DataFrame(columns=('HST_obsid', 'sourceRA', 'sourceDec', 'objectName'))
    for obsid in uniquenames[0]:
        ### index number (number of iterations) = m
#        m = uniquenames.loc[uniquenames[0] == obsid].index[0]
        m = obs.index(obsid)
        x = []
        y = []
        name = []
        for n in range(len(dic[obsid])):
            chandraSources = sourceNames.loc[sourceNames['resolvedObject'] == dic[obsid][n]]
            chandraSources = chandraSources.reset_index(drop=True)
            x_n = list(chandraSources['sourceRA'])
            x.append(x_n)
            y_n = list(chandraSources['sourceDec'])
            y.append(y_n)
            name_n = list(chandraSources['resolvedObject'])
            name.append(name_n)
        x = reduce(operator.concat, x)
        y = reduce(operator.concat, y)
        name = reduce(operator.concat, name)
        mark = overlap_index(obsTable, x, y, name, m, obsid)
        df = pd.DataFrame(mark, columns=('HST_obsid', 'sourceRA', 'sourceDec', 'objectName'))
        overlap = overlap.append(df, ignore_index=True)
    overlap.to_csv('sourceOverlapIndex.csv')
# test for if a string is a float
def is_float(string):
    """Returns True if string is a float (can be converted to a float without an error)
    """
    try:
        float(string)
        return True
    except ValueError:
        return False

          
def overlap_index(table, x, y, name, m, obsid):
    mark = []
    poly = table['s_region'][m]
    poly = poly.split()
    poly.pop(0)
    polygon = poly
    # if the item of the list after 'POLYGON' is not a float, send the coordinates through point_in_poly
    if is_float(polygon[0]) == False:
            # pop all items at the beginning of the list that are not floats
        while is_float(polygon[0]) == False:
            polygon.pop(0)
        polygon = [float(i) for i in polygon]
        for i in range(len(x)):
            x_i = x[i]
            y_i = y[i]
            n_i = name[i]
            point_in_poly(x_i, y_i, polygon)
            if point_in_poly(x_i, y_i, polygon) == True:
                mark.append([obsid, x_i, y_i, n_i])    

        # if 'POLYGON' is the only string in the list, convert the coordinates to ICRS andsend the coordinates through point_in_poly_coord
    elif is_float(polygon[0]) == True:
      #      polygon = poly
      polygon = [float(i) for i in polygon]
      poly_ra, poly_dec = poly_icrs(polygon)
      for i in range(len(x)):
          x_i = x[i]
          y_i = y[i]
          n_i = name[i]
          point_in_poly_coord(x_i, y_i, poly_ra, poly_dec)
          if point_in_poly_coord(x_i, y_i, poly_ra, poly_dec) == True:
              mark.append([obsid, x_i, y_i, n_i])
    return mark
    
            
#==============================================================================
# def create_poly(table, x, y):
#     """ marks all chandra sources in each polygon that have overlaying ra and dec """
#     m=0
#     mark = []
#     for m in range(len(table)):
#         poly = table['s_region'][m]
#         poly = poly.split()
#         poly.pop(0)
#         polygon = poly
#         # if the item of the list after 'POLYGON' is not a float, send the coordinates through point_in_poly
#         if is_float(polygon[0]) == False:
#             # pop all items at the beginning of the list that are not floats
#             while is_float(polygon[0]) == False:
#                 polygon.pop(0)
#             polygon = [float(i) for i in polygon]
#             for i in range(len(x)):
#                 x_i = x[i]
#                 y_i = y[i]
#                 point_in_poly(x_i, y_i, polygon)
#                 if point_in_poly(x_i, y_i, polygon) == True:
#                     mark.append(m)
#                     break
#         # if 'POLYGON' is the only string in the list, convert the coordinates to ICRS andsend the coordinates through point_in_poly_coord
#         elif is_float(polygon[0]) == True:
#       #      polygon = poly
#             polygon = [float(i) for i in polygon]
#             poly_ra, poly_dec = poly_icrs(polygon)
#             for i in range(len(x)):
#                 x_i = x[i]
#                 y_i = y[i]
#                 point_in_poly_coord(x_i, y_i, poly_ra, poly_dec)
#                 if point_in_poly_coord(x_i, y_i, poly_ra, poly_dec) == True:
#                     mark.append(m)
#                     break
#     return mark
#==============================================================================

def poly_icrs(poly):
    """Converts polygon coordinates from fk5 to icrs to match chandra sources """
    n = len(poly)
    poly_ra = []
    poly_dec = []
    for i in range(0,n, 2):
        a = poly[i]
        b = poly[i+1]
        c = SkyCoord(ra=a*u.degree, dec=b*u.degree, frame='fk5')
        c_icrs = c.icrs
        ra = c_icrs.ra.degree
        dec = c_icrs.dec.degree
        poly_ra.append(ra)
        poly_dec.append(dec)
    return poly_ra, poly_dec
#poly_ra, poly_dec = poly_icrs(polygon)

def point_in_poly_coord(x, y, poly_ra, poly_dec):
    """Returns True if (x,y) lies within poly and False otherwise.
    Where poly:= [ (x1,y1), (x2,y2), (x3,y3,...(xn,yn)] given by output of poly_icrs
    with ra and dec in seperate lists
    ray casting algorithim"""
    
    n = len(poly_ra)
    inside = False
    p1x = poly_ra[0]
    p1y = poly_dec[0]
    ##for i in range(2, n+1, 2)
    for i in range(n+1):
        p2x = poly_ra[i % n]
        p2y = poly_dec[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x<= max(p1x, p2x):
                    if p1y != p2y:
                        xints = ((y-p1y)*(p2x-p1x))/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
        
    return inside

def point_in_poly(x,y,poly):
    """Returns True if (x,y) lies within poly and False otherwise.
    Where poly:= [ (x1,y1), (x2,y2), (x3,y3,...(xn,yn)]
    ray casting algorithim"""
    
    n = len(poly)
    inside = False
    p1x = poly[0]
    p1y = poly[1]
    for i in range(2, n-1, 2):
        p2x = poly[i % n]
        p2y = poly[(i % n)+1]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x<= max(p1x, p2x):
                    if p1y != p2y:
                        xints = ((y-p1y)*(p2x-p1x))/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
        
    return inside

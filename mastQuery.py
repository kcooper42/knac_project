#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 11:16:19 2018

@author: kcooper
"""

from astroquery.mast import Observations
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd
##print(obsTable)
resolvedName = "MESSIER 083"
obsTable = Observations.query_criteria(objectname=resolvedName, calib_level=3, dataRights="public", obs_collection="HST", instrument_name=["ACS/WFC", "WFC3/UVIS"], em_min=[3e-07, 9.7e-07], em_max=[3e-07, 9.7e-07], filters=["%W"])
#for i in range(len(obsTable)):
 #   poly = obsTable['s_region'][n]
#poly = obsTable['s_region'][11]
#polygon = poly[13:117]
chandraSources = pd.read_csv("M83sources.csv")
x = chandraSources['sourceRA']
y = chandraSources['sourceDec']
def create_poly():
    m=0
    mark = []
    for m in range(len(obsTable)):
        poly = obsTable['s_region'][m]
        poly = poly.split()
        poly.pop(0)
        if poly[0] == "ICRS":
            polygon = poly
            polygon.pop(0)
            polygon = [float(i) for i in polygon]
            for i in range(len(x)):
                x_i = x[i]
                y_i = y[i]
                point_in_poly(x_i, y_i, polygon)
                if point_in_poly(x_i, y_i, polygon) == True:
                    mark.append(m)
                    break
        elif poly[0] == "OTHER":
            polygon = poly
            polygon.pop(0)
            polygon = [float(i) for i in polygon]
            for i in range(len(x)):
                x_i = x[i]
                y_i = y[i]
                point_in_poly(x_i, y_i, polygon)
                if point_in_poly(x_i, y_i, polygon) == True:
                    mark.append(m)
                    break
        else:
            polygon = poly
            polygon = [float(i) for i in polygon]
            poly_ra, poly_dec = poly_icrs(polygon)
            for i in range(len(x)):
                x_i = x[i]
                y_i = y[i]
                point_in_poly_coord(x_i, y_i, poly_ra, poly_dec)
                if point_in_poly_coord(x_i, y_i, poly_ra, poly_dec) == True:
                    mark.append(m)
                    break
    return mark
#poly = obsTable['s_region'][0]
#polygon = poly[8::]
#polygon = polygon.replace(" ", ",")
#polygon = np.fromstring(polygon, dtype=np.float, sep=',')
def poly_icrs(poly):
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
def point_in_poly_coord(x,y,poly_ra, poly_dec):
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

def obs_id():
    mark = create_poly()
    obs_id =[]
    for i in mark:
        obs_id.append(obsTable['obs_id'][i])
    return obs_id
def filter_products():
    obsids = obs_id()
    obs = obsTable['obsid']
    dataProductsByID = Observations.get_product_list(obs)
    dataProductsByID = Observations.filter_products(dataProductsByID, obs_id=obsids, productSubGroupDescription="DRZ")
    print(len(dataProductsByID))
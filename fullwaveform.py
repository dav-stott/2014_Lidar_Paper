# -*- coding: utf-8 -*-
"""
Created on Mon May 20 20:08:46 2013
Stuff to pull out stuff for full waveform

@author: davstott
davstott@gmail.com
"""

import os
import numpy as np
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt

import fiona
from shapely.geometry import shape
from shapely.geometry import asPoint

import csv

#******************* PATHS****************************
dirpath = os.path.dirname(os.path.abspath('...'))
print ' running in: ', dirpath
datapath = os.path.join(dirpath,'data')
print 'data located in :', datapath

outputpath = os.path.join(dirpath,'output')
if not os.path.exists(outputpath):
    os.mkdir(outputpath)
    
#poly_path = os.path.join(dirpath,'polys')


header = str('x'+','+
             'y'+','+
             'z'+','+
             'intensity'+','+
             'peak_start'+','+
             'peak_end'+','+
             'peak_location'+','+
             'peak_width'+','+
             'max_intensity'+','+
             'peak_sum'+','+
             'shoulder_location'+','+
             'shoulder_intensity'+
             '\n')



class ShpReader():
    def __init__(self,indir):
        
        print 'call to shape reader'
        os.chdir(indir)
        listindir = os.listdir(indir)
        
        self.arc = []
        self.bac = []
        self.poly = None
        
        for thing in listindir:
            item = thing[-4:]
            if item == '.shp':
                shapefile = fiona.open(thing)
                for pol in shapefile:
                    #print pol['geometry']
                    classification = pol['properties']['CLASS']
                    if classification == 'ARC':
                        #print 'ARC'
                        arc_poly = shape(pol['geometry'])
                        self.arc.append(arc_poly)
                      
                    elif classification =='BAC':
                        bac_poly = shape(pol['geometry'])
                        self.bac.append(bac_poly)
                        
                        
                shapefile.close()
        print len(self.arc),len(self.bac)     
        
        
    def classify_point(self,x,y):
        #print 'call to classify point'
        
        reading = asPoint(np.array([x,y]))
        #print reading.wkt
        
        classification = None
        
        
        for poly in self.arc:
            if reading.within(poly):
                classification = 1
                print 'arc'
                break
                
        for poly in self.bac:
            if reading.within(poly):
                classification = 0
                print 'bac'
                break
               
        #print classification
        return classification




#*******************functions*****************************

#a smoothing spline
def smoothing(waveform, kparam, sparam, weights):
    sm_x = np.arange(1,257,1)
    sm_weights = np.zeros(256)+weights
    sm_spline = interpolate.UnivariateSpline(sm_x, 
                                             waveform, 
                                             k=kparam,
                                             w=sm_weights,
                                             s=sparam)

    spline = sm_spline(sm_x)
    return spline  
    

#***************** Parameters*********************************
#for spline
kparam = 1.3
sparam = 191
weights = 1.4

#x values
x_vals = np.arange(1,257,1)

#find the data
os.chdir(datapath)
indir = os.listdir(datapath) 
#open the data in a for loop

#classes = ShpReader(poly_path)

os.chdir(datapath)
indir = os.listdir(datapath) 


for file in indir:
    print file
    reading_count = 0
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        j = 0

        for row in reader:
            #print len(row)
            if len(row)>8:
                
                
                #print row
                #print row.shape
                try: 
                    j = j+1
                    #print 'File:',file,'Line',j
                    #print row.shape
                    #osgb eastng
                    x = np.asarray(row[0],dtype=np.float64)
                    #osgb northing
                    y = np.asarray(row[1],dtype=np.float64)
                    #osgb /newlyn elevation
                    z = np.asarray(row[2],dtype=np.float64)
                    #intensity derived from LAS tools
                    intensity = np.asarray(row[3],dtype=np.int32)
                    #number of returns identified by lastools
                    #returns = row[4]
                    #and lastly pull out the waveform
                    
                    
                    waveform = np.asarray(row[9:],dtype=np.int16)
                    print 'FWE',waveform.shape, waveform.dtype
                                
                    #smooth the waveform using a univariate spline using the parameters above
                    smoothed2 = smoothing(waveform, kparam, sparam, weights)
                    
                    #print smoothed2
                   
                    
                    #identify the peaks in the smoothed waveform
                    #peaks = signal.find_peaks_cwt(smoothed2, np.arange(19,27,1))
                    
                    #first derivative of smoothed waveform
                    diff = np.diff(smoothed2, n=1)
                    #print diff
                    #second derivative of smoothed waveform
                    #diff2 = np.diff(smoothed2, n=2)
                    
                    #find the maximal value in waveform
                    max_intensity = np.argmax(waveform)
                    #print 'MAX', max_intensity
                    
                    #define the region of the returns
                    diffreg = np.logical_or(diff>1.5,diff<-0.75)
                    #print diffreg
                    #get the x values for slicing the other arrays
                    diffx = x_vals[1:]
                    regx  = diffx[diffreg]
                    #get the first value
                    reg_l = regx[0]
                    #get the last value
                    reg_r = regx[-1]
                    #print 'diffreg', reg_l, reg_r
                    
                    shoulder = np.argmin(diff[reg_l:reg_r])
                    #print 'shoulder pos', shoulder
                    #print shoulder
                            
                    peak_value = waveform[max_intensity]
                    #print peak_value
                            
                    peak_width = reg_r-reg_l
                    print peak_width
                    print reg_r, reg_l
                    print waveform.shape, waveform[reg_l:reg_r].shape
                    print waveform[reg_l:reg_r]
                    peak_sum = np.sum(waveform[reg_l:reg_r])
                    print 'peak sum', peak_sum
                    shoulder_pos = shoulder+reg_l
                    shoulder_int = waveform[shoulder_pos]
                    print type(x)
                    print y
                    vlist = [x,
                             y,
                             z,
                             intensity,
                             reg_l,
                             reg_r,
                             max_intensity,
                             peak_width,
                             peak_value,
                             peak_sum,
                             shoulder_pos,
                             shoulder_int]
                    print vlist
                             
                            
                    os.chdir(outputpath)
                    with open(file, 'a+') as outfile:
                        writer = csv.writer(outfile, delimiter=',')
                        writer.writerow(vlist)
                        outfile.close()
                        os.chdir(datapath)
                        
                    '''waveform_class = classes.classify_point(x,y)
                    
                    wv=[]
                    
                    if waveform_class == 1:
                        os.chdir(outputpath)
                        wv=smoothed2.tolist()
                        wv.insert(0,np.around(y,2))
                        wv.insert(0,np.around(x,2))
                        print wv
                        with open(file+'wv'+'a', 'a+') as outfile:
                            writer = csv.writer(outfile, delimiter=',')
                            writer.writerow(wv)
                            outfile.close()
                            os.chdir(datapath)
                    
                    elif waveform_class == 0:
                        os.chdir(outputpath)
                        wv=smoothed2.tolist()
                        wv.insert(0,np.around(y,2))
                        wv.insert(0,np.around(x,2))
                        with open(file+'wv'+'b', 'a+') as outfile:
                            writer = csv.writer(outfile, delimiter=',')
                            writer.writerow(wv)
                            outfile.close()
                            os.chdir(datapath)'''
                except:
                    continue
                
           

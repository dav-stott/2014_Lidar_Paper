# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 16:37:12 2014

@author: david
"""

'''This module is to classify features in a raster and output a classified raster
of the same array dimensions as the source raster- this makes analysis and
selection of classes within the image straight forward in Numpy ETC.
It does this by:
    * Loading the source raster with GDAL
        - Creating an empty raster of identical dimenions to the source raster
    *Importing the shapefile using OGR
        - Determine unique values for the number of classes
    *Rasterizing each class to a seperate band of the raster using gdal_rasterize
    *** NOTE- THIS IS DONE USING A OS/SYS CALL- YOU NEED TO HAVE THE GDAL BINARIES
    INSTALLLED ON YOUR SYSTEM*** 
        - To do this on Debian / Ubuntu:
             - open a terminal
             -type: 'sudo apt-get install gdal-bin' (you'll need your root password)
        - If you're using a different Linux/Unix distro or OSX consult the iterwebs
        for how to install gdal-bin
        - Not tested on windows. os.path manipulations should mean there's not
        a problem with / vs \. 
    * Other stuff:
        - buffering each feature & outputting it to a raster to map localized
        contrast

Python dependancies you'll need sicpy, numpy, ogr, gdal and shapely'''


import numpy as np
from osgeo import gdal
from osgeo import osr


import numpy.ma as ma
import os


import Image

import fiona
from shapely.geometry import shape
from shapely.geometry import asPoint
from scipy import stats


def writeimage(outpath, 
               outname,
               image,
               spatial):

    data_out = image
    print 'ROWS,COLS',image.shape
    print 'Call to write image'
    os.chdir(outpath)
    print 'OUTPATH',outpath
    print 'OUTNAME',outname
    
    #load the driver for the format of choice
    driver = gdal.GetDriverByName("Gtiff") 
    #create an empty output file
    #get the number of bands we'll need:
    
    bands = 1
    print 'BANDS OUT', bands
    
    #file name, x columns, y columns, bands, dtype
    out = driver.Create(outname, image.shape[1], image.shape[0], bands, gdal.GDT_Float32)
    #define the location using coords of top-left corner
    # minimum x, e-w pixel size, rotation, maximum y, n-s pixel size, rotation
    out.SetGeoTransform(spatial)

    srs = osr.SpatialReference()
    #get the coodrinate system using the ESPG code
    srs.SetWellKnownGeogCS("EPSG:4277")
    #set pstackedstackedstackedtojection of output file 
    out.SetProjection(srs.ExportToWkt())
    
    band = 1
    
    if bands == 1:
        out.GetRasterBand(band).WriteArray(data_out)
        #set the no data value
        out.GetRasterBand(band).SetNoDataValue(-999)
        #apend the statistics to dataset
        out.GetRasterBand(band).GetStatistics(0,1)
        
        print 'Saving %s/%s' % (band,bands)
        
    else:
        while (band<=bands):
            data = data_out[:,:,band-1]
            #write values to empty array
            out.GetRasterBand(band).WriteArray( data )    
            #set the no data value
            out.GetRasterBand(band).SetNoDataValue(-999)
            #apend the statistics to dataset
            out.GetRasterBand(band).GetStatistics(0,1)  
            print 'Saving %s/%s' % (band,bands)
            band = band+1  
    out = None     
    
    print 'Processing of %s complete' % (outname)       
    
    return outname
        

#Class to load the raster image and create an empty raster of the same dimensions
#LoadImage: loads an image using gdal
class LoadImage():
    def __init__(self,infile):
        
        # open the dataset
        self.image_name = infile[:-4]
        self.dataset = gdal.Open(infile) #GA_ReadOnly)
        # if there's nothign there print error
        #self.stacked = None
        if self.dataset is None: 
            print 'BORK: Could not load file: %s' %(infile)
       # otherwise do stuff
        else:
            #get the bit depth of the source image
           
            try:
                pillow_image = Image.open(infile)
                self.bit_depth = pillow_image.bits()
                pillow_image.close()
            except:
                print ('Cant get the bit-depth of the image with pillow')
                
            #get the format
            self.driver = self.dataset.GetDriver().ShortName
            #get the x dimension
            self.xsize = self.dataset.RasterXSize
            #get the y dimension
            self.ysize = self.dataset.RasterYSize
            #get the projection
            self.proj = self.dataset.GetProjection()
            #get the number of bands
            bands = self.dataset.RasterCount
            print 'BANDS:',bands
            #get the geotransform Returns a list object. This is standard GDAL ordering:
                #spatial[0] = top left x
                #spatial[1] = w-e pixel size
                #spatial[2] = rotation (should be 0)
                #spatial[3] = top left y
                #spatial[4] = rotation (should be 0)
                #spatial[5] = n-s pixel size
            self.spatial = self.dataset.GetGeoTransform()

            #print some stuff to console to show  we're paying attention
            print 'Found raster in %s format. Raster has %s bands' %(self.driver,bands)
            print 'Projected as %s' %(self.proj)
            print 'Dimensions: %s x %s' %(self.xsize,self.ysize)
            
            #instantiate a counter
            count = 1
            
            #OK. This is the bit that catually loads the bands in in a while loop
            # Loop through bands as long as count is equal to or less than total
            while (count<=bands):
                print 'BANDS less than COUNT'
                #show that your computer's fans are whining for a reason
                print 'Loading band: %s of %s' %(count,bands)
                #get the band
                band = self.dataset.GetRasterBand(count)
                # load this as a numpy array
                
                #mask the no data values
                data_array = band.ReadAsArray()
                data_array = ma.masked_where(data_array == 0, data_array)
                data_array = data_array.filled(-999)
                data_array = data_array.astype(np.float32, copy=False)
                # close the band object
                band = None
                #this bit stacks the bands into a combined numpy array
                #if it's the first band copy the array directly to the combined one
                if count == 1:
                    self.stacked = data_array
                #else combine these 
                else:
                    self.stacked = np.dstack((self.stacked,data_array))            
                # increment the counter
                count = count+1
                
            #self.coords_matrix = self.coords()
            #print self.coords_matrix.shape
            #print self.coords_matrix
        
    #******************** MEFFODS ******************************************
    def coords(self):
        '''This gets the geographic coordinates of each cell in the raster and 
        returns a list of tuples containing the y and x array references for each pixel
        and the geographic x and y coordinates for each pixel'''
        
        print 'call to coords'
        
        #get the shape of the array
        matrix_dims = self.stacked.shape
        print matrix_dims
        #get the number of rows
        rows = matrix_dims[0]-1
        print rows
        #get the number of columns
        columns = matrix_dims[1]-1
        print columns
        
        x_coords = np.zeros(shape=(matrix_dims[0],matrix_dims[1]))
        y_coords = np.zeros(shape=(matrix_dims[0],matrix_dims[1]))

        
        #instantiate a counter
        row = 0
        
        
        for row in range(matrix_dims[0]):
            #increment counter
            column = 0
            for column in range(matrix_dims[1]):
                #print row, column
                xgeo = self.spatial[0] + (column*self.spatial[1])+(row*self.spatial[2])
                x_coords[row,column] = xgeo
                ygeo = self.spatial[3] + (column*self.spatial[4])+(row*self.spatial[5])
                y_coords[row,column] = ygeo
                column=column+1
        
        print x_coords.shape, y_coords.shape 
        print np.min(x_coords), np.min(y_coords)
        return np.dstack((x_coords,
                          y_coords,
                          np.zeros(shape=(matrix_dims[0],matrix_dims[1]))))
        
    def coord_test(self):
        print 'call to coord test'

        print self.coords_list[0]
        print self.coords_list[-1]


#Class to load the shapefile
#Method to determine the number of unique values in the shapefile
#Map band names to a numeric index
#Method to extract each unique member of a class

#Class to take the features and rasterize them.

class RasterizeFeatures():
    def __init__(self, 
                 shapefile,         # The polygon dataset contianing the classification
                 shape_directory,   # The directory containing the shapefile
                 field_name,        # The name of the field containing the classes
                 image,
                 image_directory,    # The spatial properties of the output array       
                 output_directory):
        #*************Get the classes from the shapefile*******************
        os.chdir(shape_directory)
        #intantiate an empty list
        list_of_classes = []
        
        #load the shapefile using fiona
        shape_reader = fiona.open(shapefile)
        #iterate through the polygons in the shapefile
        for polygon in shape_reader:
            #append the classification values for each polygon to the list
            list_of_classes.append(polygon['properties'][field_name])
        self.classes = list(set(list_of_classes))
        self.classes.sort()
        
        shape_reader.close()
        
        print self.classes
        
        
        #*************** Set up the output raster ******************
        os.chdir(image_directory)
    
        source_image = LoadImage(image)
        
        os.chdir(output_directory)
        
        image_name = source_image.image_name[:-3]

        
        #****************** Rasterize stuff ********************************
                    
        ymin = source_image.spatial[3]-np.abs(source_image.spatial[5]*source_image.ysize)
        xmax = source_image.spatial[0]+source_image.spatial[1]*source_image.xsize
        
        te_string = str(str(source_image.spatial[0])+' '+
                        str(ymin)+' '+
                        str(xmax)+' '+
                        str(np.abs(source_image.spatial[3]))+' ')
        
        
        outname = 'test.tif'        
        
        class_number = 1
        for classification in self.classes:
            if not classification == None:
                print classification
                select_string = str("'"+field_name+'="'+classification+'"'+"'")
                
                print self.classes
                print classification
                print class_number
                print shapefile
                print image_name
                    
                print select_string
                
                if classification == 'ARC':
                    class_number = 1 
                    gdal_string = str('gdal_rasterize -burn '+
                                       str(class_number)+
                                       ' -where '+str(select_string)+' '+
                                       #'-b 1 '+
                                       '-te '+te_string+' '+
                                       '-ts '+str(source_image.xsize)+' '+str(source_image.ysize)+' '+
                                       os.path.join(shape_directory,shapefile)+
                                       ' '+
                                       outname)
                    print gdal_string
                    os.system(gdal_string)
                           
                elif classification == 'BAC':
                    class_number = 2
                    gdal_string = str('gdal_rasterize -burn '+
                                       str(class_number)+
                                       ' -where '+str(select_string)+' '+
                                       '-b 1 '+
                                       os.path.join(shape_directory,shapefile)+
                                       ' '+
                                       outname)
                    print gdal_string               
                    os.system(gdal_string)
                    
                    
                elif classification == 'NAT':
                    class_number = 3
                    gdal_string = str('gdal_rasterize -burn '+
                                       str(class_number)+
                                       ' -where '+str(select_string)+' '+
                                       '-b 1 '+
                                       os.path.join(shape_directory,shapefile)+
                                       ' '+
                                       outname)
                    
                    print gdal_string
                    
                    os.system(gdal_string)
                    
                    
                     
class ZonalStats():
    def __init__(self, class_image, radius=3.3):
        print ()
        self.classified_image = LoadImage(class_image)
        
        x_bounds = self.classified_image.stacked.shape[1]-1
        print 'XBOUNDS',x_bounds
        y_bounds = self.classified_image.stacked.shape[0]-1
        print 'YBOUNDS',y_bounds
        
        pixel_radius = int(radius/((np.abs(self.classified_image.spatial[1])+np.abs(self.classified_image.spatial[5]))/2))
        print ('Pixel Radius: %s ' % (pixel_radius))
        
        
        self.arc_indices = []           
        
        row_ref =0
        for row in self.classified_image.stacked:
            col_ref=0            
            for pixel in row:
                if pixel == 1:
                    min_x = col_ref - pixel_radius
                    if min_x < 0:
                        min_x = 0
                    max_x = col_ref+pixel_radius
                    if max_x > x_bounds:
                        max_x = x_bounds                
                    min_y = row_ref-pixel_radius
                    if min_y < 0:
                        min_y = 0
                    max_y = row_ref+pixel_radius
                    if max_y > y_bounds:
                        max_y = y_bounds  
                    self.arc_indices.append([row_ref,col_ref,min_y,max_y,min_x,max_x])
                    
                col_ref = col_ref+1
            row_ref =row_ref+1
        print 'ARC',len(self.arc_indices)
        
    def stats(self, image):
        
        t_image = np.zeros_like(self.classified_image.stacked)
        p_image = np.zeros_like(self.classified_image.stacked)
        
        for arc_pixel in self.arc_indices:
            print arc_pixel
            array = self.classified_image.stacked[arc_pixel[2]:arc_pixel[3],arc_pixel[4]:arc_pixel[5]]
            #print array
            
            values = image[arc_pixel[2]:arc_pixel[3],arc_pixel[4]:arc_pixel[5]]
            
            #print values
            
            arc = values[np.where(array==1)]
            bac = values[np.where(array==2)]
            nat = values[np.where(array==3)]
            
            if arc.shape[0] > 3:
                if bac.shape[0] > 3:
                    #print arc,bac
                    
                    #print 'ttest'
                    t_arc_bac = stats.ttest_ind(arc,bac, equal_var=False)
                    
                    #print t_arc_bac
                    
                    t_image[arc_pixel[0],arc_pixel[1]] = t_arc_bac[0]
                    p_image[arc_pixel[0],arc_pixel[1]] = t_arc_bac[1]
            
            '''if nat.shape[0] > 0:
                t_arc_nat = stats.ttest_ind(arc,nat, equal_var=False)
                t_bac_nat = stats.ttest_ind(nat,bac, equal_var=False)'''
                
            
        return np.dstack((t_image,p_image))
            
            
            
            
        
            
            
if __name__ == "__main__": 
    root_path = os.path.dirname(os.path.abspath('...'))
    
    print root_path
    
    for directory in os.listdir(root_path):
       print directory
       dir_path = os.path.join(root_path,directory)
       if os.path.isdir(dir_path):
            print dir_path
            image_dir = os.path.join(dir_path, 'image')
            print image_dir
            shp_dir = os.path.join(dir_path, 'shp_dir')
            print shp_dir
            outdir = os.path.join(dir_path, 'grid_output')
            print outdir
            if not os.path.exists(outdir):
                os.mkdir(outdir)
                
            shplist = os.listdir(shp_dir)
            print shplist
            for item in shplist:
                if item[-3:]=='shp':
                    print item
                    shapefile = item
                else:
                    print ('BOO- NO SHAPEFILE')
                 
            img_list = os.listdir(image_dir)
            
            for image in img_list:
                print image
                blah = RasterizeFeatures(shapefile,
                                         shp_dir,
                                         'CLASS',
                                         image,
                                         image_dir,
                                         outdir)
                os.chdir(image_dir)
                image = LoadImage(image)
                spatial=image.spatial
                image =image.stacked
                
                os.chdir(outdir)
                
                for classes in os.listdir(outdir):
                    
                    bands = image.shape[2]
                    
                    
                    
                    for band in range(bands):
                        zonal_blah = ZonalStats(classes)
                        zonal_output = zonal_blah.stats(image[:,:,band])
                        t =zonal_output[:,:,0]
                        p =zonal_output[:,:,1]
                        
                        #np.savetxt('t.txt',t, delimiter = ',')
                        #np.savetxt('p.txt',p, delimiter=',')
                        
                        writeimage(outdir, str(band)+'t.tif', t, spatial)
                        writeimage(outdir, str(band)+'p.tif', p, spatial)
                        
                        t_mean = np.mean(t[np.where(t!=0)])
                    
                        p_mean = np.mean(p[np.where(p!=0)])
                        
                        print ('T:',t_mean, 'P:',p_mean)
                        
                    
    
    
        



# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 13:29:23 2014

@author: david
"""

import os
import csv
import numpy as np
import numpy.ma as ma

import Image

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst

from scipy.ndimage import filters


class InputCSV():
    def __init__(self, indir, outdir, pointfile):
        os.chdir(indir)
        
        #points = []
        
        file_list = os.listdir(indir)
        
        file_counter = 1
        
        #point_total = 0
        
        file_total = len(file_list)
        
        points = np.zeros((0,3))
        
        for file in file_list:
            print ('Reading file %s of %s' % (file_counter, file_total))
            print ('Filename: %s' % (file))
           
            
            
            file_counter = file_counter+1   
            
            flightline = np.genfromtxt(file, delimiter=',', skiprows=1)
            
            points = np.vstack((points,flightline[:,:3]))
            
        points = points[1:,:]          
          
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        os.chdir(outdir)
            
        print ('Saving collated points to: \n %s' % (outdir)) 
        
        #self.points = np.asarray(points)
        
        self.min_x = np.min(points[:,0])
        self.max_x = np.max(points[:,0])
        self.min_y = np.min(points[:,1])
        self.max_y = np.max(points[:,1])    
                    

        np.savetxt(pointfile, points, delimiter=',')
        
        
class GridPoints():
#    Class that grids the input points using a system call to the GDAL library
    def __init__(self, 
                 outdir, #input directory
                 pointfile, #csv file containing points
                 model_dir, #output directory
                 model_name, #otput name
                 z_field, #field to grid
                 min_x, #minimum x
                 max_x, #maximum x
                 min_y, # minimum y
                 max_y, #maximum y
                 spatial_res, #desired spatial resolution
                 EPSG=4277, #EPSG code for the projection
                 method='invdist', #interpolation method
                 power=2.0, #power for interpolation
                 smoothing=0.0, #smoothing constant
                 radius1=1, #radius for interpolation x
                 radius2=1.0, #radius for interpolation y
                 angle=0.0, # angle for radius
                 max_points=9.0, #maximum points to search for
                 min_points=3.0, #minimum points to interpolate
                 nodata=-999, #no data value
                 parallelize=True ):
                
        print ('call to GridPoints')
        

          
        os.chdir(outdir)
        
        # |||||||||||||||| create header for vrt ||||||||||||||||||||||||||||||
        # we need to create a virtual dataset to load the CSV file into OGR
        #To do this we need to make a *.VRT header
        #we'll use the properties of the raster to make it.
        #In the case of files with multiple z values you can specify that field
        # using the z_field variable
        header = str('<OGRVRTDataSource>\n'+
                     '  <OGRVRTLayer name="'+pointfile[0:-4]+'">\n'+
                     '      <SrcDataSource>'+pointfile+'</SrcDataSource>\n'+
                     '      <GeometryType>wkbPoint</GeometryType>\n'+
                     '      <GeometryField encoding= "PointFromColumns" x="field_1" y="field_2" z="field_'+str(z_field)+'"/>\n'+
                     '  </OGRVRTLayer>\n'+
                     '</OGRVRTDataSource>')
                     #% (str(pointfile[0:-3]),str(pointfile), str(z_field)))
        
        print header
        
        #file name
        vrt_name = str(pointfile[0:-4]+'.vrt')
        
        #wirte it
        write_vrt = open(vrt_name, 'w')
        write_vrt.writelines(header)
        write_vrt.close()
        
        # ||||||||||||||| GET  GDAL TO DO STUFF |||||||||||||||||||||||||||||||
        
        #work out the x and y row /column dimensions
        x_size = int((max_x-min_x)/spatial_res)
        y_size = int((max_y-min_y)/spatial_res)
        
        
        print ('X size: %s Y size: %s' % (x_size, y_size))
        
        nodata = str(nodata)
        print nodata
        
        # string containing interpolation parameters
        '''interpolation_string = str(method+':'+
                                   'power='+str(power)+':'+
                                   'smoothing='+str(smoothing)+':'+
                                   'radius1='+str(radius1)+':'+
                                   'radius2='+str(radius2)+':'+
                                   'angle='+str(angle)+':'+
                                   'max_points='+str(max_points)+':'+
                                   'min_points='+str(min_points)+':'+
                                   'nodata='+nodata)'''
                                   
        interpolation_string = str(method+':'+
                                   'radius1='+str(radius1)+':'+
                                   'radius2='+str(radius2)+':'+
                                   'angle='+str(angle)+':'+
                                   'nodata='+nodata)
         
        print interpolation_string
                                   
        #string with the gdal
        gdal_grid_string = str('gdal_grid '+
                               '-ot Float64 '+
                               '-txe '+str(min_x)+' '+str(max_x)+' '+
                               '-tye '+str(min_y)+' '+str(max_y)+' '+
                               '-outsize '+str(x_size)+' '+str(y_size)+' '+
                               '-a_srs '+'EPSG:'+str(EPSG)+' '+
                               '-a '+interpolation_string+' '+
                               '-l '+pointfile[0:-4]+' '+
                               vrt_name+' '+
                               os.path.join(model_dir, model_name))
                               
        print gdal_grid_string                       
        
        #If the parallelize flag is true use all the cores                       
        if parallelize == True:
            # append the gdal config string
            command = str(gdal_grid_string+' --config GDAL_NUM_THREADS ALL_CPUS')
        else:
            # else carry on
            command = gdal_grid_string
        
        command = gdal_grid_string
        
        print ('Making call to GDAL rasterize')
        os.system(command)
        
        print ('Finished Interpolating %s' % (model_name))
                               
class GenerateCropDEMS():
    # this derives the DEMS, doing the bare soil one first and then using the 
    # properties of that to make the crop DEMS
    def __init__(self, bare_path, bare_survey, crop_path, output_dir):
        
        self.list_of_dems =[]        
        
        
        #file name
        raw_file_name = bare_survey+'_points.csv'
        
        #diectory
        bare_raw_dir = os.path.join(bare_path,bare_survey,'raw')
        
        collated = os.path.join(bare_path, bare_survey, 'collated_points')
        if not os.path.exists(collated):
            os.mkdir(collated)
        
        print ('Loading Bare Flightlines')
        
        #call to InputCSV
        bare_points = InputCSV(bare_raw_dir, collated, raw_file_name)
        
        bare_out = os.path.join(output_dir, bare_survey+'_bare')
        if not os.path.exists(bare_out):
            os.mkdir(bare_out)
            
        self.list_of_dems.append(bare_out)    
         
        #get the minima and maxima of the bare fightline
        self.bare_minx = bare_points.min_x
        self.bare_maxx = bare_points.max_x
        self.bare_miny = bare_points.min_y
        self.bare_maxy = bare_points.max_y
        

        
        
        
        
        # grid them
        print ('Gridding Points')
        bare_dem = GridPoints(collated,
                              raw_file_name,
                              bare_out,
                              bare_survey+'.tif',
                              3,
                              self.bare_minx,
                              self.bare_maxx,
                              self.bare_miny,
                              self.bare_maxy,
                              0.8)
                              
        bare_dem
        
        
        # do the crop stuff
        for survey in os.listdir(crop_path):
            crop_raw_dir = os.path.join(crop_path, survey, 'raw')
            crop_file_name = survey+'_crop.csv'
            
            collated = os.path.join(crop_path, survey,'collated_points')
            if not os.path.exists(collated):
                os.mkdir(collated)
            
            crop_dems = os.path.join(output_dir, 'crop')
            if not os.path.exists(crop_dems):
                os.mkdir(crop_dems)
            
            crop_out = os.path.join(crop_dems,survey)
            if not os.path.exists(crop_out):
                os.mkdir(crop_out)
                
            self.list_of_dems.append(crop_dems)
            
            print ('Loading Crop Flightlines')
            crop_points = InputCSV(crop_raw_dir, collated, crop_file_name)
            crop_points
            
            print ('grididng Points')
            crop_dem = GridPoints(collated,
                                  crop_file_name,
                                  crop_out,
                                  survey+'.tif',
                                  3,
                                  self.bare_minx,
                                  self.bare_maxx,
                                  self.bare_miny,
                                  self.bare_maxy,
                                  0.8)
            crop_dem
                       
             
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
        
        #fruity loops
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

                            
        
class CropHeightModel():   
    def __init__(self, bare_dir, dem_dir, outdir):
        
        #loads the GDAL derived DEMs to produce the crop height model
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        os.path.join(dem_dir,bare_dir)
        
        bare = None
        crop_model = None
        
        for bare_dem in os.listdir(bare_dir):
            os.chdir(bare_dir)
            bare = LoadImage(bare_dem)
            
        bare_filtered = filters.median_filter(bare.stacked,size=(9,9))
        bare_spatial = bare.spatial
        
        for survey in os.listdir(dem_dir):
            im_path = os.path.join(dem_dir, survey)
            os.chdir(im_path)
            for image in os.listdir(im_path):
                crop = LoadImage(image)
                crop_filtered = filters.median_filter(crop.stacked,size=(3,3))
                
                crop_model = crop_filtered-bare_filtered

                
                
                bare_mask = ma.masked_where(bare_filtered==-999, bare_filtered)
                crop_mask = ma.masked_where(crop_filtered==-999, crop_filtered)
                #combine these into a single mask so that anything masked in either dataset
                #will be masked in the combined output
                combined_mask = ma.mask_or(bare_mask.mask, crop_mask.mask)
                #use this mask t omask the crop model
                crop_model = ma.array(crop_model, mask=combined_mask)
                
                #convert back to a bog-standard numpy array and fill the masked values
                crop_model = crop_model.filled(-999)
                
                
                name = survey+'_CropHeightModel.tif'
                
                self.writeimage(outdir,
                                name,
                                crop_model.shape[0],
                                crop_model.shape[1],
                                crop_model,
                                bare_spatial)
                
    def writeimage(self,
                   outpath, 
                   outname,
                   xsize,
                   ysize,
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
           
         

if __name__ == "__main__": 
    dir_path = os.path.dirname(os.path.abspath('...'))
    
    bare_base_dir = os.path.join(dir_path,'bare')
    crop_base_dir = os.path.join(dir_path,'crop')
    
    outputdir = os.path.join(dir_path, 'outputs')
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    bare_survey = '20120323'
    #bare_path, bare_survey, crop_path, crop_survey, output_dir
    dems = GenerateCropDEMS(bare_base_dir, bare_survey, crop_base_dir, outputdir)
    dems
    
    CropHeightModel(dems.list_of_dems[0],dems.list_of_dems[1],outputdir)
 

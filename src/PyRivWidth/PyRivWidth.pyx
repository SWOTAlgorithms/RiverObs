"""
Extract river width and height from interferometric water file.

This version builds on several open source projects:

- Pyhton: http://python.org for the python language.
- Cython: http://cython.org for the fast compiled verson of python.
- Numpy: http://scipy.org for array handling.
- Scipy: http://scipy.org for convolve and the dilate and distance_transform_edt morphological operators.
- EVE: http://csb.essex.ac.uk/software/eve.html for the label_regions code, which was subsequently modified by E. Rodriguez
- mahotas: http://luispedro.org/software/mahotas for the thin algorithm

In addition, the centerline code is based on Tamlin Pavelsky's Rivwidth IDL program.
This version is based on Rivwidth_v03_r1.pro.

See the rivwidth site: 

Author: Ernesto Rodriguez

"""

import numpy as N
cimport numpy as N
from scipy.ndimage.morphology import binary_dilation, distance_transform_edt
from scipy.ndimage import convolve
#from scipy.spatial import KDTree
from scipy.spatial import cKDTree as KDTree
from scipy.interpolate import splrep, splev
from mahotas import thin
from Scientific.IO.NetCDF import NetCDFFile
from label_regions import label_regions

class RivWidthHeight:
    """Extract river width and height from interferometric water file."""

    input_params = ['mask_file','height','width','spatial_resolution',
                    'lat_min','dlat','lon_min','dlon',]

    def __init__(self,rundata,fname=None):
        """Initialize with either the input file name or a dictionary containing the run selected options."""

        self.rundata = rundata.copy()
        
        if fname != None: # read from a file
            nc = NetCDFFile(fname)
            shape = nc.variables['water_class'].shape
            self.rundata['mask_file'] = fname
            self.rundata['height'] = shape[0]
            self.rundata['width'] = shape[1]
            self.rundata['spatial_resolution'] = float(nc.spatial_resolution)
            self.rundata['lat_min'] = float(nc.lat_min)
            self.rundata['dlat'] = float(nc.dlat)
            self.rundata['lon_min'] = float(nc.lon_min)
            self.rundata['dlon'] = float(nc.dlon)
            nc.close()
            
        for param in self.input_params:
            if param not in self.rundata:
                raise Exception("Parameter: '%s' not in input data"%param)

        self.mask_file  =self.rundata['mask_file']
        self.height  =self.rundata['height']
        self.width  =self.rundata['width']
        self.spatial_resolution  =self.rundata['spatial_resolution']
        self.lat_min  =self.rundata['lat_min']
        self.dlat  =self.rundata['dlat']
        self.lon_min  =self.rundata['lon_min']
        self.dlon  =self.rundata['dlon']

    def make_channel_mask(self,output_file=None):
        """Develop a channel mask for the desired river system.

        Output:

        chan_maskr: original channel mask (1=water, 0=land)
        chan_mask: channel mask with small water regions erased
        """

        cdef:
            int i, dcnt
            int width, height, resolution, braidflag
            int y
            
        width = self.rundata['width']
        height = self.rundata['height']
        braidflag = self.rundata['braidnarrowflag']
                
        # 1.3.1) Reads in a classified binary water mask mask: this mask can contain water bodies
        # outside of the desired river, but these bodies must be separated from the river by at
        # least one pixel.

        nc = NetCDFFile(self.rundata['mask_file'],'r+')
        chan_mask = nc.variables['water_class'][:]
 
        chan_maskr = chan_mask.copy()

        # Find all connected regions and process only the largest one
        # 1.3.2) Removes all water bodies except for the largest one.
        # If there is a body of water in the image of greater extent than the desired river network,
        # it will have to be removed manually prior to this stage.
        
        if braidflag != 1: # use only 4 nearest neighbors, remove it for simple channel rivers.
            chan_mask = self.erase_small_regions(chan_mask,con8=False)
        else: # use all 8 neighbors for braiding, remove it for simple channel rivers.
            chan_mask = self.erase_small_regions(chan_mask,con8=True)
  
        # 1.3.3) Outputs the final channel mask.  If the original channel mask is already clean
        # (i.e. like the input for prior version of RivWidth), then it will be identical to the
        # input water mask.

        if output_file != None:
            if not 'chan_mask' in nc.variables:
                nc.createVariable('chan_mask','b',('lat','lon'))
            nc.variables['chan_mask'][:] = chan_mask[:]
            nc.sync()
            nc.close()
  
        # 1.3.4) Optional section that exports the original water mask with boundary pixels
        # highlighted around all water bodies.
        # imgbnd=chan_maskr
        # operator=[[1,1,1],[1,0,1],[1,1,1]]
        # imgbnd(where(chan_maskr eq 0 and convol(chan_maskr,operator) gt 0))=2
        # 
        # openw, 2, "/Volumes/Data/files/SWOT/Paris/Ohio/Ohio_bound_image"
        # writeu, 2, imgbnd
        # close, 2

        return chan_maskr, chan_mask

    def make_river_mask(self,chan_maskr,chan_mask,int max_iter = 5,output_file=None):
        """
        1.4)  Removes all islands from within the mask and uses a series of erosions and dilations
        to develop the river mask, from which the river center line will be computed.

        Inputs:

        chan_maskr: 2D array with 0 labeling water and 1 labeling land, but with small regions remaining
        chan_mask: 2D array with 0 labeling water and 1 labeling land, with all small regions removed

        Outputs:

        riv_mask: 2D byte array containing 1 for water and 0 for land, with islands removed

        """

        cdef:
            int i, p, dcnt, height, width, braidflag, max_label

        braidflag = self.rundata['braidnarrowflag']
        height = chan_maskr.shape[0]
        width = chan_maskr.shape[1]

        # this is a 4-connect operator that is used for dilations and erosions
        
        operator=N.array([[0,1,0],[1,0,1],[0,1,0]]) 
        
        # 1.4.1)  Begins with the original water mask and removes all water bodies except the largest,
        # which is presumed to be the river/river network of interest.

        if braidflag == 1: # If we are in a braided system, then we perform 2 dilations before removing islands
                           # to avoid removing key channel pixels that may be disconnected from the main river
            for i in range(0,2):
                chan_maskr = binary_dilation(chan_maskr,structure=operator).astype(chan_maskr.dtype)

            # erase small regions
            chan_maskr = self.erase_small_regions(chan_maskr,con8=True)
            
        else: # not a braided, take the chan_mask, which has small regions erased
            chan_maskr = chan_mask

        # flipim is a version of the binary mask in which water=0 and land=1.
        # It is used to label and remove islands in the river.

        flipim = N.where(chan_maskr == 0, 1, 0).astype(N.int8)
        print "flipim calculated"
    
        labelflip = label_regions(flipim, con8=False)

        print "labelflip calculated"

        # 1.4.2)  Loop that first detects all land areas that intersect the edge of the image.
        # These are retained as "land" in the river mask.
        # All other land polygons are considered to be islands in the river and their values are changed
        # accordingly.  With each loop iteration, we dilate the river area using a 4-connect operator in
        # order to fully bound islands that are incompletely bounded by water in the original 
        # channel mask.  At the moment, we do this process max_iter times, but this number is arbitrary
        # and can be modified on input.

        p = 1
        dcnt = 0
  
        while p == 1:

            print "Iteration:",dcnt

            # Label all the land areas

            max_label = labelflip.max()
            print 'max_label',max_label
##            labelhist = N.histogram(labelflip,  bins=N.arange(max_label+1))[0]

            # Find the labels that occur at the edges (-1 pixel to avoid edge effects)

            side1 = labelflip[0:height,1]
            side1hist = N.histogram(side1, bins=N.arange(max_label+2) )[0] 

            side2 = labelflip[0:height,width-2]
            side2hist = N.histogram(side2, bins=N.arange(max_label+2) )[0] 

            side3 = labelflip[1,0:width]
            side3hist = N.histogram(side3, bins=N.arange(max_label+2) )[0]

            side4 = labelflip[height-2,0:width]
            side4hist = N.histogram(side4, bins=N.arange(max_label+2) )[0]

            print "histograms calculated"

##            if dcnt == 1: return labelflip,(side1,side3,side3,side4),(side1hist,side3hist,side3hist,side4hist)
            
            # This array contains all the region labels that occur at the edges
    
            sidespoly = N.unique(
                N.concatenate( (N.flatnonzero(side1hist != 0),
                                N.flatnonzero(side2hist != 0),
                                N.flatnonzero(side3hist != 0),
                                N.flatnonzero(side4hist != 0) ) )
                )

            print "sidespoly",sidespoly

            # Turn all the land in one of the regions that intersect the edges into -1
    
            for y in range(0,len(sidespoly)):
                if  sidespoly[y] != 0:  # make sure you are not labeling the water
                    index = labelflip == sidespoly[y]
                    if index.any():
                        labelflip[index] = -1

            print "edge land found"

            # Make sure that all of the pixels are either water, or land that intersects the edges
            
            index0 = labelflip > 0 # this should be the index of the islands 
            if index0.any(): # if any islands left?

                print 'some islands still left'
                
                # Turn labelflip into a binary mask for land
                
                labelflip[index0] = 0 # turn those into water
                labelflip[labelflip == -1] = 1 # label the remaining pixels as land

                riv_maskfl = labelflip.copy() # Save the flipped river mask for later use
                riv_mask = N.abs(labelflip-1).astype(N.int8) # invert land and water to get a new river mask

                # Grow the water

                riv_maskd = binary_dilation(riv_mask,structure=operator).astype(riv_mask.dtype)
                flipim = N.where(riv_maskd == 0, 1, 0).astype(N.int8) # This is the new land mask
                labelflip = label_regions(flipim, con8=False) # and these are the new land regions
                dcnt += 1
            else:
                p = 0

            # Stop after the maximum number of iterations has been reached

            if dcnt > max_iter: p = 0


        # This looks like just another dilation, but in fact we are dilating the opposite direction
        # (ie eroding "river" areas) because we are applying the dilations to a 
        # flipped version of the image where water=0 and land=1.  We do as many "erosions" as we have
        # done dilations prior to this point.

        if dcnt > 2:
            for i in range(dcnt-2):
                riv_maskfl = binary_dilation(riv_maskfl,structure=operator).astype(riv_maskfl.dtype)
        ## else:
        ##     raise Exception("Error in computing river mask.")
  
        # Rivmask is the final river mask to be used in calculating the river center line.

        riv_mask = N.where(riv_maskfl == 0, 1, 0).astype(N.int8)

        # 1.5) Outputs the final river mask and final channel mask to file locations specified
        # in the input parameter file.
        
        if output_file != None:
            nc = NetCDFFile(self.rundata['mask_file'],'r+')
            if not 'riv_mask' in nc.variables:
                nc.createVariable('riv_mask','b',('lat','lon'))
            nc.variables['riv_mask'][:] = riv_mask[:]
            nc.sync()
            nc.close()

        return riv_mask

    def erase_small_regions(self,mask,con8=True,copy=False):
        """Given a binary mask, find all connected regions, and keep only the largest one.

        mask: input binary mask
        con8: if True, use all 8 neighbors, otherwise, use only 4 nearest neighbors
        copy: if True, a new array is returned. Otherwise, the original array is modified in place.
        """

        label = label_regions(mask, con8=con8)
        labelhist = N.histogram(label, bins=N.arange(label.max()+2))[0]
        maxhist = labelhist[1:].max()
        maxhist_index = N.flatnonzero(labelhist == maxhist)[0]
        eraser = label != maxhist_index
        if copy:
            return N.where(eraser,0,mask)
        
        mask[eraser] = 0

        return mask
        
    def center_line_calc(self,img,deriv_threshold=0.85,output_file=None,use_thin=False):
        """
        center_line_calc is a function, that calculates the initial centerline for the river channel.
        This is based either on edge-detection methods detailed in Pavelsky and Smith (2008, IEEE GRSL),
        or on image thinning.
        """

        if output_file != None:
            nc = NetCDFFile(self.rundata['mask_file'],'r+')
            
        if use_thin:
            center_line_mask = thin(img)
        else:
            width = img.shape[0]
            height = img.shape[1]

            # 1) Create an image (imgbnd) where the nonriver pixels that are adjacent to river pixels have a
            # value of 2, and all other pixels have a value of 0.

            operator = N.array([[1,1,1],[1,0,1],[1,1,1]])
            imgbnd = N.empty(img.shape,dtype=img.dtype)

            bnd_index = ( img == 0 ) & ( convolve(img,operator) > 2 )
            imgbnd[bnd_index] = 2
            imgbnd[imgbnd < 2] = 0

            # 2) Calculate the Euclidean distance map using the scipy  distance_transform_edt

            mask_distance = distance_transform_edt(img) # Exact Euclidean distance in python


            # 6) With the true euclidean distance map computed, we can now calculate the centerline.  

            # 6.1) The horizontal and vertical laplacian operators
            # NOTE: THIS COULD BE SPEEDED UP WITH A LOOP OF 1D CONVOLUTIONS.

            operator_hor = N.array([[0.,0.,0.],
                                   [-0.5,0,0.5],
                                   [0.,0.,0.]])

            operator_vert = N.array([[0.,0.5,0.],
                                    [0.,0,0.],
                                    [0.,-0.5,0.]])


            # 6.2) The two operators above are convolved with the distance map to produce a laplacian map.
            # River values will be close to 1 except near the river centerline, where they will be near 0.

            mask_distance_derivative = ( convolve(mask_distance, operator_hor)**2 +
                                         convolve(mask_distance, operator_vert)**2 )

            # 6.3) The following steps extract the centerline from the raw laplacian image and refine it by
            # removing any centerline pixels not attached to the main centerline.  The value 0.85 can be edited 
            # to provide a more or less lenient standard for selecting centerline pixels.
            # Editing it downward (e.g. 0.80) will exclude more pixels, and increasing it (e.g. 0.90) will
            # include more.  At the moment, everything that doesn't connect to the largest centerline segment
            # will be erased during the cleanup process.  As such, if you find that your centerline stops at
            # an incorrect point, I recommend revising this value up slightly and rerunning.

            MD_Cline = mask_distance_derivative.copy()
            index = ( MD_Cline <= deriv_threshold ) & ( imgbnd == 2)
            MD_Cline[index] = 1
            index = ( MD_Cline > deriv_threshold ) | ( img == 0)
            MD_Cline[index] = 2
            MD_Cline[MD_Cline <= deriv_threshold] = 1.
            MD_Cline[MD_Cline == 2] = 0
            center_line_mask = self.erase_small_regions(MD_Cline.astype(N.int8), con8=True,copy=True)
            center_line_mask[center_line_mask != 0] = 1

            # Save the results, if desired

            if output_file != None:
                if not 'center_line_mask' in nc.variables:
                    nc.createVariable('center_line_mask','b',('lat','lon'))
                nc.variables['center_line_mask'][:] = center_line_mask.astype(N.int8)

                if not 'mask_distance' in nc.variables:
                    nc.createVariable('mask_distance','f',('lat','lon'))
                nc.variables['mask_distance'][:] = mask_distance.astype(N.float32)

                if not 'mask_distance_derivative' in nc.variables:
                    nc.createVariable('mask_distance_derivative','f',('lat','lon'))
                nc.variables['mask_distance_derivative'][:] = mask_distance_derivative.astype(N.float32)

        if output_file != None:
            if not 'center_line_mask' in nc.variables:
                nc.createVariable('center_line_mask','b',('lat','lon'))
            nc.variables['center_line_mask'][:] = center_line_mask.astype(N.int8)

            nc.sync()
            nc.close()

        # 7) Returns the finished centerline
        
        return center_line_mask

    def dynamic_vector_v3(self, input_image, int srow, int scol,  int erow, int ecol):
        """
        Begin Dynamic_Vector Code:  Dynamic vector uses a minimum cost search algorithm to determine
        the shortest path from the start point (srow,scol) to the end point (erow,ecol) contained within
        the initial centerline.  This should not be edited except in very unusual circumstances.
        """
        cdef:
            int count, cost, signal, x, ec
            int width, height
            N.ndarray[N.int8_t, ndim=2] step
            N.ndarray[N.int16_t, ndim=3] cim
            N.ndarray[N.int8_t, ndim=2] out_image
            N.ndarray[N.int16_t, ndim=2] queue
            N.ndarray[N.int16_t, ndim=1] cur
            N.ndarray[N.int16_t, ndim=2] center_line_indexes
            N.ndarray[N.int16_t, ndim=1] larray

        # Get the inputs

        width = int(self.rundata['width'])
        height = int(self.rundata['height'])

        # Get the points closest to the starting guess

        srow,scol,erow,ecol = self.search_nearest_mask(input_image,srow,scol,erow,ecol,mask_value=1)
        print 'srow',srow,'scol',scol,'input_image[srow,scol]',input_image[srow,scol]
        print 'srow',erow,'scol',ecol,'input_image[erow,ecol]',input_image[erow,ecol]

        ## # Sanity check: have to start at at valid connect water pixel

        ## if input_image[srow,scol] != 1:
        ##     print 'srow',srow,'scol',scol,'input_image[srow,scol]',input_image[srow,scol]
        ##     raise Exception('Starting or ending guess is not in water')

        # "I can't remember why I have to do this, but it turns out right if I do and doesn't work if I don't.
        # That's the problem with writing code over the course of 8 years..." Tamlin Pavelsky. FIX THIS!!
        
        ## erow = erow-1
        ## ecol = ecol-1

        ## if input_image[erow,ecol] != 1:
        ##     print 'erow',erow,'ecol',ecol,'input_image[erow,ecol]',input_image[erow,ecol]
        ##     raise Exception('Ending guess is not water')

        # 1) Variables and arrays initialized that will be used in the loop below:
        
        count = 0 # initializes number of values in queue
        cost = 0 #initializes the cost value
        signal = 0 # a signal variable that is switched to 1 when the destination pixel is reached
        # an array that stores all centerline pixels that have already been searched
        step = N.zeros((height,width),dtype=N.int8)
        # an array thot holds the x/y coordinates of the adjacent centerline pixel from
        # which the pixel in question was reached.
        cim = N.zeros((height,width,2),dtype=N.int16)
        # the image that will contain the final, 1-pixel centerline
        out_image = N.zeros((height,width),dtype=N.int8)
        # the queue that contains the pixels to be searched for (erow,ecol) 
        queue = N.zeros((height*2,5),dtype=N.int16) 

        # 2) Initialize the queueu by adding the start pixel

        # initializes the first pixel in the queue with the start pixel
        queue[count,:] = N.array([srow,scol,cost,0,0],dtype=queue.dtype)  
        step[srow,scol] = 1
        count = count + 1

        # 3) Each iteration of the following loop extracts a pixel from the queue and determines if
        # it is the end pixel (erow, ecol). If it is, then we can output our final centerline.
        # If not, then we find all of the pixels in an eight-connect neighborhood that have not
        # already been added to the queue and add them to the queue.  

        print 'start search', signal, count
        while ( signal != 1 )  and ( count != 0 ):
  
            cur = queue[0,:].copy() # extracts the first pixel in the queue)
            print 'first pixel in queue',cur

            # keeps track of the prior pixel coordinates so we can determine our path back to (srow,scol)
            cim[cur[0],cur[1],0] = cur[3]  
            cim[cur[0],cur[1],1] = cur[4]
            
            # resets the queue by removing the pixel we are about to consider and shifting all other
            # pixels up one slot   
            queue[0:count,:] = queue[1:count+1,:] 
            count = count-1 # keeps track of the queue length
            print 'first pixel in queue',cur
            print 'elements in queue',count

            # If we have reached our end pixel (erow,ecol), then this runs
            if ( cur[0] == erow ) and ( cur[1] == ecol ):

                print 'reached end'
                
                signal = 1  # sets signal variable to stop our while loop
                # assigns a value of 1 to our end pixel (erow,ecol) in the output image
                out_image[cur[0],cur[1]] = 1
                # initializes our output array, which just contains the x,y coordinates
                # of the centerline pixels
                center_line_indexes = N.zeros((cur[2],2),dtype=N.int16)

                # initializes larray, which allows us to go back through our search and determine
                # the fastest route back to (srow,scol) from (erow,ecol)
                larray=N.array([cur[3],cur[4]],dtype=N.int16)
                
                for x in range(cur[2]):
                    center_line_indexes[cur[2]-x-1,0] = larray[0]
                    center_line_indexes[cur[2]-x-1,1] = larray[1]
                    out_image[larray[0],larray[1]] = 1
                    larray = N.array([cim[larray[0],larray[1],0],cim[larray[0],larray[1],1]],
                                     dtype=N.int16)
                    
            else: # If we have not yet reached our end pixel (erow,ecol), then this runs

                print 'have not reached end'

                # determines the number of initial centerline pixels in an
                # 8-connect neighborhood around our current pixel
                print 'image neighborhood around',cur[0],cur[1]
                print input_image[cur[0]-1:cur[0]+2,cur[1]-1:cur[1]+2]
                ec = self.eight_vector(cur[0],cur[1],input_image, width, height)
                print 'ec',ec
                print 'cur',cur
                
                if ( ( ec >= 1 ) and ( cur[0] != 0 ) and ( cur[1] != 0 ) and
                     ( cur[0] != height-1 ) and ( cur[1] != width-1 ) ):

                    # gets the coordinates of all pixels in an 8-connect neighborhood with values of 1
                    neigh = self.getpoints(input_image,ec,cur[0],cur[1], cur[2], width, height)
                    print 'neigh',neigh

                    for x in range(ec):
                        print 'x',x
                        print 'step[neigh[x,0],neigh[x,1]]',step[neigh[x,0],neigh[x,1]]
                        # if it hasn't already been added, add each pixel obtained using getpoints to the queue
                        if step[neigh[x,0],neigh[x,1]] == 0:
                            print 'in if'
                            queue[count,:] = N.array([neigh[x,0],neigh[x,1],neigh[x,2],neigh[x,3],neigh[x,4]],
                                                     dtype=N.int16)
                            step[neigh[x,0],neigh[x,1]] = 1
                            count=count+1
                    print 'out of for'
                print 'out of if'
            print 'out of else'
            print 'signal',signal,'count',count
        print 'out of while'

        # 4) Return the final centerline array

        print 'ready to return'
        
        return center_line_indexes

    def eight_vector(self, int srow, int scol, im, int width, int height):
        """determines the number of initial centerline pixels in an
        8-connect neighborhood around our current pixel"""

        if ( srow == 0 ) or ( srow == height-1 ) or ( scol == 0 ) or ( scol == width-1 ):
            eight_conn = 0
            if srow > 0:
                eight_conn += im[srow-1,scol]
                if scol < width - 1:
                    eight_conn += im[srow-1,scol+1]
                if scol > 0:
                    eight_conn += im[srow-1,scol-1]
            if srow < height - 1:
                eight_conn += im[srow+1,scol]
                if scol < width - 1:
                    eight_conn += im[srow+1,scol+1]
                if scol > 0:
                    eight_conn += im[srow+1,scol-1]
            if scol < width - 1:
                eight_conn += im[srow,scol+1]
            if scol > 0:
                eight_conn += im[srow,scol-1]
        else:
            eight_conn = ( im[srow,scol+1] + im[srow,scol-1] + im[srow-1,scol] + im[srow+1,scol] +
                           im[srow-1,scol+1] + im[srow+1,scol+1] + im[srow+1,scol-1] + im[srow-1,scol-1])
        return eight_conn

    def getpoints(self, im, int ec, int cx, int cy, int cost, int width, int height):
        """gets the coordinates of all pixels in an 8-connect neighborhood with values of 1"""

        cdef:
            int n
            
        outpoints = N.zeros(shape=(ec+1,5),dtype=N.int16)
        n=0

        if ( im[cx+1,cy] == 1 ) and ( cx < height - 1 ):
            outpoints[n,:] = N.array([cx+1,cy,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if ( im[cx-1,cy] == 1 ) and ( cx > 0 ):
            outpoints[n,:] = N.array([cx-1,cy,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if ( im[cx,cy+1] == 1 ) and ( cy < width - 1 ):
            outpoints[n,:] = N.array([cx,cy+1,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if im[cx,cy-1] == 1:
            outpoints[n,:] = N.array([cx,cy-1,cost+1,cx,cy],dtype=outpoints.dtype) 
            n += 1

        if ( im[cx+1,cy-1] == 1 ) and  ( cx < height - 1 ):
            outpoints[n,:] = N.array([cx+1,cy-1,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if ( im[cx+1,cy+1] == 1 ) and ( cx < height - 1 ) and ( cy < width - 1 ):
            outpoints[n,:] = N.array([cx+1,cy+1,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if ( im[cx-1,cy-1] == 1 ) and ( cx > 0 ) and ( cy > 0):
            outpoints[n,:] = N.array([cx-1,cy-1,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        if ( im[cx-1,cy+1] == 1 ) and ( cx > 0 ) and ( cy < width - 1 ):
            outpoints[n,:] = N.array([cx-1,cy+1,cost+1,cx,cy],dtype=outpoints.dtype)
            n += 1

        return outpoints

    def search_nearest_mask(self,mask,srow,scol,erow,ecol,mask_value=1):
        """Given a starting guess for points."""

        # Build the array of input data
        
        data = N.array(N.nonzero((mask==mask_value))).transpose()

        # Get the kdtree

        kdt = KDTree(data)

        # Get the nearest points

        distance, index = kdt.query(((srow,scol),(erow,ecol)))

        # get the points

        srow, scol = kdt.data[index[0]]
        erow, ecol = kdt.data[index[1]]

        return srow,scol,erow,ecol

    def get_center_line_parameters(self,center_line_indexes,  ds):
        """Turn the center line data into physical distances from the center tile by using
        the spatial resolution. Initialize the KDTrees structure, the spline and vector fields."""

        cdef:
            int i, nx, ns
            double spatial_resolution
        ##     N.ndarray[N.float64_t, ndim=1] x
        ##     N.ndarray[N.float64_t, ndim=1] y
        ##     N.ndarray[N.float64_t, ndim=1] s

        spatial_resolution = self.rundata['spatial_resolution']
        self.center_line = spatial_resolution*center_line_indexes
        
        # Get the distance along the centerline

        x = self.center_line[:,0].copy()
        y = self.center_line[:,1].copy()
        nx = len(x)
        s = N.zeros((nx,),dtype=x.dtype)

        for i in range(1,nx):
            s[i] = s[i-1] + N.sqrt((x[i] - x[i-1])**2 + (y[i] - y[i-1])**2 )
        
        # Get a spline representation of the center line

        self.x_tck = splrep(s,x)
        self.y_tck = splrep(s,y)

        # Now get the array of desired s values (rather than sampled svalues)

        ns = int( (s[-1] - s[0])/ds + 1)
        self.s_interpolated_center_line = N.zeros((ns,),dtype=x.dtype)
        for i in range(1,ns):
            self.s_interpolated_center_line[i] = self.s_interpolated_center_line[i-1] + ds

        # Get the interpolated center line and derivatives

        self.interpolated_center_line = N.zeros((ns,2),dtype=x.dtype)
        self.interpolated_center_line[:,0] = splev(self.s_interpolated_center_line,self.x_tck)
        self.interpolated_center_line[:,1] = splev(self.s_interpolated_center_line,self.y_tck)
        
        dx_ds = splev(self.s_interpolated_center_line,self.x_tck,der=1)
        dy_ds = splev(self.s_interpolated_center_line,self.y_tck,der=1)
        norm = N.sqrt(dx_ds**2 + dy_ds**2)
        dx_ds /= norm
        dy_ds /= norm

        # Get the tangent and normal vectors

        self.tangent_vector = N.zeros((ns,2),dtype=dx_ds.dtype)
        self.tangent_vector[:,0] = dx_ds
        self.tangent_vector[:,1] = dy_ds

        ## self.normal_vector = N.zeros((ns,2),dtype=dx_ds.dtype)
        ## self.normal_vector[:,0] = -dy_ds
        ## self.normal_vector[:,1] = dx_ds

    def center_line_to_file(self,center_line_file):
        """Write the center line information into a netcdf file."""

        nc = NetCDFFile(center_line_file,'w')

        nraw = self.center_line.shape[0]
        n  = self.interpolated_center_line.shape[0]
        
        nc.createDimension('raw_scoord',nraw)
        nc.createDimension('scoord',n)

        rcenter_line_x = nc.createVariable('raw_center_line_x','d',('raw_scoord',))
        rcenter_line_x[:] = self.center_line[:,0].astype(N.float64)

        rcenter_line_y = nc.createVariable('raw_center_line_y','d',('raw_scoord',))
        rcenter_line_y[:] = self.center_line[:,1].astype(N.float64)

        rlat, rlon = self.xy_to_latlon(self.center_line[:,0],self.center_line[:,1])

        rcenter_line_lat = nc.createVariable('raw_center_line_lat','d',('raw_scoord',))
        rcenter_line_lat[:] = rlat.astype(N.float64)

        rcenter_line_lon = nc.createVariable('raw_center_line_lon','d',('raw_scoord',))
        rcenter_line_lon[:] = rlon.astype(N.float64)
        
        center_line_x = nc.createVariable('center_line_x','d',('scoord',))
        center_line_x[:] = self.interpolated_center_line[:,0].astype(N.float64)

        center_line_y = nc.createVariable('center_line_y','d',('scoord',))
        center_line_y[:] = self.interpolated_center_line[:,1].astype(N.float64)

        S = nc.createVariable('scoord','d',('scoord',))
        S[:] = self.s_interpolated_center_line[:].astype(N.float64)

        tangent_vector_x = nc.createVariable('tangent_vector_x','d',('scoord',))
        tangent_vector_x[:] = self.tangent_vector[:,0].astype(N.float64)

        tangent_vector_y = nc.createVariable('tangent_vector_y','d',('scoord',))
        tangent_vector_y[:] = self.tangent_vector[:,1].astype(N.float64)

        lat, lon = self.xy_to_latlon(self.interpolated_center_line[:,0],
                                     self.interpolated_center_line[:,1])

        center_line_lat = nc.createVariable('center_line_lat','d',('scoord',))
        center_line_lat[:] = lat.astype(N.float64)

        center_line_lon = nc.createVariable('center_line_lon','d',('scoord',))
        center_line_lon[:] = lon.astype(N.float64)

        nc.sync()
        nc.close()

    def get_center_line_kdtree(self):
        """Compute the centerline kdtree based on an interpolated centerline."""
        
        self.center_line_kdtree = KDTree(self.interpolated_center_line)

    def get_s_c_coordinates(self,N.ndarray[N.float64_t, ndim=2] xy):
        """Given an array of (x,y) points, get their coordinates in the river reference frame."""

        cdef:
            int i, ii, nxy
            
        distance, index = self.center_line_kdtree.query(xy)
        nxy = index.shape[0]

        S = N.empty(nxy,dtype=xy.dtype)
        C = N.empty(nxy,dtype=xy.dtype)
        
        get_S_C_coordinates(self.interpolated_center_line,
                            self.tangent_vector,
                            self.s_interpolated_center_line,
                            xy, S, C, index)

        return S, C, index, distance

    def latlon_to_xy(self,
                     N.ndarray[N.float64_t, ndim=1] lat,
                     N.ndarray[N.float64_t, ndim=1] lon):
        """Get the local projection x,y coordinates of a lat/lon array.

        The result is returned as a 2D array suitable for input into get_s_c_coordinates.
        """

        cdef:
            double spatial_resolution
            
        n = len(lat)
        
        xy = N.empty((n,2),dtype=N.float64)

        spatial_resolution = self.rundata['spatial_resolution']
        xy[:,0] = ( (lat - self.lat_min )/self.dlat )*self.spatial_resolution
        xy[:,1] = ( (lon - self.lon_min )/self.dlon )*self.spatial_resolution

        return xy

    def xy_to_latlon(self,
                     N.ndarray[N.float64_t, ndim=1] x,
                     N.ndarray[N.float64_t, ndim=1] y):
        """Go from local projection coordinates to lat/lon."""

        n = len(x)
        
        lat = self.lat_min + (x/self.spatial_resolution)*self.dlat
        lon = self.lon_min + (y/self.spatial_resolution)*self.dlon

        return lat, lon

    def process_file(self,center_line_file,data_file_name,lat_name='lat',lon_name='lon'):
        """Process all the data in a data file and update the file with the estimated
        S and C coordinates."""

        # Read the center line parameters
        
        nc_cl = NetCDFFile(center_line_file)
        clv = nc_cl.variables

        ncl = clv['center_line_x'].shape[0]
        dtype = N.float64 

        self.interpolated_center_line = N.empty((ncl,2),dtype=dtype)
        self.interpolated_center_line[:,0] = clv['center_line_x'][:]
        self.interpolated_center_line[:,1] = clv['center_line_y'][:]

        self.tangent_vector = N.empty((ncl,2),dtype=dtype)
        self.tangent_vector[:,0] = clv['tangent_vector_x'][:]
        self.tangent_vector[:,1] = clv['tangent_vector_y'][:]

        self.s_interpolated_center_line = clv['scoord'][:]

        nc_cl.close()

        print "centerline read"

        # Initialize the kdtree for center line distance

        self.get_center_line_kdtree()

        print "KDTree initialized"

        # Open the output file for updates and read the lat,lon data

        nc = NetCDFFile(data_file_name,'r+')
        var = nc.variables

        lat = var[lat_name][:]
        lon = var[lon_name][:]

        print "lat and lon read"
        
        # Project them into the local coordinate system

        xy = self.latlon_to_xy(lat, lon)

        print "data projected to local coordinates"

        # Get their S and C coordinates

        S, C, index, distance = self.get_s_c_coordinates(xy)

        print "S, C, index calculated"

        # Now output the results

        if not 's_coord' in var:
            nc.createVariable('s_coord','d',('record',))
        var['s_coord'][:] = S[:].astype(N.float64)
        
        if not 'c_coord' in var:
            nc.createVariable('c_coord','d',('record',))
        var['c_coord'][:] = C[:].astype(N.float64)

        if not 'center_line_index' in var:
            nc.createVariable('center_line_index','i',('record',))
        var['center_line_index'][:] = index[:].astype(N.int32)

        if not 'center_line_distance' in var:
            nc.createVariable('center_line_distance','d',('record',))
        var['center_line_distance'][:] = distance[:].astype(N.float64)

        nc.sync()
        nc.close()
        
        print "Data output finished"
        
cdef get_S_C_coordinates(N.ndarray[N.float64_t, ndim=2] center_line,
                         N.ndarray[N.float64_t, ndim=2] tangent_vector,
                         N.ndarray[N.float64_t, ndim=1] S0,
                         N.ndarray[N.float64_t, ndim=2] xy,
                         N.ndarray[N.float64_t, ndim=1] S,
                         N.ndarray[N.float64_t, ndim=1] C,
                         N.ndarray[N.int32_t, ndim=1] index                         
                         ):
    cdef:
        int i, ii, nxy
        double x0, y0, dx_ds, dy_ds, x, y
        double s, c, s0

    nxy = index.shape[0]
    
    for i in range(nxy):
        ii = index[i]
        x0 = center_line[ii,0]
        y0 = center_line[ii,1]
        dx_ds = tangent_vector[ii,0]
        dy_ds = tangent_vector[ii,1]
        x = xy[i,0]
        y = xy[i,1]
        s, c = get_s_c(x0,y0,dx_ds,dy_ds, x, y)
        S[i] = S0[ii] + s
        C[i] = c
        
cdef get_s_c(double x0, double y0, double dx_ds, double dy_ds, double x, double y):
    """Given a local origin for a coordinate axis (x0, y0) and a normal defining the
    the s axis, return the coordinates of a target point (x,y) in the (s,c)
    coordinate system, where c is perpendicular to s."""

    cdef:
        double s, c, dx, dy

    dx = (x - x0)
    dy = (y - y0)
    s = dx_ds*dx + dy_ds*dy
    c = -dy_ds*dx + dx_ds*dy

    return s, c
        

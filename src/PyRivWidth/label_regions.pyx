"""
Given a segmented image, return an image with its regions labelled.
This code is adapted from the label_regions code written by
Adrian F. Clark <alien@essex.ac.uk> for EVE

http://csb.essex.ac.uk/software/eve.html

and
 
http://csb.essex.ac.uk/software/eve.py.txt

and distributed entirely freely.

Modifications made by Ernesto Rodriguez:

1) Translated the code to cython to improve efficiency.
2) Removed 3rd (unused) dimension in inputs and outputs
3) Added the possibility of treating all pixels labeled by 0 as belonging to the same region.
This last behavior matches the behavior of the IDL LABEL_REGION function.
"""

import numpy
cimport numpy

def  label_regions(im, con8=True, zero_regions=False):
    """
    Given a segmented image, return an image with its regions labelled.

    Arguments:
    im -- image to be labelled
    con8 -- if True, consider all 8 nearest neighbours
            if False, consider only 4 nearest neighbours
    zero_regions -- if False, all 0 values are considered to be in the 0 region (IDL behavior),
                    else, 0 values can form different regions.

    Returns: labels in a two-dimensional int32 array
    """

    if con8: int_con8 = 1
    else: int_con8 = 0
    if zero_regions: int_zero_regions = 1
    else: int_zero_regions = 0
        
    # call fast cython routines, based on data type

    if im.dtype == numpy.int8:
        return label_regions_int8 (im, int_con8,int_zero_regions)
    if im.dtype == numpy.uint8:
        return label_regions_uint8 (im, int_con8,int_zero_regions)
    elif im.dtype == numpy.int32:
        return label_regions_int32 (im, int_con8, int_zero_regions)
    else:
        try:
            return label_regions_uint8 (im.astype(numpy.unit8), int_con8,int_zero_regions)
        except:
            raise Exception("label_regions not implemented for dtype: %s"%im.dtype)

    
cdef label_regions_uint8 (numpy.ndarray[numpy.uint8_t, ndim=2] im, int con8,int zero_regions):
    """
    Given a segmented image, return an image with its regions labelled.

    Arguments:
    im -- image to be labelled
    con8 -- if True, consider all 8 nearest neighbours
            if False, consider only 4 nearest neighbours

    Returns: labels in a two-dimensional int32 array            
    """
    cdef:
        int nx, ny, x, y
        int lastcol, inreg
        numpy.ndarray[numpy.int32_t, ndim=2] lab
        numpy.ndarray[numpy.int32_t, ndim=1] vals
        numpy.ndarray[numpy.int32_t, ndim=1] labs

    # For fast evaluation of if statements by cython
    
    ny = im.shape[0]
    nx = im.shape[1]
    lab = numpy.zeros((ny, nx), dtype=numpy.int32)
    vals = numpy.array([0, 0, 0, 0])
    labs = numpy.array([1, 0, 0, 0])

    # The upper left pixel is in region zero.
    lastlabel = 0
    equiv = [lastlabel]
    if zero_regions == 0:
        if im[0,0] != 0:
            lastlabel += 1
            equiv.append (lastlabel)
            lab[0,0] = lastlabel
    else:
        lab[0,0] = 0
        
    # Process the rest of the first row of the image.
    y = 0
    if ( zero_regions == 0 ):
        for x in range (1, nx):
            if ( im[y,x] != 0 ):
                if im[y,x] != im[y,x-1]:
                    lastlabel += 1
                    equiv.append (lastlabel)
                lab[y,x] = lastlabel
    else:
        for x in range (1, nx):
            if im[y,x] != im[y,x-1]:
                lastlabel += 1
                equiv.append (lastlabel)
            lab[y,x] = lastlabel

    # Process the first column of the image.
    x = 0
    if ( zero_regions == 0 ):
        for y in range (1, ny):
            if im[y,x] != 0:
                if im[y,x] == im[y-1,x]:
                    lv = lab[y-1,x]
                else:
                    lastlabel += 1
                    equiv.append (lastlabel)
                    lv = lastlabel
                lab[y,x] = lv
    else:
        for y in range (1, ny):
            if im[y,x] == im[y-1,x]:
                lv = lab[y-1,x]
            else:
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            lab[y,x] = lv
                
    # Process the remainder of the image.
    for y in range (1, ny):
        y1 = y - 1
        for x in range (1, nx):
            if con8 == 1: nv = 4
            else:    nv = 2
            x1 = x - 1
            x2 = x + 1
            if ( x2 >= nx - 1 ) and ( con8 == 1 ): lastcol = 1; nv -= 1
            else: lastcol = 0

            # Do not process if this is a 0 pixel (goes to label 0)
            
            if ( im[y,x] == 0 ) and (zero_regions == 0): continue
            
            val = im[y,x]
            # Get the four neighbours' values and labels, taking care not
            # to index off the end of the image.
            vals[0] = im[y, x1]; labs[0] = lab[y, x1]
            vals[1] = im[y1,x]; labs[1] = lab[y1,x]
            if con8 == 1:
                vals[2] = im[y1,x1]; labs[2] = lab[y1,x1]
                if lastcol == 0:
                    vals[3] = im[y1,x2]; labs[3] = lab[y1,x2]
            inreg = False
            for i in range (0, nv):
                if val == vals[i]: inreg = 1
            if inreg == 0:
                # We're in a new region.
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            else:
                # We must be in the same region as a neighbour.
                matches = []
                for i in range (0, nv):
                    if val == vals[i]: matches.append (labs[i])
                matches.sort ()
                lv = int(matches[0])
                for v in matches[1:]:
                    if equiv[v] > lv:
                        equiv[v] = lv
                    elif lv > equiv[v]:
                        equiv[lv] = equiv[v]
            lab[y,x] = lv

    # Tidy up the equivalence table.
    remap = list()
    nc = -1
    for i in range (0, len(equiv)):
        if equiv[i] == i:
            nc += 1
            v = nc
        else:
            v = i
            while equiv[v] != v:
                v = equiv[v]
            v = remap[v]
        remap.append (v)

    # Make a second pass through the image, re-labelling the regions, then
    # return the labelled image.
    for y in range (0, ny):
        for x in range (0, nx):
            lab[y,x] = remap[lab[y,x]]
    return lab

cdef label_regions_int8 (numpy.ndarray[numpy.int8_t, ndim=2] im, int con8,int zero_regions):
    """
    Given a segmented image, return an image with its regions labelled.

    Arguments:
    im -- image to be labelled
    con8 -- if True, consider all 8 nearest neighbours
            if False, consider only 4 nearest neighbours

    Returns: labels in a two-dimensional int32 array            
    """
    cdef:
        int nx, ny, x, y
        int lastcol, inreg
        numpy.ndarray[numpy.int32_t, ndim=2] lab
        numpy.ndarray[numpy.int32_t, ndim=1] vals
        numpy.ndarray[numpy.int32_t, ndim=1] labs

    # For fast evaluation of if statements by cython
    
    ny = im.shape[0]
    nx = im.shape[1]
    lab = numpy.zeros((ny, nx), dtype=numpy.int32)
    vals = numpy.array([0, 0, 0, 0])
    labs = numpy.array([1, 0, 0, 0])

    # The upper left pixel is in region zero.
    lastlabel = 0
    equiv = [lastlabel]
    if zero_regions == 0:
        if im[0,0] != 0:
            lastlabel += 1
            equiv.append (lastlabel)
            lab[0,0] = lastlabel
    else:
        lab[0,0] = 0
        
    # Process the rest of the first row of the image.
    y = 0
    if ( zero_regions == 0 ):
        for x in range (1, nx):
            if ( im[y,x] != 0 ):
                if im[y,x] != im[y,x-1]:
                    lastlabel += 1
                    equiv.append (lastlabel)
                lab[y,x] = lastlabel
    else:
        for x in range (1, nx):
            if im[y,x] != im[y,x-1]:
                lastlabel += 1
                equiv.append (lastlabel)
            lab[y,x] = lastlabel

    # Process the first column of the image.
    x = 0
    if ( zero_regions == 0 ):
        for y in range (1, ny):
            if im[y,x] != 0:
                if im[y,x] == im[y-1,x]:
                    lv = lab[y-1,x]
                else:
                    lastlabel += 1
                    equiv.append (lastlabel)
                    lv = lastlabel
                lab[y,x] = lv
    else:
        for y in range (1, ny):
            if im[y,x] == im[y-1,x]:
                lv = lab[y-1,x]
            else:
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            lab[y,x] = lv
                
    # Process the remainder of the image.
    for y in range (1, ny):
        y1 = y - 1
        for x in range (1, nx):
            if con8 == 1: nv = 4
            else:    nv = 2
            x1 = x - 1
            x2 = x + 1
            if ( x2 >= nx - 1 ) and ( con8 == 1 ): lastcol = 1; nv -= 1
            else: lastcol = 0

            # Do not process if this is a 0 pixel (goes to label 0)
            
            if ( im[y,x] == 0 ) and (zero_regions == 0): continue
            
            val = im[y,x]
            # Get the four neighbours' values and labels, taking care not
            # to index off the end of the image.
            vals[0] = im[y, x1]; labs[0] = lab[y, x1]
            vals[1] = im[y1,x]; labs[1] = lab[y1,x]
            if con8 == 1:
                vals[2] = im[y1,x1]; labs[2] = lab[y1,x1]
                if lastcol == 0:
                    vals[3] = im[y1,x2]; labs[3] = lab[y1,x2]
            inreg = False
            for i in range (0, nv):
                if val == vals[i]: inreg = 1
            if inreg == 0:
                # We're in a new region.
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            else:
                # We must be in the same region as a neighbour.
                matches = []
                for i in range (0, nv):
                    if val == vals[i]: matches.append (labs[i])
                matches.sort ()
                lv = int(matches[0])
                for v in matches[1:]:
                    if equiv[v] > lv:
                        equiv[v] = lv
                    elif lv > equiv[v]:
                        equiv[lv] = equiv[v]
            lab[y,x] = lv

    # Tidy up the equivalence table.
    remap = list()
    nc = -1
    for i in range (0, len(equiv)):
        if equiv[i] == i:
            nc += 1
            v = nc
        else:
            v = i
            while equiv[v] != v:
                v = equiv[v]
            v = remap[v]
        remap.append (v)

    # Make a second pass through the image, re-labelling the regions, then
    # return the labelled image.
    for y in range (0, ny):
        for x in range (0, nx):
            lab[y,x] = remap[lab[y,x]]
    return lab

cdef label_regions_int32 (numpy.ndarray[numpy.int32_t, ndim=2] im, int con8,int zero_regions):
    """
    Given a segmented image, return an image with its regions labelled.

    Arguments:
    im -- image to be labelled
    con8 -- if True, consider all 8 nearest neighbours
            if False, consider only 4 nearest neighbours

    Returns: labels in a two-dimensional int32 array            
    """
    cdef:
        int nx, ny, x, y
        int lastcol, inreg
        numpy.ndarray[numpy.int32_t, ndim=2] lab
        numpy.ndarray[numpy.int32_t, ndim=1] vals
        numpy.ndarray[numpy.int32_t, ndim=1] labs

    # For fast evaluation of if statements by cython
    
    ny = im.shape[0]
    nx = im.shape[1]
    lab = numpy.zeros((ny, nx), dtype=numpy.int32)
    vals = numpy.array([0, 0, 0, 0])
    labs = numpy.array([1, 0, 0, 0])

    # The upper left pixel is in region zero.
    lastlabel = 0
    equiv = [lastlabel]
    if zero_regions == 0:
        if im[0,0] != 0:
            lastlabel += 1
            equiv.append (lastlabel)
            lab[0,0] = lastlabel
    else:
        lab[0,0] = 0
        
    # Process the rest of the first row of the image.
    y = 0
    if ( zero_regions == 0 ):
        for x in range (1, nx):
            if ( im[y,x] != 0 ):
                if im[y,x] != im[y,x-1]:
                    lastlabel += 1
                    equiv.append (lastlabel)
                lab[y,x] = lastlabel
    else:
        for x in range (1, nx):
            if im[y,x] != im[y,x-1]:
                lastlabel += 1
                equiv.append (lastlabel)
            lab[y,x] = lastlabel

    # Process the first column of the image.
    x = 0
    if ( zero_regions == 0 ):
        for y in range (1, ny):
            if im[y,x] != 0:
                if im[y,x] == im[y-1,x]:
                    lv = lab[y-1,x]
                else:
                    lastlabel += 1
                    equiv.append (lastlabel)
                    lv = lastlabel
                lab[y,x] = lv
    else:
        for y in range (1, ny):
            if im[y,x] == im[y-1,x]:
                lv = lab[y-1,x]
            else:
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            lab[y,x] = lv
                
    # Process the remainder of the image.
    for y in range (1, ny):
        y1 = y - 1
        for x in range (1, nx):
            if con8 == 1: nv = 4
            else:    nv = 2
            x1 = x - 1
            x2 = x + 1
            if ( x2 >= nx - 1 ) and ( con8 == 1 ): lastcol = 1; nv -= 1
            else: lastcol = 0

            # Do not process if this is a 0 pixel (goes to label 0)
            
            if ( im[y,x] == 0 ) and (zero_regions == 0): continue
            
            val = im[y,x]
            # Get the four neighbours' values and labels, taking care not
            # to index off the end of the image.
            vals[0] = im[y, x1]; labs[0] = lab[y, x1]
            vals[1] = im[y1,x]; labs[1] = lab[y1,x]
            if con8 == 1:
                vals[2] = im[y1,x1]; labs[2] = lab[y1,x1]
                if lastcol == 0:
                    vals[3] = im[y1,x2]; labs[3] = lab[y1,x2]
            inreg = False
            for i in range (0, nv):
                if val == vals[i]: inreg = 1
            if inreg == 0:
                # We're in a new region.
                lastlabel += 1
                equiv.append (lastlabel)
                lv = lastlabel
            else:
                # We must be in the same region as a neighbour.
                matches = []
                for i in range (0, nv):
                    if val == vals[i]: matches.append (labs[i])
                matches.sort ()
                lv = int(matches[0])
                for v in matches[1:]:
                    if equiv[v] > lv:
                        equiv[v] = lv
                    elif lv > equiv[v]:
                        equiv[lv] = equiv[v]
            lab[y,x] = lv

    # Tidy up the equivalence table.
    remap = list()
    nc = -1
    for i in range (0, len(equiv)):
        if equiv[i] == i:
            nc += 1
            v = nc
        else:
            v = i
            while equiv[v] != v:
                v = equiv[v]
            v = remap[v]
        remap.append (v)

    # Make a second pass through the image, re-labelling the regions, then
    # return the labelled image.
    for y in range (0, ny):
        for x in range (0, nx):
            lab[y,x] = remap[lab[y,x]]
    return lab

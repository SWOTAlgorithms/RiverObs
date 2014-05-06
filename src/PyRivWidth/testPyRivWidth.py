import sys
from os.path import join
import numpy as N
from numpy.ma import masked_array
from pylab import *

## rivwidth_lib = '../'
## sys.path.append(rivwidth_lib)
from RivWidthHeight import RivWidthHeight
from label_regions import label_regions

input_data_dir = './'

test_chan_mask = N.array([
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
    [1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1],
    [1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,1],
    [1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
    [1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1],
    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1],
    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1],
    [0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ], dtype=N.uint8)

def test_file_init(fname='mask.nc',braidflag=0):
    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    print rw.rundata

    return rw

def test_make_channel_mask(fname='mask.nc',braidflag=0):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    figure(figsize=(12,8))
    imshow(chan_maskr,cmap=cm.gray)    
    title(r'chan_maskr, braidflag: %d'%braidflag)

    figure(figsize=(12,8))
    imshow(chan_mask,cmap=cm.gray)    
    title(r'chan_mask, braidflag: %d'%braidflag)
    
def test_make_river_mask(fname='mask.nc',braidflag=0):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=True)
    
    figure(figsize=(12,8))
##    imshow(riv_mask[::subsample,::subsample],cmap=cm.gray)
    imshow(riv_mask,cmap=cm.gray,interpolation='nearest')    
    title(r'riv_mask, braidflag: %d'%braidflag)

    return riv_mask

def test_center_line_calc(fname='mask.nc',braidflag=0,deriv_threshold=0.85,use_thin=True):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=True)
    
    center_line_mask = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold,use_thin=use_thin,output_file=True)
    
    figure(figsize=(12,8))
    imshow(center_line_mask,cmap=cm.gray,interpolation='nearest')    
    title(r'center_line, braidflag: %d'%braidflag)

    return center_line_mask

def test_dynamic_vector_v3(fname='mask.nc',braidflag=0,deriv_threshold=0.85,use_riv_mask=False,use_thin=True):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=True)

    if use_riv_mask:
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(riv_mask, srow, scol,  erow, ecol)
    else:
        center_line_mask = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold,output_file=True,use_thin=use_thin)
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(center_line_mask, srow, scol,  erow, ecol)
        
##     figure(figsize=(12,8))
## ##    imshow(riv_mask[::subsample,::subsample],cmap=cm.gray)
##     imshow(riv_mask,cmap=cm.gray,interpolation='nearest')    
##     title(r'riv_mask, braidflag: %d'%braidflag)

    return center_line

def test_get_center_line_parameters(fname='mask.nc',braidflag=1,deriv_threshold=0.85,use_riv_mask=False,
                                    use_thin=True,ds=30.):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=True)

    if use_riv_mask:
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(riv_mask, srow, scol,  erow, ecol)
    else:
        center_line_mask = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold,output_file=True,use_thin=use_thin)
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(center_line_mask, srow, scol,  erow, ecol)

    print "get_center_line_parameters"
    rw.get_center_line_parameters(center_line,  ds)

    return rw

def test_get_s_c_coordinates(fname='mask.nc',braidflag=1,deriv_threshold=0.85,use_riv_mask=False,
                             use_thin=True,ds=30.,xshift=10,yshift=10):

    rundata = {}
    rundata['braidnarrowflag'] = braidflag

    rw = RivWidthHeight(rundata,fname)
    chan_maskr, chan_mask = rw.make_channel_mask(output_file=True)

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=True)

    if use_riv_mask:
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(riv_mask, srow, scol,  erow, ecol)
    else:
        center_line_mask = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold,output_file=True,use_thin=use_thin)
        srow = 210
        scol = 1
        erow = 106
        ecol = 310
        center_line = rw.dynamic_vector_v3(center_line_mask, srow, scol,  erow, ecol)

    print "get_center_line_parameters"
    rw.get_center_line_parameters(center_line,  ds)

    # Make a test line

    print 'Getting S,C'
    
    nxy = rw.interpolated_center_line.shape[0]
    xy = N.zeros((nxy,2),dtype=N.float64)

    xy[:,0] = rw.interpolated_center_line[:,0] + xshift
    xy[:,1] = rw.interpolated_center_line[:,1] + yshift

    S, C, index = rw.get_s_c_coordinates(xy)

    return S,C, index

def simple_test_make_river_mask(braidflag=0):

    rundata = {}
    rundata['height'] = test_chan_mask.shape[0]
    rundata['width'] = test_chan_mask.shape[1]
    rundata['braidnarrowflag'] = braidflag
    
    rw = RivWidthHeight(rundata=rundata)

    print "making channel masks"
    chan_maskr, chan_mask = test_chan_mask, test_chan_mask

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=None)
    
    figure(figsize=(12,8))
##    imshow(riv_mask[::subsample,::subsample],cmap=cm.gray)
    imshow(riv_mask,cmap=cm.gray)    
    title(r'riv_mask, braidflag: %d'%braidflag)

    return riv_mask

def simple_test_center_line_calc(braidflag=0,deriv_threshold=0.85):

    rundata = {}
    rundata['height'] = test_chan_mask.shape[0]
    rundata['width'] = test_chan_mask.shape[1]
    rundata['braidnarrowflag'] = braidflag
    
    rw = RivWidthHeight(rundata=rundata)

    print "making channel masks"
    chan_maskr, chan_mask = test_chan_mask, test_chan_mask

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=None)

    center_line = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold)
    
    figure(figsize=(12,8))
##    imshow(riv_mask[::subsample,::subsample],cmap=cm.gray)
    imshow(center_line,cmap=cm.gray)    
    title(r'center_line, braidflag: %d'%braidflag)

    return center_line

def simple_test_dynamic_vector_v3(braidflag=0,deriv_threshold=0.85):

    rundata = {}
    rundata['height'] = test_chan_mask.shape[0]
    rundata['width'] = test_chan_mask.shape[1]
    rundata['braidnarrowflag'] = braidflag

    rundata['spx'] = 1
    rundata['spy'] = 8
    rundata['epx'] = test_chan_mask.shape[1]-2
    rundata['epy'] = 9

    print 'starting point: (%d,%d)'%(1,8)
    print 'end point: (%d,%d)'%(test_chan_mask.shape[1]-2,8)
    
    rw = RivWidthHeight(rundata=rundata)

    print "making channel masks"
    chan_maskr, chan_mask = test_chan_mask, test_chan_mask

    print "making river mask"
    riv_mask = rw.make_river_mask(chan_maskr,chan_mask,max_iter = 5,output_file=None)

    center_line = rw.center_line_calc(riv_mask,deriv_threshold=deriv_threshold)

    center_line_array,center_line_image = rw.dynamic_vector_v3(center_line)

    figure(figsize=(12,8))
##    imshow(riv_mask[::subsample,::subsample],cmap=cm.gray)
    imshow(center_line_image,cmap=cm.gray)
    title(r'center_line')

    return riv_mask

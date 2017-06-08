"""
Take GWDLR files and turn them into shape files. These programs need to be run in
a grass shell.
"""

from __future__ import absolute_import, division, print_function

import os
from os.path import join
from subprocess import call
import shlex
from .GWDLR import GWDLR

class GWDLR2shape:
    """Take GWDLR files and turn them into shape files. These programs need to be run in
    a grass shell."""

    def __init__(self,gwdlr_data_dir,output_dir):
        """Initialize with input and output data specifications.

        gwdlr_data_dir: contains .bin and .ctl data
        output_dir: write shapefiles and mask files to this directory"""

        # Check that we are working in a grass shell (for the moment)

        if not 'GISBASE' in os.environ:
            raise Exception('Need to be working inside a grass shell')

        self.gwdlr_data_dir = gwdlr_data_dir
        self.output_dir = output_dir

    def process_tile(self,rootname,min_width,create_location=True):
        """Process a tile based on selecting all widths > min_width."""

        # Make the mask file

        gwdlr = GWDLR(rootname,data_dir=self.gwdlr_data_dir)
        gwdlr.to_mask(min_width,overwrite=True,thin=True)

        root = rootname.split('_')[0]
        mask_file = join(self.output_dir,root+'_mask_width%d.tif'%min_width)
        gwdlr.to_gdal(mask_file)

        # Start a new location with the mask data

        if create_location:
            command = 'g.proj -c georef=%(mask_file)s location=%(root)s --verbose'%locals()
            self.exec_command(command)

        # Switch to the new location

        command= 'g.mapset mapset=PERMANENT location=%(root)s'%locals()
        self.exec_command(command)

        # Read the mask data into the location

        command = 'r.in.gdal input=%(mask_file)s output=center_line_mask_width%(min_width)d'%locals()
        self.exec_command(command)

        # Thin the mask to make sure that r.to.vect does not crash

        command = 'r.thin input=center_line_mask_width%(min_width)d output=center_line_mask_thin_width%(min_width)d --verbose'%locals()
        self.exec_command(command)

        # Make a vector from the thinned data

        command = 'r.to.vect input=center_line_mask_thin_width%(min_width)d output=center_line_width%(min_width)d feature=line --verbose'%locals()
        self.exec_command(command)

        # Export to a shapefile

        dsn = join(self.output_dir,root+'_center_lines_width%d'%min_width)
        command = 'v.out.ogr input=center_line_width%(min_width)d dsn=%(dsn)s'%locals()
        self.exec_command(command)

    def exec_command(self,command):
        """Execute a grass command and catch errors."""

        print('Executing command: %s'%command)
        args = shlex.split(command)
        status = call(args)

        if status:
            print('Command: "%s" exited with status: %s'%(command,status))

        return status

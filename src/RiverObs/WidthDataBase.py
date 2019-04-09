"""
Query a pandas HDFStore width data base for rivers.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import pandas


class WidthDataBase:
    """Query a pandas HDFStore width data base for rivers.

    Parameters
    ----------

    db_file : str
        pandas HDFStore file containing river and reach data frames.
    river_df_name : str
        Name of the pandas DataFrame containing point
        width information (default: 'river').
    reach_df_name : str
        Name of the pandas DataFrame containing
        reach statistical information (default: 'reach').
    mode: str
        how to open the file (default: 'r' read only).
    reach_index_kwd: str
        column name uniquely identifying the reach.
        (default: 'reach_index', but could be the reach name, etc.)
    """

    def __init__(self,
                 db_file,
                 river_df_name='river',
                 reach_df_name='reach',
                 mode='r',
                 reach_index_kwd='reach_index'):
        """
        """

        self.h5 = pandas.HDFStore(db_file, mode=mode)

        self.river_df = self.h5[river_df_name]
        self.reach_df = self.h5[reach_df_name]

    def get_river(self,
                  reach_index,
                  columns=None,
                  asarray=False,
                  transpose=False,
                  lat_kwd='lat',
                  lon_kwd='long',
                  bounding_box=None,
                  clip_buffer=0):
        """Return selected information for the reach index by reach_index.

        Parameters
        -----------

        reach_index : int
            index identifying the reach.
        columns : list
            if None, all of the iformation associated with the reach is
            returned. Otherwise, pass either a column name or list of
            column names;e.g., columns='width' orcolumns=['long','lat','width'].
        asarray : bool
            if True, returns a numpy ndarray rather than a pandas DataFrame.
        transpose: bool
            if True and asarray=True, return the transpose array.
            This is useful to unpack arrays easily.
        lat_kwd : str
            latitude column name in the data base (default 'lat')
        lon_kwd : str
            latitude column name in the data base (default 'long')
        bounding_box: tuple
            (lonmin, latmin, lonmax, latmax). If not None, the
            only lon/lat in the bounding box + clip_buffer are returned.
        """

        if bounding_box != None:
            lon, lat, inbbox = self.get_lon_lat(
                reach_index,
                lat_kwd=lat_kwd,
                lon_kwd=lon_kwd,
                bounding_box=bounding_box,
                clip_buffer=clip_buffer)
        else:
            inbbox = None

        # Select the DataFrame for this river

        df = self.river_df[self.river_df['reach_index'] == reach_index]

        # Select the desired columns

        if columns != None:
            df = df[columns]

        # If a bounding box has been specified, extract the appropriate records

        if type(inbbox) != type(None):
            df = df.iloc[inbbox]

        # Return the desired columns in the desired format

        if not asarray:
            return df
        else:
            if transpose:
                return np.asarray(df).T
            else:
                return np.asarray(df)

    def get_lon_lat(self,
                    reach_index,
                    lat_kwd='lat',
                    lon_kwd='long',
                    bounding_box=None,
                    clip_buffer=0):
        """Return the latitude and longitude associated with a reach and
        a clip box.

        Parameters
        ----------

        lat_kwd : str
            latitude column name in the data base (default 'lat')
        lon_kwd : str
            latitude column name in the data base (default 'long')
        bounding_box: tuple
            (lonmin, latmin, lonmax, latmax). If not None, the
            only lon/lat in the bounding box + clip_buffer are returned.

        Returns
        -------

        Returns lon, lat numpy arrays. If bounding_box != None, also returns
        an index array for the good data.
        """

        lon, lat = self.get_river(
            reach_index,
            columns=[lon_kwd, lat_kwd],
            asarray=True,
            transpose=True)

        if bounding_box != None:
            inbbox = ((lon >= bounding_box[0] - clip_buffer) &
                      (lat >= bounding_box[1] - clip_buffer) &
                      (lon <= bounding_box[2] + clip_buffer) &
                      (lat <= bounding_box[3] + clip_buffer))
            lon = lon[inbbox]
            lat = lat[inbbox]
            return lon, lat, inbbox

        return lon, lat

    def get_xy(self,
               reach_index,
               proj,
               lat_kwd='lat',
               lon_kwd='long',
               bounding_box=None,
               clip_buffer=0):
        """Given a projection function (e.g., from pyproj.Proj ) return x,y.

        Parameters
        ----------

        proj : function
            x,y = proj(lon,lat)
        lat_kwd: str
            latitude column name in the data base (default 'lat')
        lon_kwd: str
            latitude column name in the data base (default 'long')

        Returns
        -------

        Returns x,y numpy arrays.
        """

        if bounding_box == None:
            lon, lat = self.get_lon_lat(
                reach_index,
                lat_kwd=lat_kwd,
                lon_kwd=lon_kwd,
                bounding_box=bounding_box,
                clip_buffer=clip_buffer)
        else:
            lon, lat, inbbox = self.get_lon_lat(
                reach_index,
                lat_kwd=lat_kwd,
                lon_kwd=lon_kwd,
                bounding_box=bounding_box,
                clip_buffer=clip_buffer)

        return proj(lon, lat)

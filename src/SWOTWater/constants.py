"""
Module for swot constants
"""

PIXC_CLASSES = {
    'invalid': -1,
    'land': 1,
    'land_near_water': 2,
    'water_near_land': 3,
    'open_water': 4,
    'dark_water': 5,
    'low_coh_water_near_land': 6,
    'low_coh_water': 7,
    'land_near_dark_water': 22,# legacy/depreciated
    'dark_water_near_land': 23,# legacy/depreciated
    'dark_water_legacy': 24 # legacy/depreciated
    }

AGG_CLASSES = {
    'interior_water_klasses':[
        PIXC_CLASSES['open_water'],
        #PIXC_CLASSES['low_coh_water'],
        ],
    'water_edge_klasses':[
        PIXC_CLASSES['water_near_land'],
        ],
    'land_edge_klasses':[
        PIXC_CLASSES['land_near_water'],
        ],
    'dark_water_klasses':[
        PIXC_CLASSES['dark_water'],
        PIXC_CLASSES['land_near_dark_water'],
        PIXC_CLASSES['dark_water_near_land'],
        PIXC_CLASSES['dark_water_legacy'],
        ]
    }

GDEM_PIXC_CLASSES = {
    'open_water': 4, 'dark_water': 24,
    'open_water_lake': 34, 'dark_water_lake': 54}

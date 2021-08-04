"""
Module for swot constants
"""

PIXC_CLASSES = {
    'land': 1, 'land_near_water': 2, 'water_near_land': 3, 'open_water': 4,
    'land_near_dark_water': 22, 'dark_water_edge': 23, 'dark_water': 24,
    'dark_water_klasses':[5,23,24], 'low_coherence_klasses':[6,7]}  # includes legacy values

GDEM_PIXC_CLASSES = {
    'open_water': 4, 'dark_water': 24,
    'open_water_lake': 34, 'dark_water_lake': 54}

'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import textwrap
import numpy as np
from collections import OrderedDict as odict

from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

class L2PIXCVector(Product):
    UID = "l2_hr_pixcvector"
    ATTRIBUTES = []
    VARIABLES = []

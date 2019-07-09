'''
Copyright (c) 2019-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
from SWOTRiver.products.rivertile import L2HRRiverTile

class L2HRRiverSP(L2HRRiverTile):
    """
    Class for River SP data product
    """
    UID = "l2_hr_riversp"
    @classmethod
    def from_rivertiles(cls, rivertiles):
        """
        Builds a River SP product from a list of rivetile data products
        """
        klass = cls()
        # It does somthing like this
        for rivertile in rivertiles:
           klass += rivertile
        return klass

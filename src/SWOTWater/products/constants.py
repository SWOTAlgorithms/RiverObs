'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore

'''
FILL_VALUES = {
    'double': 9.9692099683868690e+36,
    'single': 9.96921e+36,
    'int64': 2**63-1, 'uint64': 2**64-1,
    'int32': 2**31-1, 'uint32': 2**32-1,
    'int16': 2**15-1, 'uint16': 2**16-1,
    'int8': 2**7-1, 'uint8': 2**8-1,
    'char': '*',
}
FILL_VALUES['c16'] = FILL_VALUES['double'] + 1j*FILL_VALUES['double']
FILL_VALUES['c8'] = FILL_VALUES['single'] + 1j*FILL_VALUES['single']
FILL_VALUES['f8'] = FILL_VALUES['double']
FILL_VALUES['f4'] = FILL_VALUES['single']
FILL_VALUES['i8'] = FILL_VALUES['int64']
FILL_VALUES['i4'] = FILL_VALUES['int32']
FILL_VALUES['i2'] = FILL_VALUES['int16']
FILL_VALUES['i1'] = FILL_VALUES['int8']
FILL_VALUES['u8'] = FILL_VALUES['uint64']
FILL_VALUES['u4'] = FILL_VALUES['uint32']
FILL_VALUES['u2'] = FILL_VALUES['uint16']
FILL_VALUES['u1'] = FILL_VALUES['uint8']
FILL_VALUES['S1'] = FILL_VALUES['char']

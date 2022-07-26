#!/usr/bin/env python
'''
Copyright(c) 2017-, California Institute of Technology("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Brent Williams
'''
import argparse
import pdb

import SWOTRiver.analysis.tabley


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('table_file1', type=str)
    parser.add_argument('table_file2', type=str)
    parser.add_argument('--output_file', '-o', type=str, default=None)
    args = parser.parse_args()
    print(args)

    # read in each file
    table1 = SWOTRiver.analysis.tabley.Table.from_file(args.table_file1)
    table2 = SWOTRiver.analysis.tabley.Table.from_file(args.table_file2)
    n_missing = 0
    #print(table1)
    #print(table2)
    # match-up lines
    if not (table1.headers == table2.headers):
        print("tables are not comparable, headers dont match")
        return
    #print(table1.headers)
    print('Total reaches in', args.table_file1, 'is', len(table1.data))
    print('Total reaches in', args.table_file2, 'is', len(table2.data))

    tile_ind = -1
    reach_ind = -1
    for k, hdr in enumerate(table1.headers):
        if 'tile' in hdr:
            tile_ind = k
        if 'reach' == hdr:
            reach_ind = k
    outdata = []
    for line1 in table1.data:
        #print(line1)
        tile1 = None
        reach1 = None
        matchlinenum = -1
        if tile_ind > 0:
            tile1 = line1[tile_ind]
        if reach_ind > 0:
            reach1 = line1[reach_ind]
        for k2,line2 in enumerate(table2.data):
            tile2 = None
            reach2 = None
            if tile_ind > 0:
                tile2 = line2[tile_ind]
            if reach_ind > 0:
                reach2 = line2[reach_ind]
            if (tile1==tile2) and (reach1==reach2) and (tile1 is not None):
                matchlinenum = k2
                #print (tile1, tile2, reach1, reach2)
        if matchlinenum > 0:
            matchline = table2.data[matchlinenum]
            outline = []
            for k, item in enumerate(line1):
                if isinstance(item, str):
                    outitem = item
                else:
                    # difference the matched line
                    outitem = item - matchline[k]
                outline.append(outitem)
            outdata.append(outline)
            #print("##",matchline)
            #print("***",outline)
        else:
            # no matching line in table2 for this table1 line
            print('No matching line in', args.table_file2,
                  'for', tile1, reach1, 'in', args.table_file1)
            n_missing += 1
    # create the difference table
    # setup preamble
    print('total missing reaches', n_missing)
    preamble = 'Difference table: %s - %s'%(args.table_file1, args.table_file2)
    # TODO: setup passfail dictionary for difference
    SWOTRiver.analysis.tabley.print_table(
        outdata, headers=table1.headers, style=None,
        precision=table1.precision, width=table1.fixed_width,
        passfail={}, fname=args.output_file, preamble=preamble)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass


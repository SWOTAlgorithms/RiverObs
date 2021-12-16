#!/usr/bin/env python
import argparse
import inspect
import io
import logging
import sys
import warnings

import numpy as np


def print_table(*args, **kwargs):
    """Create a Table() from arguments and print it."""
    table = Table(*args, **kwargs)
    print('table.fname:',table.fname)
    if table.fname is not None:
        # write to file if fname is set
        #f = open(table.fname,'w')
        #f.write(table.__str__())
        #f.close()
        table.to_file()
    else:
        # print to terminal
        print(table)


class Table():
    """Turn 2D data into a pretty table, suitable for printing or markdown.

    Arguments:
        data: a 2D[row, column] structure, or a dict of {header: column} pairs
        headers: a 1D[column] structure of table column headers.
        style: either None or 'pipe'
        precision: floating point decimal precision, defaults
            {:[width].[precision]g}
        width: a minimum width for each column

    This class will fit a minimum width for each column >= width, and can be
    printed.
    """
    def __init__(
            self, data, headers=None, style=None, precision=4, width=0, passfail={},
            fname=None, preamble=None):
        if isinstance(data, dict):
            headers = list(data.keys())
            # data = np.column_stack(data.values())
            tmp = data
            data = None
            for column in tmp.values():
                if data is None:
                    data = [[] for item in column]
                for i_row, row_item in enumerate(column):
                    data[i_row].append(row_item)
        self.data = data
        self.headers = headers
        self.precision = precision
        self.fixed_width = width
        self.widths = None
        self.formats = None
        self.passfail = passfail
        self.fname = fname
        self.preamble = preamble
        self.fit_widths()
        if style == 'pipe':
            self.border = '|'
            self.separator = '|'
        else:
            self.border = ''
            self.separator = ' '

    def fit_widths(self):
        """Fit minimum width for each column using self.data and self.width.

        This is called automatically in __init__()."""
        n_columns = 0
        if len(self.data)>0:
            n_columns = len(self.data[0])
        self.widths = np.zeros(n_columns, dtype='i4')
        self.formats = ['']*n_columns
        for i in range(n_columns):
            if self.headers is None:
                width_header = 0
            else:
                width_header = len(self.headers[i])
            if isinstance(self.data[0][i], str):
                width_data = 0
                for j in range(len(self.data)):
                    width_data = max(width_data, len(self.data[j][i]))
            else:
                width_data = self.precision + 6
            self.widths[i] = max(width_header, width_data)
            if self.widths[i] < self.fixed_width:
                self.widths[i] = self.fixed_width
            if isinstance(self.data[0][i], str):
                self.formats[i] = '{{:>{0}}}'.format(self.widths[i])
            else:
                # check if this data has a pass/fail condition to set colors
                self.formats[i] = '{{:#{0}.{1}g}}'.format(
                            self.widths[i], self.precision)
                #if self.headers[i] in self.passfail:
                #    passes, fails = self.passfail[self.headers[i]]
                #    if np.abs(self.data[0][i]) > fails:
                #        print("\n\nRED\n\n")
                #        # print in red
                #        self.formats[i] = '\033[91m {{:#{0}.{1}g}}\033[00m'.format#(
                #            self.widths[i], self.precision)
                #    elif np.abs(self.data[0][i]) < passes:
                #        print("\n\nGREEN\n\n")
                #        # print in green
                #        self.formats[i] = '\033[92m{{:#{0}.{1}g}}\033[00m'.format(
                #            self.widths[i], self.precision)
                #    print("\n\nBLACK\n\n")

            logging.debug('%d %d %s', i, self.widths[i], self.formats[i])

    def _header_line(self):
        string = ''
        if self.headers is not None:
            string = self.border
            for header, width in zip(self.headers, self.widths):
                string = string + '{0:>{1}}'.format(header, width)
                string = string + self.separator
        return string[:-1] + self.border

    def _horizontal_rule(self):
        string = self.border
        if self.separator == '|':
            aligner = ':'
        else:
            aligner = '-'
        for width in self.widths:
            string = string + '-'*(width-1) + aligner + self.separator
        return string[:-1] + self.border

    def _data_lines(self):
        string = ''
        for row in self.data:
            string = string + self.border
            for i, (this_format, item) in enumerate(zip(self.formats, row)):
                my_format = this_format
                if self.headers[i] in self.passfail:
                    passes, fails = self.passfail[self.headers[i]]
                    if not(isinstance(item, str)):
                        if np.abs(item) > fails:
                            my_format = '\033[91m'+this_format+'\033[00m'
                        elif np.abs(item) < passes:
                            my_format = '\033[92m'+this_format+'\033[00m'
                string = string + my_format.format(item) + self.separator
            string = string[:-1] + self.border + '\n'
        return string

    def __str__(self):
        string = ''
        if self.preamble is not None:
            string = self.preamble + '\n'
        if self.headers is not None:
            string = string + self._header_line() + '\n'
        string = string + self._horizontal_rule() + '\n'
        string = string + self._data_lines()
        if self.headers is None:
            string = string + self._horizontal_rule() + '\n'
        return string

    def to_file(self, fname=None):
        if fname is None:
            if self.fname is None:
                print('No filename given so not writting to file')
                return
            else:
                fname = self.fname
        f = open(fname,'w')
        f.write(self.__str__())
        f.close()

    @classmethod
    def from_file(cls, fname):
        # create a Table object and populate the contents from the file
        # read the table
        f = open(fname,'r')
        lines = f.readlines()
        f.close()
        
        # decode the lines
        # search for the dashed-line
        dash_line_num = 0
        for k, line in enumerate(lines):
            if '---' in line:
                dash_line_num = k
        dash_line = lines[dash_line_num]
        # header is the line before the dashed line
        header = lines[dash_line_num-1]
        # preamble is the line before the header
        if dash_line_num<2:
            preamble = None
        else:
            preamble = lines[dash_line_num-2]
        #print('dash_line_num', dash_line_num)
        #print(preamble)
        #print(header)
        #print(dash_line)
        # set some defaults
        style = None
        precision=4
        width=0
        passfail = {} # don't do pass-fail when reading, tha ANSI escape codes should already be there
        separator = None
        # handle pipe style
        if len(dash_line.split('|')) > 1:
            style = 'pipe'
            separator = '|'
        dashes = dash_line.split(separator)
        #print(dashes)
        # use the dash line to decode the header line, which may have spaces
        l_accum = 0
        headers = []
        for dash_word in dashes:
            l = len(dash_word)
            headers.append(header[l_accum:l_accum+l].lstrip().rstrip())
            #print("dash_word[%d:%d]"%(l_accum, l))
            l_accum = l_accum+l+1
        #print(headers)
        data = {}
        # init the data dictionary with empty lists 
        for key in headers:
            data[key] = []
        for line in lines[dash_line_num+1:]:
            parts0 = line.split(separator)
            # merge back some escape sequences
            parts = []
            pass_fail_line = []
            kp = 0
            for k, h in enumerate(headers):
                part = parts0[kp]
                pf = 'none'
                dat_str = part.strip()
                if (dat_str.startswith('(') and dat_str.endswith(')')) and k>0:
                    # this is a unit separated by space that should
                    # go with to the previous item
                    #print('part:',parts[-1],dat_str)
                    parts[-1] = parts[-1] + ' ' + dat_str
                    kp=kp+1
                    part = parts0[kp]
                # handle ANSI escape color codes
                if ('\033[91m' in part or '\033[92m' in part) and not (
                        '\033[00m' in part):
                    part = part + parts0[kp+1]
                    if '\033[91m' in part:
                        pf = 'pass'
                    else:
                        pf = 'fail'
                    kp=kp+1
                    # strip out the ANSI escape codes now
                    #print('part:', part, part[5:-5])
                    part = part[5:-5]
                dat_str = part.strip()
                #print(dat_str, k)
                dat = dat_str
                # if it has a decimal point try to cast to float
                # also detect the max precision
                if '.' in dat_str:
                    this_precision = -1
                    try:
                        dat = float(dat_str)
                        this_precision = len(dat_str.split('e')[0].strip('-'))
                    except ValueError:
                        pass
                    if precision < this_precision:
                        # set precision to max precision found of any float
                        precision = this_precision
                parts.append(dat)
                pass_fail_line.append(pf)
                kp=kp+1
            #print("len(parts), len(headers):",len(parts), len(headers))
            #print(parts)
            for key, part in zip(headers, parts):
                data[key].append(part)
        # TODO: detect 'width', and estimate pass-fail criteria based on data...
        # create instance
        table = cls(data, headers, style, precision, width, passfail,
            fname=None, preamble=preamble)
        return table

def main():
    log_format = '%(levelname)s:%(module)s.%(funcName)s:%(message)s'
    logging.basicConfig(format=log_format, level=logging.INFO)
    parser = argparse.ArgumentParser()
    signature = inspect.signature(Table)
    parser.add_argument('-s', '--split', default=' ')
    for parameter in signature.parameters.values():
        if parameter.name in ('data', 'headers'):
            continue
        parser.add_argument(
            '--{}'.format(parameter.name), type=type(parameter.default),
            default=parameter.default)
    parser.add_argument('filename', nargs='?')
    args = parser.parse_args()

    table_args = {}
    for arg, value in vars(args).items():
        if arg in ('split', 'filename'):
            continue
        table_args[arg] = value

    if args.filename is None:
        lines = sys.stdin.read()
    else:
        with open(args.filename) as handle:
            lines = handle.read()
    lines = lines.strip('\n')
    tables = lines.split('\n\n')
    for table in tables:
        logging.debug('\n%s', table)
        this_table = np.genfromtxt(
            table.split('\n'), dtype=None, encoding=None)
        logging.debug('\n%s', this_table)
        print_table(this_table, **table_args)


if __name__ == '__main__':
    main()

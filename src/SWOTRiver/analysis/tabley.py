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
            self, data, headers=None, style=None, precision=4, width=0):
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
                self.formats[i] = '{{:#{0}.{1}g}}'.format(
                    self.widths[i], self.precision)
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
            for this_format, item in zip(self.formats, row):
                string = string + this_format.format(item) + self.separator
            string = string[:-1] + self.border + '\n'
        return string
        return string

    def __str__(self):
        string = ''
        if self.headers is not None:
            string = self._header_line() + '\n'
        string = string + self._horizontal_rule() + '\n'
        string = string + self._data_lines()
        if self.headers is None:
            string = string + self._horizontal_rule() + '\n'
        return string


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

#! /usr/bin/env python
"""
A class for parsing, storing, and writing files of the format:

    keyword (units) = value
"""

from __future__ import absolute_import, division, print_function

import re  #regular expressions for parsing
#import string #string library
import types
import locale


class RDF:
    """A class for parsing, storing, and writing files of the format:

    keyword (units) = value

    Parameters
    ----------

    separator : char
        Separates keyword (unit) from value (default '='(
    continuation : str
        Line continuation.
    comment: char
        Anything after this is ignored in a line (default '!')
    eof : str
        For RDF embedded in a binary file, this starts the binary content.
    """

    #set the RDF comment and continuation characters
    # rdfComm = "!"
    # rdfCont = "\\\\"
    # rdfSep = "="

    def __init__(self,
                 separator="=",
                 continuation="\\\\",
                 comment="!",
                 eof="eof"):
        #declare the return dictionaries
        self.value = {}  #the value dictionary
        self.units = {}  #the units dictionary
        self.key_list = []
        self.rdfSep = separator
        self.rdfCont = continuation
        self.rdfComm = comment
        self.eof = eof

    def __getitem__(self, key):
        return self.value[key]

    def __setitem__(self, key, value):
        if type(value) != bytes: value = repr(value)
        key = str.lower(key)
        self.value[key] = value
        self.units[key] = '-'
        if key not in self.key_list:
            self.key_list.append(key)

    def __delitem__(self, key):
        key = str.lower(key)
        del self.value[key]
        del self.units[key]
        self.key_list.remove(key)

    def __len__(self):
        return len(self.value)

    def keys(self):
        """To comply with dictionary interface."""
        return list(self.value.keys())

    def values(self):
        """To comply with dictionary interface."""
        return list(self.value.values())

    def set(self, key, value, units='-'):
        key = str.lower(key)
        if type(value) != bytes: value = repr(value)
        self.value[key] = value
        self.units[key] = units
        if key not in self.key_list:
            self.key_list.append(key)

    def stripComments(self, line):
        """Strip comments from a line"""
        m = re.match("(.*)" + self.rdfComm + "(.*)", line)
        while m:
            line = m.group(1)
            c = str.strip(m.group(2))
            m = re.match("(.*)" + self.rdfComm + "(.*)", line)
        return line

    def rdfParse(self, files):

        #make sure a list was passed
        if type(files) != list:
            files = [files]

        #Loop over all RDF files

        for file in files:

            #open the file
            try:
                fin = open(file, "r")
            except:
                raise Exception("Cannot open file: %s" % file)

            #get a line and parse it
            while 1:

                #read the line
                line = fin.readline()
                if not line: break

                # strip comments

                line = self.stripComments(line)

                #while continuation line, read next line and append
                m = re.match("(.*)" + self.rdfCont, line)  #regexp match
                if m: line = m.group(1)
                while (m != None):
                    nextLine = self.stripComments(fin.readline())
                    m = re.match("(.*)" + self.rdfCont, nextLine)
                    if m: line = line + m.group(1)
                    else: line = line + nextLine

                #extract keyword and value
                m = re.match("([^=]+)" + self.rdfSep + "(.*)", line)
                if m:  #an RDF line

                    k = m.group(1)
                    v = m.group(2)
                    v = str.strip(v)

                    #units key separation
                    m = re.match("(.*)\((.*)\).*", k)
                    if m:
                        u = str.strip(m.group(2))
                        if u == "": u = "-"

                        k = str.lower(str.strip(m.group(1)))
                        self.units[k] = u
                    else:
                        k = str.lower(str.strip(k))
                        self.units[k] = "-"

                    self.value[k] = v
                    if k not in self.key_list: self.key_list.append(k)

                    # if the keyword eof is set, the end of the header of
                    #  an RDF/binary file has been reached. Exit parsing loop

                    if re.match(self.eof, k): break

            #close RDF file
            fin.close()
        return self

    ## def rdfParseString(self,rdfString):
    ##     """Parse a string containing RDF lines separated by \\n."""

    ##     for line in rdfStr.split('\n'):

    ##         # strip comments

    ##         line = self.stripComments(line)

    ##         # while continuation line, read next line and append
    ##         m = re.match("(.*)"+self.rdfCont,line) #regexp match
    ##         if m: line = m.group(1)
    ##         while ( m != None ):
    ##             nextLine = self.stripComments(fin.readline())
    ##             m = re.match("(.*)"+self.rdfCont,nextLine)
    ##             if m: line = line + m.group(1)
    ##             else: line = line + nextLine

    ##         # extract keyword and value
    ##         m = re.match("([^=]+)"+self.rdfSep+"(.*)",line)
    ##         if m: #an RDF line

    ##             k = m.group(1)
    ##             v = m.group(2)
    ##             v = str.strip(v)

    ##             # units key separation
    ##             m = re.match("(.*)\((.*)\).*",k)
    ##             if m:
    ##                 u = str.strip(m.group(2))
    ##                 if u == "": u = "-"

    ##                 k = str.lower(str.strip(m.group(1)))
    ##                 self.units[k] = u
    ##             else:
    ##                 k = str.lower(str.strip(k))
    ##                 self.units[k] = "-"

    ##             self.value[k] = v
    ##             if k not in self.key_list: self.key_list.append(k)

    ##             # if the keyword eof is set, the end of the header of
    ##             #  an RDF/binary file has been reached. Exit parsing loop

    ##             if re.match(self.eof,k): break

    ##     return self

    def float(self, key):
        x = list(map(locale.atof, str.split(self.value[key])))
        if len(x) == 1:
            return x[0]
        else:
            return x

    def int(self, key):
        x = list(map(locale.atoi, str.split(self.value[key])))
        if len(x) == 1:
            return x[0]
        else:
            return x

    def long(self, key):
        x = list(map(locale.atol, str.split(self.value[key])))
        if len(x) == 1:
            return x[0]
        else:
            return x

    def strings(self, key, separator=None):
        """Return a stripped list of strings, originally separated by
        separator (default blank)"""

        if separator:
            strings = self.value[key].split(separator)
        else:
            strings = self.value[key].split()
        for i in range(len(strings)):
            strings[i] = strings[i].strip()
        return strings

    def printRDF(self):
        """Print RDF structure to stdout."""

        for k in self.key_list:
            print("%-30s%10s %s %s" % (k, '(' + self.units[k] + ')',
                                       self.rdfSep, self.value[k]))

    def writeRDF(self, fileOut, exclude=None):
        """Print RDF structure to file."""

        fout = open(fileOut, "w")
        for k in self.key_list:
            if exclude:
                if self.value[k] != exclude:
                    fout.write("%-30s%10s %s %s\n" %
                               (k, '(' + self.units[k] + ')', self.rdfSep,
                                self.value[k]))
            else:
                fout.write("%-30s%10s %s %s\n" % (k, '(' + self.units[k] + ')',
                                                  self.rdfSep, self.value[k]))
        fout.close()

    def writeTemplate(self, fileOut):
        """Print RDF structure to file."""

        fout = open(fileOut, "w")
        for k in self.key_list:
            fout.write("%-30s%10s %s \n" % (k, '(' + self.units[k] + ')',
                                            self.rdfSep))

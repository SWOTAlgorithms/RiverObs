#! /usr/bin/env python

from RDF import RDF

class RDF_to_class:
    """Convert the data in an RDF file or class to a class whose
    members are assigned desired variable names and corresponting
    values.

    To initialize the class a dictionary with the following structure

    d = {
    keyword1: (variableName1,format1,[default1]),
    keyword2: (variableName2,format2,[default2]),
    ...
    keywordN: (variableNameN,formatN,[defaultN])
    }

    where keywordI is the RDF keyword for the ith varaible (string)
    varaiableNameI is the name to be assigned to the ith variable (string)

    formatN is "s" or "f" or "d" depending whether the value returned
    is a string, float, or integer. (This also applies to arrays.)

    defaultI is an (optional) default value (string) which will be used
    to initialize the value. If no default value is specified, the
    value will be set to None until read from the RDF instance. If the
    RDF keyword is not present, no exception will be raised, and the
    value will be set to the default value

    As an example, if one entries is for class A is

    "x_keyword":("x","f","4.5")

    or
    
    "x_keyword":("x","f")
    
    the A.x will return the floating value of rdf["x_keyword"]

    Parameters
    ----------

    d : dictionary
        Dictionary containing the defintions, as above.
    file : str
        Path to a file containing RDF inputs.
    rdf : RDF instance
        Initialize from an RDF instance.

    Notes
    -----

    Intialize the variables to None. If rdf instance is given,
    load the values from the rdf structure. If a file name is given,
    the rdf is read from the file. file takes precedence over rdf.
    """

    def __init__(self,d, file=None,  rdf=None):
        self.d = d
        for tuple in d.values():
            if len(tuple) < 3:
                name,format = tuple
                exec("self.%s = None"%(name))
            else:
                name,format,default = tuple
                exec("self.%s = %s"%(name,default))

        if rdf != None: self.fromRDF(rdf)
        if file != None: self.fromFile(file)

    def fromRDF(self,rdf):
        """Read the values from an rdf instance."""

        for key in self.d.keys():
            format = self.d[key][1]
            try:
                if format == "s":
                    exec("self.%s = rdf[key]"%(self.d[key][0]))
                elif format == "d":
                    exec("self.%s = rdf.int(key)"%(self.d[key][0]))    
                elif format == "f":
                    exec("self.%s = rdf.float(key)"%(self.d[key][0]))
            except:
                name = self.d[key][0]
                if eval("self.%s"%name) == None:
                    raise ValueError("Cannot find keyword: %s"%key)
                

    def fromFile(self,file):
        """Read from a file."""

        self.fromRDF(RDF().rdfParse(file))

    
def test():

    rdf = RDF()
    rdf["key1"] = "1.2 2.3 3.4"
    rdf["key2"] = "1 2 3"
    rdf["key3"] = "1.2 2.3 3.4"

    d = {
        "key1":("x","f","None"),
        "key2":("y","d","None"),
        "key3":("z","s","None")
        }

    a = RDF_to_class(d)

    print a.__dict__

    a.fromRDF(rdf)

    print a.__dict__

    rdf["key1"] = "2.2 3.3 4.4"
    rdf["key2"] = "2 3 4"
    rdf["key3"] = "2.2 3.3 4.4"

    rdf.writeRDF("junk.dat")

    a.fromFile("junk.dat")

    print a.__dict__
    
if __name__ == '__main__': test()


            

        

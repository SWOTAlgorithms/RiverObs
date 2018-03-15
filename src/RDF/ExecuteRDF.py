"""
A python wrapper for calling a program that takes as input
an RDF file. This is supposed to be subclassed for different
programs.
"""

from __future__ import absolute_import, division, print_function

from threading import Thread
from subprocess import Popen, PIPE
from string import Template
from distutils.spawn import find_executable
from tempfile import NamedTemporaryFile


class ExecuteRDF(Thread):
    """A class meant to be subclassed to run a program that takes as an input
    an RDF file in a separate thread."""

    rdf_template = Template(
        "")  # This will need to be defined for derived classes

    def __init__(self, executable):
        """Intialize the thread. As an option, pass the path to the executable.
        If not passed, it will be searched for in the PATH."""

        self.stdout = None
        self.stderr = None
        Thread.__init__(self)

        if find_executable(executable) == None:
            raise Exception('Cannot find executable: %s' % executable)

        self.executable = executable

    def set_params(self, rdf=None, **kwargs):
        """Set the parameters in the template by passing an RDF
        object and/or kwargs. kwargs override the rdf object."""

        if rdf != None:
            for k in list(rdf.keys()):
                #exec('self.%s = rdf["%s"]'%(k,k))
                setattr(self, k, k)

        for k in kwargs:
            #exec('self.%s = kwargs["%s"]'%(k,k))
            setattr(self, k, kwargs[k])

    def set_template(self, rdf_template):
        """Set the rdf template by passing a string with $variable
        substitutions corresponding to the class variables set
        using set_params."""

        self.rdf_template = Template(rdf_template)

    def write_rdf(self, rdf_file=None, dir=None, delete=True):
        """Open an RDF file and write the rdf inputs. If rdf_file = None,
        a temporary file will be created. It is assumed that all appropriate
        parameters have been defined in set_params."""

        if rdf_file == None:
            self.fout_rdf = NamedTemporaryFile(delete=delete, dir=dir)
            self.rdf_name = self.fout_rdf.name
        else:
            self.fout_rdf = open(rdf_file, 'w')
            self.rdf_name = rdf_file

        self.fout_rdf.write(self.rdf_template.substitute(self.__dict__))
        self.fout_rdf.flush()

    def run(self, args=[], args_start=True):
        """This is the thread definition. args is a list of optional arguments
        to be passed to the executable by Popen. If args_start == True, the
        call sequence is 'executable args rdf_file'. Otherwise, it is
        'executable args rdf_file'."""

        command = [self.executable]
        if args != [] and args_start:
            command += args
        command.append(self.rdf_name)
        if args != [] and not args_start:
            command += args

        print(command)
        p = Popen(command, shell=False, stdout=PIPE, stderr=PIPE)

        self.stdout, self.stderr = p.communicate()

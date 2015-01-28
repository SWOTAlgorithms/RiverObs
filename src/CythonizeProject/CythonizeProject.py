"""
Make a parallel project tree, make desired files into C files, and copy other
files over as Python files.
"""

import os, os.path
from subprocess import call
from Cython.Build import cythonize

class CythonizeProject:
    """Make a parallel project tree, make desired files into C files, and copy other
    files over as Python files.

    Initialize Cython project tree.

    Parameters
    ----------

    included_dirs : list
        List of included directory paths to search for code.
    build_dir : str
        Where the sources are compiled. 
    ignored_files : list, default ['__init__.py','version.py']
        Files not to turn into C but copy as python.
    excluded_dirs : list, default ['scripts']
        List of relative directory names to exclude
    """
    def __init__(self,included_dirs,build_dir,ignored_files=['__init__.py','version.py'],
                 excluded_dirs=['scripts']):

        self.included_dirs = included_dirs
        self.build_dir = build_dir
        self.ignored_files = ignored_files
        self.excluded_dirs = excluded_dirs

    def make_extensions(self,force_compile=True):
        """Return a list with all of the extensions for this package."""

        modules = []
        for dir_ in self.included_dirs:
            for dirname, dirnames, filenames in os.walk(dir_):
                for name in self.excluded_dirs:
                    if name in dirnames:
                        dirnames.remove(name)
                for filename in filenames:
                    file_ = os.path.join(dirname, filename)
                    stripped_name = os.path.split(file_)[-1] #os.path.relpath(file_, code_python_path)
                    file_name, extension = os.path.splitext(stripped_name)
                    print file_
                    if (extension == '.py') or (extension == '.pyx'):
                        if filename not in self.ignored_files:
                            modules += cythonize(file_,
                                            force=force_compile,
                                            build_dir=self.build_dir)
                            
                        else: # Copy the file to the source directory
                            file_copy = os.path.join(self.build_dir,file_)
                            cp_dir = os.path.split(file_copy)[0]
                            if not os.path.exists(cp_dir):
                                os.makedirs(cp_dir)
                            call(['cp',file_,file_copy])

        return modules
                            
                            

                            
  

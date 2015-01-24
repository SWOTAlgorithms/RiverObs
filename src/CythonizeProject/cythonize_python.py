"""
Based on Julia Sanchez Simon:

http://blog.biicode.com/bii-internals-compiling-your-python-application-with-cython/
"""

import os.path
from Cython.Build import cythonize

def cythonize_python(included_dirs,code_python_path,src_dir,
                     build_dir=None,
                     ignored_files=['__init__.py'],
                     excluded_dirs=['scripts'],force_compile=True):
    """
    Traverse your source tree, for every file cython will analyze its
    dependencies tree and recompile only if needed.

    Creates C files from your source python.

    Parameters
    ----------

    included_dirs : list
        List of included directory paths.
    code_python_path : str
        Path to python code base directory
    src_dir : str
        Path to where the C files end up.
    build_dir : str, defaut None
        Where the sources are compiled. If None, then use src_dir
    ignored_files : list, default ['__init__.py']
        Files not to turn into C.
    excluded_dirs : list, default ['scripts']
        List of relative directory names to exclude
    force_compile : boolean, default True
        if true compiles regardeless of whether the file has changed or not.

    Returns
    -------
        list of c files relative to pkg_path
    """

    if build_dir == None:
        build_dir = src_dir

    if not os.path.exists(src_dir):
        os.makedirs(src_dir)
    
    c_files = []
    for dir_ in included_dirs:
        for dirname, dirnames, filenames in os.walk(dir_):
            for name in excluded_dirs:
                if name in dirnames:
                    dirnames.remove(name)
            for filename in filenames:
                file_ = os.path.join(dirname, filename)
                stripped_name = os.path.relpath(file_, code_python_path)
                file_name, extension = os.path.splitext(stripped_name)
                if extension == '.py':
                    target_file = os.path.join(src_dir, file_name + '.c')
                    if filename not in ignored_files:
                        c_files.append(stripped_name.replace('.py', '.c'))
                        file_dir = os.path.dirname(target_file)
                        if not os.path.exists(file_dir):
                            os.makedirs(file_dir)
                        print stripped_name,build_dir
                        extension = cythonize(file_, #stripped_name,
                                                force=force_compile,
                                                build_dir=build_dir)
    return c_files


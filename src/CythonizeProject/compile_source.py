"""
Now that we have the list of C files we want to compile them into python extensions but before
we may want to ignore python asserts and remove debug flags.
"""

import os.path
import platform
import os.path
from distutils import sysconfig

def compile_source(abs_path_c_files):
    """Depending on the system, compile the C files created by cythonize.

    Parameters
    ----------

    abs_path_c_files : list
        List of paths to the C files.    
    """
    
    modules = []
    for c_file in abs_path_c_files:
        relfile = os.path.relpath(c_file, src_dir)
        filename = os.path.splitext(relfile)[0]
        extName = filename.replace(os.path.sep, ".")
        extension = Extension(extName,
                            sources=[c_file],
                            define_macros=[('PYREX_WITHOUT_ASSERTIONS',
                            None)] # ignore asserts in code
                            )
        modules.append(extension)
        
    if platform.system() != 'Windows':
        cflags = sysconfig.get_config_var('CFLAGS')
        opt = sysconfig.get_config_var('OPT')
        sysconfig._config_vars['CFLAGS'] = cflags.replace(' -g ', ' ')
        sysconfig._config_vars['OPT'] = opt.replace(' -g ', ' ')
    if platform.system() == 'Linux':
        ldshared = sysconfig.get_config_var('LDSHARED')
        sysconfig._config_vars['LDSHARED'] = ldshared.replace(' -g ', ' ')
    elif platform.system() == 'Darwin':
        #-mno-fused-madd is a deprecated flag that now causes a hard error
        # but distuitls still keeps it
        # it was used to disable the generation of the fused multiply/add instruction
        for flag, flags_line in sysconfig._config_vars.iteritems():
            if ' -g' in str(flags_line):
                sysconfig._config_vars[flag] = flags_line.replace(' -g', '')
        for key in ['CONFIG_ARGS', 'LIBTOOL', 'PY_CFLAGS', 'CFLAGS']:
            value = sysconfig.get_config_var(key)
            if value:
                sysconfig._config_vars[key] = value.replace('-mno-fused-madd', '')
                sysconfig._config_vars[key] = value.replace('-DENABLE_DTRACE',  '')
                

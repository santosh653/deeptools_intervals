#!/usr/bin/env python
from setuptools import setup, Extension, find_packages
from distutils import sysconfig
import subprocess
import glob
import sys

srcs = [x for x in glob.glob("*.c")]

libs=["z"]
if sysconfig.get_config_vars('BLDLIBRARY') is not None:
    #Note the "-l" prefix!
    for e in sysconfig.get_config_vars('BLDLIBRARY')[0].split():
        if e[0:2] == "-l":
            libs.append(e[2:])
elif(sys.version_info[0] >= 3 and sys.version_info[1] >= 3) :
    libs.append("python%i.%im" % (sys.version_info[0], sys.version_info[1]))
else :
    libs.append("python%i.%i" % (sys.version_info[0], sys.version_info[1]))

additional_libs = [sysconfig.get_config_var("LIBDIR"), sysconfig.get_config_var("LIBPL")]

try:
    foo, _ = subprocess.Popen(['pcre-config', '--libs'], stdout=subprocess.PIPE).communicate()
except:
    sys.exit("Either libpcre isn't installed, it didn't come with pcre-config, or pcre-config isn't in your $PATH. This must be corrected before installing deeptoolsinterals!\n")

foo = foo.strip().split()
for v in foo:
    if(v[0:2] == "-L") :
        additional_libs.append(v[2:])
    elif(v[0:2] == "-l") :
        libs.append(v[2:])

module1 = Extension('deeptoolsintervals',
                    sources = srcs,
                    libraries = libs,
                    library_dirs = additional_libs, 
                    include_dirs = [sysconfig.get_config_var("INCLUDEPY")])

setup(name = 'deeptoolsintervals',
       version = '0.1.0',
       description = 'A python module creating/accessing GTF-based interval trees with associated meta-data',
       author = "Devon P. Ryan",
       author_email = "ryan@ie-freiburg.mpg.de",
       url = "https://github.com/dpryan79/deeptools_intervals",
       keywords = ["bioinformatics", "GTF"],
       classifier = ["Development Status :: 2 - Pre-Alpha",
                     "Environment :: Console",
                     "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
                     "Intended Audience :: Developers",
                     "Programming Language :: C",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 2",
                     "Programming Language :: Python :: 2.7",
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.3",
                     "Programming Language :: Python :: 3.4",
                     "Programming Language :: Python :: 3.5",
                     "Programming Language :: Python :: Implementation :: CPython",
                     "Operating System :: POSIX",
                     "Operating System :: Unix",
                     "Operating System :: MacOS",
                     "Topic :: Scientific/Engineering"],
       packages = find_packages(),
       include_package_data=True,
       ext_modules = [module1])

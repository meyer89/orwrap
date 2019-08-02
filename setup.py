import os
import sys
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import ortools


# the compilation is not compatible with windows
if not sys.platform.startswith('win'):
    extensions = [
        Extension("orwrap.matrix_constraints", ["orwrap/matrix_constraints.pyx"],
            include_dirs=[os.path.abspath("include")],
            libraries=["ortools"],
            library_dirs=[os.path.join(ortools.__path__[0], ".libs")],
            language="c++"
        ),
    ]
else:
    extensions = []


# read the long readme file
try:
    with open('README.rst', 'rb') as f:
        long_description = f.read().decode('utf-8')
except FileNotFoundError:
    long_description = None


# call the setup function for the compilation
setup(
    name="orwrap",
    version="19.8.2",
    author="Marian Meyer",
    ext_modules=cythonize(extensions),
    install_requires=["ortools==" + ortools.__version__, "numpy", "scipy"],
    description="Advanced wrapper around ortools",
    long_description=long_description,
    packages=find_packages()
)

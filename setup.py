import os
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import ortools

extensions = [
    Extension("ortools_ext.matrix_constraints", ["ortools_ext/matrix_constraints.pyx"],
        include_dirs=[os.path.abspath("include")],
        libraries=["ortools"],
        library_dirs=[os.path.join(ortools.__path__[0], ".libs")],
        language="c++"
    ),
]

# call the setup function for the compilation
setup(
    name="ortools_ext",
    version=ortools.__version__,
    author="Marian Meyer",
    ext_modules=cythonize(extensions),
    install_requires=["ortools==" + ortools.__version__],
    packages=find_packages()
)

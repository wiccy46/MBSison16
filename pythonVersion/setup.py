
from distutils.core import setup, Extension
import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext


"""
Usage:
python setup.py build_ext --inplace
1. use the first setup
2. the send setup
"""

# setup(
#     ext_modules=cythonize("Trajectory.pyx"),
#     include_path=[numpy.get_include()]
# )


setup(
    ext_modules=[
        Extension("Trajectory", ["Trajectory.c"],
                  include_dirs=[numpy.get_include()]),
    ],
)



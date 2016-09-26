
from distutils.core import setup
import numpy, os
from Cython.Build import cythonize

# To Build in Terminal: python setup.py build_ext --inplace


os.environ["CPPFLAGS"] = os.getenv("CPPFLAGS", "")  + "-I" + numpy.get_include()

setup(
    ext_modules=cythonize("Trajectory.pyx")
)





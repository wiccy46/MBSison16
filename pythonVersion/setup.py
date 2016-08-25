
from distutils.core import setup
import numpy, os
from Cython.Build import cythonize


os.environ["CPPFLAGS"] = os.getenv("CPPFLAGS", "")  + "-I" + numpy.get_include()

setup(
    ext_modules=cythonize("Trajectory.pyx")
)





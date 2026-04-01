from setuptools import setup, Extension
import pybind11


ext_modules = [
    Extension(
        'starmatcher',            # your module name
        ['starmatcher.cpp'],      # just the filename because it's in the same folder
        include_dirs=[pybind11.get_include()],
        language='c++'
    ),
]

setup(
    name='starmatcher',
    version='0.1',
    ext_modules=ext_modules,
)

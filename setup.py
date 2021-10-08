from setuptools import setup, find_packages
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
from distutils import ccompiler
import os
import pysam
import glob
import sys


# Note building htslib for OSX version might need to be set: make CXXFLAGS="-mmacosx-version-min=10.09"

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


# This was stolen from pybind11
# https://github.com/pybind/python_example/blob/master/setup.py
# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler, flags):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    for flag in flags:
        if has_flag(compiler, flag):
            return flag


def get_extra_args():
    compiler = ccompiler.new_compiler()
    extra_compile_args = []
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']
    f = cpp_flag(compiler, flags)
    if not f:
        return ['-std=c++11']  # raise RuntimeError("Invalid compiler")
    extra_compile_args.append(f)
    flags = ['--stdlib=libc++']
    f = cpp_flag(compiler, flags)
    if f:
        extra_compile_args.append(f)

    return extra_compile_args


extras = get_extra_args() + ["-Wno-sign-compare", "-Wno-unused-function",
                             "-Wno-unused-result", '-Wno-ignored-qualifiers',
                             "-Wno-deprecated-declarations"
                             ]

ext_modules = list()

root = os.path.abspath(os.path.dirname(__file__))
# Try and link dynamically to htslib
htslib = None
if "--htslib" in sys.argv:
    idx = sys.argv.index("--htslib")
    h = sys.argv[idx + 1]
    if h and os.path.exists(h):
        if any("libhts" in i for i in glob.glob(h + "/*")):
            print("Using --htslib at {}".format(h))
            htslib = h
            if htslib[-1] == "/":
                htslib = htslib[:-1]
            sys.argv.remove("--htslib")
            sys.argv.remove(h)
    else:
        raise ValueError("--htslib path does not exists")

#
# if os.getenv("PATH", None):
#     for i in os.environ["PATH"].split(":"):
#         if "htslib" in os.path.basename(i):
#             print(glob.glob(i))
#             if any("libhts" in i for i in glob.glob(i + "/*")):
#                 htslib = i
#                 break

if htslib is None:
    print("Using packaged htslib")
    htslib = os.path.join(root, "dysgu/htslib")


libraries = [f"{htslib}/hts"]  # Library name for libhts.so
library_dirs = [htslib, numpy.get_include(), f"{htslib}/htslib"] + pysam.get_include()
include_dirs = [numpy.get_include(), root, # htslib,
                f"{htslib}/htslib", f"{htslib}/cram"] + pysam.get_include()
runtime_dirs = [htslib]  # os.path.join(root, "dysgu/htslib")

# extra_link_args = ["-L ./dysgu/htslib", "-L ./dysgu/htslib/htslib"]

print("Libs", libraries)
print("Library dirs", library_dirs)
print("Include dirs", include_dirs)
print("Runtime dirs", runtime_dirs)
print("Extras compiler args", extras)

for item in ["sv2bam", "io_funcs", "graph", "coverage", "assembler", "call_component",
             "map_set_utils", "cluster", "post_call_metrics", "sv_category"]:

    ext_modules.append(Extension(f"dysgu.{item}",
                                 [f"dysgu/{item}.pyx"],
                                 libraries=libraries,
                                 library_dirs=library_dirs,
                                 include_dirs=include_dirs,
                                 runtime_library_dirs=runtime_dirs,
                                 extra_compile_args=extras,
                                 # extra_link_args=extra_link_args,
                                 define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                                 language="c++"))


print("Found packages", find_packages(where="."))
setup(
    name="dysgu",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    url="https://github.com/kcleal/dysgu",
    description="Structural variant calling",
    license="MIT",
    version='1.3.0',
    python_requires='>=3.7',
    install_requires=[  # runtime requires
            'cython',
            'click',
            'numpy',
            'pandas',
            'pysam',
            'networkx>=2.4',
            'scikit-learn',
            'ncls',
            'scikit-bio',
            'sortedcontainers',
            'lightgbm',
            'edlib',

        ],
    setup_requires=[
            'cython',
            'click',
            'numpy',
            'pandas',
            'pysam',
            'networkx>=2.4',
            'scikit-learn',
            'ncls',
            'scikit-bio',
            'sortedcontainers',
            'lightgbm',
            'edlib'],  # setup.py requires at compile time
    packages=["dysgu", "dysgu.tests"],
    ext_modules=cythonize(ext_modules),

    include_package_data=True,
    zip_safe=False,
    entry_points='''
        [console_scripts]
        dysgu=dysgu.main:cli
    ''',
)

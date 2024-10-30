from setuptools import setup, Extension
import setuptools
from Cython.Build import cythonize
import numpy
import os
import sys
import glob
import pysam
import sysconfig

cfg_vars = sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if isinstance(value, str):
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

def has_flag(compiler, flagname):
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

def cpp_flag(compiler, flags):
    for flag in flags:
        if has_flag(compiler, flag):
            return flag

def get_extra_args():
    from distutils import ccompiler
    compiler = ccompiler.new_compiler()
    extra_compile_args = []
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']
    f = cpp_flag(compiler, flags)
    if not f:
        return ['-std=c++11']
    extra_compile_args.append(f)
    flags = ['-stdlib=libc++']
    f = cpp_flag(compiler, flags)
    if f:
        extra_compile_args.append(f)
    return extra_compile_args

extras = get_extra_args() + ["-Wno-sign-compare", "-Wno-unused-function",
                             "-Wno-unused-result", '-Wno-ignored-qualifiers',
                             "-Wno-deprecated-declarations", "-fpermissive",
                             "-Wno-unreachable-code-fallthrough",
                             ]

def get_extension_modules():
    ext_modules = []

    # root = os.path.abspath(os.path.dirname(__file__))
    root = os.path.dirname(__file__)
    libraries, library_dirs, include_dirs, runtime_dirs = [], [], [], []

    if "--conda-prefix" in sys.argv or os.getenv('PREFIX'):
        prefix = os.getenv('PREFIX') if not "--conda-prefix" in sys.argv else sys.argv[sys.argv.index("--conda-prefix") + 1]
        if prefix and os.path.exists(prefix):
            if any("libhts" in i for i in glob.glob(prefix + "/lib/*")):
                print(f"Using htslib at {prefix}")
                if prefix[-1] == "/":
                    prefix = prefix[:-1]
            else:
                raise ValueError(f"libhts not found at {prefix}/lib/*")
        else:
            raise ValueError("prefix path does not exists")
        libraries = ["hts"]
        library_dirs = [f"{prefix}/lib", numpy.get_include()] + pysam.get_include()
        include_dirs = [numpy.get_include(), ".",
                        f"{prefix}/include/htslib", f"{prefix}/include"] + pysam.get_include()
        runtime_dirs = [f"{prefix}/lib"]
    else:
        htslib = os.getenv('HTSLIB_DIR') if not "--htslib" in sys.argv else sys.argv[sys.argv.index("--htslib") + 1]
        if htslib and os.path.exists(htslib):
            if any("libhts" in i for i in glob.glob(htslib + "/*")):
                print(f"Using --htslib at {htslib}")
                if htslib[-1] == "/":
                    htslib = htslib[:-1]
            else:
                raise ValueError("--htslib path does not exists")
        else:
            print("Using packaged htslib")
            htslib = os.path.join(root, "dysgu/htslib")
        libraries = ["hts"]
        library_dirs = [htslib, numpy.get_include(), f"{htslib}/htslib"] + pysam.get_include()
        include_dirs = [numpy.get_include(), ".",
                        f"{htslib}/htslib", f"{htslib}/cram"] + pysam.get_include()
        runtime_dirs = [htslib]

    ext_modules.append(Extension("dysgu.scikitbio._ssw_wrapper",
                                 ["dysgu/scikitbio/_ssw_wrapper.pyx", "dysgu/scikitbio/ssw.c"],
                                 include_dirs=["dysgu/scikitbio", numpy.get_include()],
                                 extra_compile_args=["-Wno-deprecated-declarations", '-std=c99', '-I.'],
                                 language="c"))

    ext_modules.append(Extension("dysgu.edlib.edlib",
                                 ["dysgu/edlib/edlib.pyx", "dysgu/edlib/src/edlib.cpp"],
                                 include_dirs=["dysgu/edlib", numpy.get_include()],
                                 extra_compile_args=["-O3", "-std=c++11"],
                                 language="c++"))

    ext_modules.append(Extension("dysgu.sortedintersect.sintersect",
                                 ["dysgu/sortedintersect/sintersect.pyx"],
                                 extra_compile_args=["-O3", "-std=c++11"],
                                 language="c++"))

    for item in ["sv2bam", "io_funcs", "graph", "coverage", "consensus", "call_component",
                 "map_set_utils", "cluster", "sv_category", "extra_metrics", "merge_svs"]:
        ext_modules.append(Extension(f"dysgu.{item}",
                                     [f"dysgu/{item}.pyx"],
                                     libraries=libraries,
                                     library_dirs=library_dirs,
                                     include_dirs=include_dirs,
                                     runtime_library_dirs=runtime_dirs,
                                     extra_compile_args=extras,
                                     define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                                     language="c++"))

    return cythonize(ext_modules)

setup(
    name="dysgu",
    ext_modules=get_extension_modules(),
)
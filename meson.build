project(
  'dysgu', ['c', 'cpp'],
  meson_version  : '>=0.64',
  default_options: [
    'warning_level=1',
    'cpp_std=c++17',
    'c_std=c99'
  ]
)

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
py_mod = import('python')

cython = find_program('cython')
if run_command(cython, '--version',
                check: true
                ).stdout().strip().split()[-1].version_compare('<3.0.0')
  error('Cython ≥3.0.0 is required (needed for NumPy2 support)')
endif

# ---------------------------------------------------------------------------
# Python interpreter and numpy / pysam discovery
# ---------------------------------------------------------------------------
py3     = py_mod.find_installation(pure : false)   # ← the interpreter that is running Meson
py3_dep = py3.dependency()

python_info = run_command(
               py3,
               ['-c', 'import sys; print(f"Using Python {sys.version} at {sys.executable}")'],
               check : true
             ).stdout().strip()
message(python_info)

numpy_include   = run_command(py3,
                    ['-c', 'import numpy, sys; print(numpy.get_include())'],
                    check : true).stdout().strip()

pysam_includes  = run_command(py3,
                    ['-c', 'import pysam, sys; print(":".join(pysam.get_include()))'],
                    check : true).stdout().strip().split(':')


# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------
extra_args  = [
  '-Wno-sign-compare', '-Wno-unused-function',
  '-Wno-unused-result', '-Wno-ignored-qualifiers',
  '-Wno-deprecated-declarations'
]
cpp_defines = ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']

# ---------------------------------------------------------------------------
# htslib
# ---------------------------------------------------------------------------

htslib_prefix = get_option('htslib_prefix')

if htslib_prefix != ''
  # Use the system-provided htslib
  hts_lib = declare_dependency(
    include_directories : include_directories(
                            join_paths(htslib_prefix, 'include'),
                            join_paths(htslib_prefix, 'include/htslib')),
    link_args           : ['-L' + join_paths(htslib_prefix, 'lib'), '-lhts'])
else
  # Try to get libdeflate, but make it optional
  libdeflate_dep = dependency('libdeflate', required: false)

  # Set htslib options based on whether libdeflate was found
  htslib_options = []
  if libdeflate_dep.found()
    message('libdeflate found, enabling it for htslib')
    htslib_options += 'libdeflate=true'
  else
    message('libdeflate not found, building htslib without it')
    htslib_options += 'libdeflate=false'
  endif

  # Build htslib as a dependency
  hts_lib = dependency('htslib',
                       fallback : ['htslib', 'htslib_dep'],
                       required : true,
                       default_options: htslib_options)
endif


# ---------------------------------------------------------------------------
# Common include directories
# ---------------------------------------------------------------------------
includes = [
  include_directories('.'),
  include_directories(numpy_include),
]

foreach inc : pysam_includes
  includes += include_directories(inc)
endforeach

includes += include_directories('dysgu')

# ---------------------------------------------------------------------------
# Additional source files
# ---------------------------------------------------------------------------
extra_headers = [
  'find_reads.hpp',
  'graph_objects.hpp'
]

foreach header : extra_headers
  install_headers(join_paths('dysgu', header), subdir : 'dysgu')
endforeach

# ---------------------------------------------------------------------------
# Cython handling and module definitions
# ---------------------------------------------------------------------------

# Define the main Cython modules (those in the main dysgu directory)
main_modules = [
  'sv2bam', 'io_funcs', 'graph', 'coverage', 'consensus',
  'call_component', 'map_set_utils', 'cluster',
  'sv_category', 'extra_metrics', 'merge_svs'
]

# Generate Cython targets for all modules
cython_targets = {}

# Special modules first
# edlib module
cython_targets += {
  'edlib' : custom_target(
    'cythonize_edlib',
    input           : join_paths('dysgu', 'edlib', 'edlib.pyx'),
    output          : 'edlib.cpp',
    command         : [
      cython,
      '--cplus',
      '-3',
      '-I', meson.current_source_dir(),
      '@INPUT@',
      '-o', '@OUTPUT@'
    ],
  )
}

# ssw_wrapper module
cython_targets += {
  '_ssw_wrapper' : custom_target(
    'cythonize_ssw_wrapper',
    input           : join_paths('dysgu', 'scikitbio', '_ssw_wrapper.pyx'),
    output          : '_ssw_wrapper.c',
    command         : [
      cython,
      '-3',
      '-I', meson.current_source_dir(),
      '@INPUT@',
      '-o', '@OUTPUT@'
    ],
  )
}

# Main modules
foreach module : main_modules
  cython_targets += {
    module : custom_target(
      'cythonize_' + module,
      input           : join_paths('dysgu', module + '.pyx'),
      output          : module + '.cpp',
      command         : [
        cython,
        '--cplus',
        '-3',
        '-I', meson.current_source_dir(),
        '-I', join_paths(meson.current_source_dir(), 'dysgu'),
        '@INPUT@',
        '-o', '@OUTPUT@'
      ],
    )
  }
endforeach

# ---------------------------------------------------------------------------
# Extension modules
# ---------------------------------------------------------------------------
# Special case modules
py3.extension_module(
  'edlib',
  [cython_targets['edlib'], join_paths('dysgu', 'edlib', 'src', 'edlib.cpp')],
  include_directories : [includes, include_directories('dysgu/edlib')],
  cpp_args            : ['-O3', '-std=c++17'] + cpp_defines,
  install             : true,
  install_dir         : py3.get_install_dir() / 'dysgu/edlib'
)

py3.extension_module(
  '_ssw_wrapper',
  [cython_targets['_ssw_wrapper'], join_paths('dysgu', 'scikitbio', 'ssw.c')],
  include_directories : [includes, include_directories('dysgu/scikitbio')],
  c_args              : ['-Wno-deprecated-declarations'] + cpp_defines,
  install             : true,
  install_dir         : py3.get_install_dir() / 'dysgu/scikitbio'
)

# Main modules
foreach module : main_modules
  py3.extension_module(
    module,
    [cython_targets[module]],
    include_directories : includes,
    dependencies        : [py3_dep, hts_lib],
    cpp_args            : extra_args + cpp_defines,
    install             : true,
    install_dir         : py3.get_install_dir() / 'dysgu'
  )
endforeach

# ---------------------------------------------------------------------------
# Install Python source files
# ---------------------------------------------------------------------------
# Absolute path where Python packages go (e.g. site‑packages)
py_pkg_root = py3.get_install_dir() / 'dysgu'

python_sources = [
  '__init__.py', 'main.py', 'python_api.py',
  'sites_utils.py', 'filter_normals.py', 'post_call.py',
  're_map.py', 'view.py'
]

# Install Python sources
foreach source : python_sources
  py3.install_sources(join_paths('dysgu', source), subdir : 'dysgu')
endforeach

# Install data files
install_data(
  join_paths('dysgu', 'dysgu_model.1.pkl.gz'),
  install_dir : py_pkg_root
)

# Install sub‑packages
foreach pkg : ['tests', 'scikitbio', 'edlib', 'gfftool']
  # copies dysgu/<pkg>/… into <site‑packages>/dysgu/<pkg>/…
  install_subdir(
    join_paths('dysgu', pkg),
    install_dir : py_pkg_root,
    exclude_files : ['*.pyx', '*.c', '*.cpp', '*.h', '*.hpp']
  )
endforeach
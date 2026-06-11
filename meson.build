project('guptri_py', 'c',
  meson_version: '>=1.1.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran', native: false)

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_project_arguments(
    '-fallow-argument-mismatch', # Downgrades Type/Rank mismatch errors
    '-w',                        # Suppresses all warnings (-Wcompare-reals, etc.)
    language: 'fortran'
  )
endif

py = import('python').find_installation(pure: false)

# NumPy >=2.0 provides include dirs via dependency()
np_dep = dependency('numpy')

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true
).stdout().strip()

# f2py include dir (for fortranobject.h) is not in dependency('numpy'),
# so add it separately
f2py_dep = declare_dependency(
  include_directories : incdir_f2py,
)

# Generate the f2py wrappers
guptri_source = custom_target('fguptrimodule.c',
  input : ['guptri/fguptri.f'],
  output : ['fguptrimodule.c', 'fguptri-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'fguptri', '--lower']
)

py.extension_module('fguptri',
  ['guptri/fguptri.f', 'guptri/guptribase.f', 'guptri/zguptri.f', guptri_source, incdir_f2py / 'fortranobject.c'],
  dependencies : [np_dep, f2py_dep],
  subdir: 'guptri_py',
  install : true,
)

# Install the Python package files
py.install_sources(
  'guptri_py/__init__.py',
  subdir: 'guptri_py',
)

project('GLOW', 'fortran',
  version : '0.982',
  meson_version : '>= 0.50',
  default_options : ['default_library=static', 'buildtype=release'])

fc = meson.get_compiler('fortran')
fargs_old = []
fargs_new = []
if fc.get_id() == 'gcc'
  if fc.version().version_compare('>=8')
    fargs_old = ['-w', '-std=legacy']
    fargs_new += '-std=f2018'
  else
    fargs_old += '-w'
  endif
elif fc.get_id() == 'intel'
  fargs_old += '-w'
  fargs_new = ['-stand f18', '-warn declarations', '-W1']
elif fc.get_id() == 'intel-cl'
  fargs_old += '/w'
  fargs_new = ['/stand:f18', '/warn:declarations', '/W1']
elif fc.get_id() == 'pgi'
  fargs_old += '-w'
  fargs_new +=  '-Mdclchk'
elif fc.get_id() == 'flang'
  fargs_old += '-w'
endif


cglow_src = files('glow.f90', 'cglow.f90', 'fieldm.f', 'solzen.f90', 'ssflux.f90', 'rcolum.f90', 'qback.f90', 'qproton.f90', 'etrans.f90', 'exsect.f', 'gchem.f90', 'bands.f90', 'ephoto.f90', 'egrid.f90', 'maxt.f90')
snoem_src = files('snoem.f90', 'snoemint.f90', 'geomag.f90')
utils_src = files('mzgrid.f90', 'conduct.f90', 'nrlmsise00.f', 'iri90.f')

cglow = library('cglow', cglow_src, snoem_src, utils_src, fortran_args : fargs_old)

#if meson.is_subproject()
  subdir_done()
#endif

# ------  Python Matlab command line driver -----------
exe_name = 'glow.bin'
glowexe = executable(exe_name,
  'glowpython.f90', 'utils.f90', 'fsutils.f90',
  link_with : cglow,
  fortran_args : fargs_new)

# -- hack: an always dirty target to be like RUNTIME_OUTPUT_DIRECTORY on CMake
if host_machine.system() == 'windows'
  exe_name += '.exe'
endif

custom_target('glow_runtime_output_directory',
  build_by_default: true,
  output: 'dummy',
  command: ['python', '-c', 'import shutil; shutil.copy("' + meson.current_build_dir() / exe_name + '", "' + meson.current_build_dir() / '..' + '")'],
  depends: glowexe)
# -- end hack

mpi = dependency('mpi', language : 'fortran', required : false)
netcdf = dependency('netcdf', required : false)

if mpi.found() and netcdf.found()
mpiglow = executable('mpi_glow.bin',
  'glowdriver.f90', 'readtgcm.f90', 'output.f90', 'tzgrid.f90','fsutils.f90',
  dependencies : [netcdf, mpi],
  link_with : cglow,
  fortran_args : fargs_sold)
endif


# Octave and Matlab tests can be squirrely, make user request to test them.
if meson.is_subproject() or not get_option('testmatlab')
  subdir_done()
endif

margs = {'matlab': '-batch', 'octave': '--eval'}

foreach m : ['matlab', 'octave']

  mexe = find_program(m, required : false)
  if not mexe.found()
    continue
  endif

  test('Simple-' + m, mexe,
    args : [margs[m], 'Simple'],
    workdir : meson.current_source_dir() / 'matlab',
    timeout: 60)

endforeach

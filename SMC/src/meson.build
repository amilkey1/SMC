project('new-smc', 'cpp',
	default_options : ['cpp_std=c++11','buildtype=release','prefix=/home/aam21005/test-new-smc'],
	version : '1.0')
cpp = meson.get_compiler('cpp')

# These lines specify the locations of the NCL, Boost compiled libraries
lib_system = cpp.find_library('boost_system', dirs: ['/home/aam21005/lib/static'], required: true)
lib_filesystem = cpp.find_library('boost_filesystem', dirs: ['/home/aam21005/lib/static'], required: true)
lib_program_options = cpp.find_library('boost_program_options', dirs: ['/home/aam21005/lib/static'], required: true)
lib_ncl = cpp.find_library('ncl', dirs: ['/home/aam21005/lib/static'], required: true)

# These lines specify the locations of header files for the NCL, Boost
incl_ncl = include_directories('/home/aam21005/include')
incl_boost = include_directories('/home/aam21005/Documents/libraries/boost_1_77_0')

# This line creates the executable file
executable('new-smc', 'main.cpp', link_args:'-pthread', install: true, install_dir: '.', dependencies: [lib_ncl,lib_program_options,lib_system,lib_filesystem], include_directories: [incl_ncl,incl_boost])

# These lines just copy files to the install directory
#install_data('sim_network_data.nex', install_dir: '.')
install_data('s.sh', install_dir: '.')
install_data('test.sh', install_dir: '.')

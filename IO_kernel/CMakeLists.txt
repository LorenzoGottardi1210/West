# ==============================================================================
# WEST - IO_kernel/CMakeLists.txt
# ==============================================================================

# 1. Creating the target
add_library(west_io_kernel STATIC
    cubefile.f90
    function3d.f90
    mod_west_io.f90
    pdep_db.f90
    pdep_io.f90
    plep_db.f90
    plep_io.f90
    qbox_interface.f90
    wbse_io.f90
    wfreq_db.f90
    wfreq_io.f90
)

# 2. apply the Fortran preprocessor to the target
set_target_properties(west_io_kernel PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting the include directories for the target
target_include_directories(west_io_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. Linking of the target with the libraries of WEST, QE already inside the WEST build tree
target_link_libraries(west_io_kernel PUBLIC 
    west_modules
    west_tools
    west_fft_kernel
    west_para_kernel
)
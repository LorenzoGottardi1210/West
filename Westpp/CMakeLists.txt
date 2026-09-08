# ==============================================================================
# WEST - Westpp/CMakeLists.txt
# ==============================================================================

set(WESTPP_LIB_SOURCES
    do_dip.f90
    do_eigenpot2.f90
    do_exc_comp.f90
    do_exc_drho.f90
    do_exc_resp.f90
    do_exc_spin.f90
    do_loc.f90
    do_rho.f90
    do_sxx.f90
    do_wann.f90
    do_wfc2.f90
    dump_r.f90
    westpp_setup.f90
    write_wfc_1d_r.f90
    write_wfc_spav.f90
)

# 1. static library west_westpp 
add_library(west_westpp STATIC ${WESTPP_LIB_SOURCES})
set_target_properties(west_westpp PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "westpp"
)
target_include_directories(west_westpp PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# linking the library with all the other libraries of WEST, so that it can use all the functionalities of WEST 
target_link_libraries(west_westpp PUBLIC 
    west_wbse
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)

# 2. executables westpp
add_executable(westpp westpp.f90)
set_target_properties(westpp PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "westpp"
    SUFFIX ".x"
)

# executable linking 
target_link_libraries(westpp PRIVATE 
    west_westpp
    west_wbse
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)
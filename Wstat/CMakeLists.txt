# ==============================================================================
# WEST - Wstat/CMakeLists.txt
# ==============================================================================
#setting the sources for Wstat
set(WSTAT_LIB_SOURCES
    apply_operator.f90
    davidson_diago.f90
    davidson_restart.f90
    wstat_memory_report.f90
    wstat_setup.f90
    wstat_tools.f90
)

#library of Wstat since it is used by Wbse
add_library(west_wstat STATIC ${WSTAT_LIB_SOURCES})
set_target_properties(west_wstat PROPERTIES Fortran_PREPROCESS ON OUTPUT_NAME "wstat")
target_include_directories(west_wstat PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(west_wstat PUBLIC 
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)

# executable wstat.x 
add_executable(wstat wstat.f90)
set_target_properties(wstat PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "wstat"
    SUFFIX ".x"
)
#target_include_directories
target_link_libraries(wstat PRIVATE 
    west_wstat
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)

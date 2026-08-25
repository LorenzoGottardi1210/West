# ==============================================================================
# WEST - Wbse/CMakeLists.txt
# ==============================================================================
#setting the sources for Wbse 
set(WBSE_LIB_SOURCES
    bse_kernel.f90
    build_rhs_zvector_eq.f90
    calc_tau.f90
    check_ovl_wfc.f90
    hybrid_kernel.f90
    lanczos_restart.f90
    solve_zvector_eq_cg.f90
    td_liouville_oper.f90
    wbse_bgrp.f90
    wbse_calc_dens.f90
    wbse_calc_forces.f90
    wbse_davidson_diago.f90
    wbse_dot.f90
    wbse_dv.f90
    wbse_init_setup.f90
    wbse_lanczos_diago.f90
    wbse_localization.f90
    wbse_memory_report.f90
    wbse_setup.f90
    wbse_solve_e_psi.f90
    wbse_tools.f90
)
#library of Wbse since it shared between the two executables wbse.x and wbse_init.x 
add_library(west_wbse STATIC ${WBSE_LIB_SOURCES})
set_target_properties(west_wbse PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "wbse"
)
target_include_directories(west_wbse PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(west_wbse PUBLIC 
    west_wstat
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
    qe_lr_modules
)

# linking of wstat inside wbse library, as it is used by the two executables wbse.x and wbse_init.x 
target_link_libraries(west_wbse PUBLIC 
    west_wstat
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
    qe_lr_modules
)

# executable wbse.x 
add_executable(wbse wbse.f90)
set_target_properties(wbse PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "wbse"
    SUFFIX ".x"
)
target_link_libraries(wbse PRIVATE 
    west_wbse
    west_wstat
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
    qe_lr_modules
)

# executable wbse_init.x 
add_executable(wbse_init wbse_init.f90)
set_target_properties(wbse_init PROPERTIES 
    Fortran_PREPROCESS ON
    OUTPUT_NAME "wbse_init"
    SUFFIX ".x"
)
target_link_libraries(wbse_init PRIVATE 
    west_wbse
    west_wstat
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
    qe_lr_modules
)
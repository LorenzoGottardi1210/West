# ==============================================================================
# WEST - Wfreq/CMakeLists.txt
# ==============================================================================

# setting the sources for Wfreq 
set(WFREQ_LIB_SOURCES
    calc_corr.f90
    calc_exx2.f90
    calc_vxc.f90
    chi_invert.f90
    diago_lanczos.f90
    get_brak_hyper_parallel.f90
    set_freqlists.f90
    solve_deflated_lanczos_w_full_ortho.f90
    solve_eri.f90
    solve_gfreq.f90
    solve_h1e.f90
    solve_hf.f90
    solve_qp.f90
    solve_wfreq.f90
    wfreq_memory_report.f90
    wfreq_setup.f90
)
#adding the library of Wfreq since it is used by the executable wfreq.x 
add_library(west_wfreq STATIC ${WFREQ_LIB_SOURCES})
set_target_properties(west_wfreq PROPERTIES Fortran_PREPROCESS ON OUTPUT_NAME "wfreq")
target_include_directories(west_wfreq PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(west_wfreq PUBLIC 
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)

#Creating the executable wfreq.x 
add_executable(wfreq wfreq.f90)
set_target_properties(wfreq PROPERTIES Fortran_PREPROCESS ON OUTPUT_NAME "wfreq" SUFFIX ".x")

#linking of the target with the libraries of WEST, QE already inside the WEST build tree 
target_link_libraries(wfreq PRIVATE 
    west_wfreq
    west_io_kernel
    west_dfpt_kernel
    west_hamiltonian_kernel
    west_coulomb_kernel
)
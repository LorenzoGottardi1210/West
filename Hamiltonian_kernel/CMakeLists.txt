# ==============================================================================
# WEST - Hamiltonian_kernel/CMakeLists.txt
# ==============================================================================

# 1. creating target
add_library(west_hamiltonian_kernel STATIC
    apply_alpha_pa_to_m_wfcs.f90
    apply_alpha_pc_to_m_wfcs.f90
    apply_alpha_pv_to_m_wfcs.f90
    apply_hqp_to_m_wfcs.f90
    commut_Hx_psi.f90
    glbrak.f90
    k_psi.f90
)

# 2. apply the Fortran preprocessor to the target
set_target_properties(west_hamiltonian_kernel PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting the include directories for the target 
target_include_directories(west_hamiltonian_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. Linking of the target with the libraries of WEST, QE already inside the WEST build tree 
target_link_libraries(west_hamiltonian_kernel PUBLIC  west_modules west_tools west_fft_kernel)

# ==============================================================================
# WEST - DFPT_kernel/CMakeLists.txt
# ==============================================================================

# 1. Creating the target
add_library(west_dfpt_kernel STATIC
    apply_sternheimerop_to_m_wfcs.f90
    dfpt_module.f90
    linsolve_commut_Hx.f90
    linsolve_sternheimer_m_wfcts.f90
    precondition_m_wfcts.f90
)

# 2. apply the Fortran preprocessor to the target
set_target_properties(west_dfpt_kernel PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting the include directories for the target
target_include_directories(west_dfpt_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. Linking of the target with the libraries of WEST, QE already inside the WEST build tree
target_link_libraries(west_dfpt_kernel PUBLIC 
    west_modules
    west_tools
    west_fft_kernel
    west_para_kernel
)
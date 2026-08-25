# ==============================================================================
# WEST - Para_kernel/CMakeLists.txt
# ==============================================================================

# 1. Creating the target 
add_library(west_para_kernel STATIC
    class_idistribute.f90
    distribution_center.f90
    parallel_distributed_diago.f90
    west_mp.f90
)

# 2. apply the Fortran preprocessor to the target
set_target_properties(west_para_kernel PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting the include directories for the target 
target_include_directories(west_para_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. Linking of the target with the libraries of WEST, QE already inside the WEST build tree 
target_link_libraries(west_para_kernel PUBLIC west_tools west_modules)
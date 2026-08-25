#==============================================================================
# WEST - Coulomb_kernel/CMakeLists.txt
# ==============================================================================

# 1. creating the target 
add_library(west_coulomb_kernel STATIC
    class_coulomb.f90
    types_coulomb.f90
)

# 2. apply the Fortran preprocessor to the target
set_target_properties(west_coulomb_kernel PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting the include directories for the target
target_include_directories(west_coulomb_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. Linking of the target with the libraries of WEST, QE already inside the WEST build tree
target_link_libraries(west_coulomb_kernel PUBLIC west_tools west_modules)
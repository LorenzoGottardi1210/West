# ==============================================================================
#  CMakeLists.txt for Forpy Library 
# ==============================================================================

# adding static library forpy_mod.f90 
add_library(forpy STATIC forpy_mod.f90)

# adding the current binary directory to the include directories of forpy
target_include_directories(forpy PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# setting the properties of the forpy library to use free format and enable preprocessing 
set_target_properties(forpy PROPERTIES Fortran_FORMAT FREE)
set_target_properties(forpy PROPERTIES Fortran_PREPROCESS ON)
# linking the libraries found in the Master file, anyone in WEST linking "forpy" will automatically inherit 
# the link to Python C-API.

if(Python3_FOUND)
    target_link_libraries(forpy PUBLIC ${Python3_LIBRARIES})
    target_include_directories(forpy PUBLIC ${Python3_INCLUDE_DIRS})
endif()
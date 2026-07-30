# ==============================================================================
#  CMakeLists.txt for Json Library 
# ==============================================================================

add_library(json STATIC 
    json_kinds.f90
    json_parameters.f90
    json_string_utilities.f90
    json_value_module.f90
    json_file_module.f90
    json_module.f90
)

# setting the properties of the json library to use free format and enable preprocessing
set_target_properties(json PROPERTIES Fortran_FORMAT FREE Fortran_PREPROCESS ON )

target_compile_options(json PRIVATE -O0)
target_include_directories(json PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
# ==============================================================================
#  CMakeLists.txt for Base64 Library 
# ==============================================================================

add_library(base64 STATIC 
    cbase64.c 
    cbase64.h
    base64module.f90
)

# setting the properties of the base64 library to use free format and enable preprocessing
set_target_properties(base64 PROPERTIES Fortran_FORMAT FREE Fortran_PREPROCESS ON )


target_include_directories(base64 PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(base64 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
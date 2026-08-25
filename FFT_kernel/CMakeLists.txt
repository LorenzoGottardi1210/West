# ==============================================================================
# WEST - FFT_kernel/CMakeLists.txt
# ==============================================================================

add_library(west_fft_kernel STATIC
    fft_at_gamma.f90
    fft_at_k.f90
    fft_interpolation.f90
)

# Applichiamo il preprocessing a tutto il target in una riga sola
set_target_properties(west_fft_kernel PROPERTIES Fortran_PREPROCESS ON)

# Linking e include
target_link_libraries(west_fft_kernel PUBLIC west_modules)
target_include_directories(west_fft_kernel PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
# ==============================================================================
# WEST - Tools/CMakeLists.txt
# ==============================================================================

# 1. creating the target
add_library(west_tools STATIC
    bar.f90
    class_bz_grid.f90
    clean_scratchfiles.f90
    conversions.f90
    destroy_pw_arrays.f90
    do_setup.f90
    exx_go.f90
    exx_ungo.f90
    fetch_input.f90
    get_alpha_pv.f90
    heapsort.f90
    human_readable_time.f90
    init_pw_arrays.f90
    io_push.f90
    linear_algebra_kernel.f90
    my_mkdir.f90
    parse_command_arguments.f90
    pw_memory_report.f90
    set_dirs.f90
    set_eprec.f90
    set_nbndocc.f90
    set_npwq.f90
    types_bz_grid.f90
    wannier.f90
    west_environment.f90
    west_readin.f90
    west_print_clocks.f90
)

# 2. setting the target properties 
set_target_properties(west_tools PROPERTIES Fortran_PREPROCESS ON)

# 3. exporting of the modules for the next steps 
target_include_directories(west_tools PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# 4. linking of the libraries needed by west_tools, including the Quantum ESPRESSO libraries 
target_link_libraries(west_tools PUBLIC west_modules base64 json qe_xclib)
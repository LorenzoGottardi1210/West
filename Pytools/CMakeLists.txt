# ==============================================================================
#  CMakeLists.txt for Pytools  
# ==============================================================================

set(PYTOOLS_SRC
    read_json.py
    setup.py
    west_clientserver.py
    west_fetch_input.py
    west_function3d.py
    west_read_bse.py
    west_utils.py
)
set(PYTOOLS_ABS_SRC "")
foreach(src ${PYTOOLS_SRC})
    list(APPEND PYTOOLS_ABS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/${src})
endforeach()
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/install_timestamp
    # Il tuo $(PYT) -m pip install --user .
    COMMAND ${Python3_EXECUTABLE} -m pip install  .

    # Il tuo @touch $@ (usando la versione cross-platform integrata in CMake)
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/install_timestamp
    
    # Se uno di questi file cambia, CMake sa che deve rifare l'installazione
    DEPENDS ${PYTOOLS_ABS_SRC}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "installing pytools "
)
add_custom_target(west_pytools ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/install_timestamp
)
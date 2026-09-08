# ==============================================================================
# WEST - test-suite/CMakeLists.txt
# ==============================================================================

# directory for the test-suite build 
enable_testing()
#finding the directory of the test-suite for adding rapidly next times 
file(GLOB TEST_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "test*")
list(SORT TEST_DIRS)

# ATOMPAW - prevent build in-place

if (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})

  message(FATAL_ERROR "In-source builds not allowed. Please run CMake from top directory and specify a build directory (e.g., cmake -H. -Bbuild).")

endif()

# ATOMPAW - Execute bootstrap.sh script if configure does not exist

if (NOT EXISTS "${CMAKE_SOURCE_DIR}/configure")
  execute_process(
    VERBATIM
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR} ./bootstrap.sh 2>&1 >/dev/null
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()

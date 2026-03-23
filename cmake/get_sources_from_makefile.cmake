# ATOMPAW - Get source file list from Makefile.am
#  FUNCTION get_sources_from_makefile
#  Arguments:
#   - MAKEFILE_AM_PATH = path to Makefile.am file
#   - ATP_SRCS_STRING = string containing source files in Makefile.am
#   - SOURCE_FILES = output, list of source files

function(get_sources_from_makefile MAKEFILE_AM_PATH ATP_SRCS_STRING SOURCE_FILES)

  file(STRINGS "${MAKEFILE_AM_PATH}" LINES)

  set(sources)

  foreach(line ${LINES})
    if (line MATCHES "^[ \t]*${ATP_SRCS_STRING}[ \t]*=")
      string(REGEX MATCHALL "[a-zA-Z0-9_]+\\.(c|F90|f90)" files ${line})
      list(APPEND sources ${files})
    endif()
  endforeach()

  list(REMOVE_DUPLICATES sources)

  set(${SOURCE_FILES} ${sources} PARENT_SCOPE)

endfunction()

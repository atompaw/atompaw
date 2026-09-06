# ATOMPAW - Create dist target (to build tarball)

message(CHECK_START "Creating dist target ")
set(CMAKE_REQUIRED_QUIET TRUE)

# Temporary destination directory
set(DIST_DIR atompaw-${PROJECT_VERSION})

# Files to distribute in the source package
set(DIST_PATHS
    src
    doc
    config/m4
    example
    CHANGELOG
    Makefile.am
    README.md
    AUTHORS
    configure.ac
    COPYING
    INSTALL
    README
    VERSION
    cmake
    CMakeLists.txt
    configure
    config.h.in
    config/gnu
    aclocal.m4
    Makefile.in
)

# Patterns that must never be shipped, in any directory
set(DIST_EXCLUDE_REGEX
"(~|\\.bak|\\.old|\\.orig|\\.rej|\\.swp)$|(^|/)(tmp|tmpsave|core)[^/]*$|(^|/)#[^/]*#$|_old\\.[Ff]90$|(^|/)\\.[^/]+(/|$)")

# For the src/ directory, only real source and build files are distributed
# (Fortran/C sources, configure templates, Makefiles, CMake files). Anything
# else living in src/ is left out of the tarball
set(SRC_INCLUDE_PATTERNS
    "*.F90" "*.f90" "*.c" "*.h" "*.in"          # sources + templates configure
    "Makefile.am" "Makefile.in" "CMakeLists.txt")

# For the example/ directory, drop the run-time output files and directories,
# keeping only the inputs and reference results worth shipping
set(EXAMPLE_EXCLUDE_REGEX
"(^|/)(explore|dummy|NC|OCCWFN|rvf|rVx|vloc|checkvxc|hattest|tp|compare\\.abinit)(/|$)|(^|/)(den|pot)[^/]*$|(^|/)([0-9]+\\.)?wfn00[^/]*$|(^|/)[^/]*(AE0|SC1)$")

# Extract a full file list from previous list
set(DIST_FILES)
foreach(path ${DIST_PATHS})
  if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${path})
    if (path STREQUAL "src")
      # Whitelisted, non-recursive glob: keep only wanted file types
      set(_globs)
      foreach(pat ${SRC_INCLUDE_PATTERNS})
        list(APPEND _globs "${CMAKE_SOURCE_DIR}/${path}/${pat}")
      endforeach()
      file(GLOB files RELATIVE "${CMAKE_SOURCE_DIR}/${path}" ${_globs})
    else()
      # Everything else: full recursive content of the directory
      file(GLOB_RECURSE files
           RELATIVE "${CMAKE_SOURCE_DIR}/${path}"
           "${CMAKE_SOURCE_DIR}/${path}/*")
    endif()
    foreach(file ${files})
      # Drop backup / scratch files wherever they are
      if ("${file}" MATCHES "${DIST_EXCLUDE_REGEX}")
        continue()
      endif()
      # Drop run-time output from the example/ tree
      if (path STREQUAL "example" AND "${file}" MATCHES "${EXAMPLE_EXCLUDE_REGEX}")
        continue()
      endif()
      list(APPEND DIST_FILES "${path}/${file}")
    endforeach()
  else()
    list(APPEND DIST_FILES ${path})
  endif()
endforeach()
list(REMOVE_DUPLICATES DIST_FILES)

# Create a temporary command to copy files
set(COPY_COMMANDS)
foreach(file ${DIST_FILES})
    list(APPEND COPY_COMMANDS
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_SOURCE_DIR}/${file}
        ${CMAKE_CURRENT_BINARY_DIR}/${DIST_DIR}/${file}
    )
endforeach()

# Add dist target (with possible execution of bootstrap/wipeout scripts)
add_custom_target(dist
    COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DSTAMP_FILE=${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp
            -P ${CMAKE_SOURCE_DIR}/cmake/run_bootstrap.cmake
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DIST_DIR}
    ${COPY_COMMANDS}
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR} tar czf ${DIST_DIR}.tar.gz ${DIST_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR}
    COMMAND ${CMAKE_COMMAND}
            -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
            -DSTAMP_FILE=${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp
            -P ${CMAKE_SOURCE_DIR}/cmake/run_wipeout.cmake
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp
               ${CMAKE_CURRENT_BINARY_DIR}/bootstrap_tmp.stamp
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Creating ${DIST_DIR}.tar.gz tarball"
    VERBATIM
)

# add_custom_target(dist
#     COMMAND test -f ${CMAKE_SOURCE_DIR}/configure && exit 0 ||
#     ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp && 
#     ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR} ./bootstrap.sh >/dev/null 2>&1
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR}
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${DIST_DIR}
#     ${COPY_COMMANDS}
#     COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR} tar czf ${DIST_DIR}.tar.gz ${DIST_DIR}
# #   COMMAND ${CMAKE_COMMAND} -E chdir ${DIST_DIR} tar czf ../${DIST_DIR}.tar.gz .
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR}
#     COMMAND test -f ${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp && 
#            ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR} ./wipeout.sh >/dev/null 2>&1 &&
#            ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/configure~ &&
#            ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp || 
#            exit 0
#     BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/bootstrap.stamp
#                ${CMAKE_CURRENT_BINARY_DIR}/bootstrap_tmp.stamp
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     COMMENT "Creating ${DIST_DIR}.tar.gz tarball"
#     VERBATIM
# )

unset(CMAKE_REQUIRED_QUIET)
message(CHECK_PASS "OK")

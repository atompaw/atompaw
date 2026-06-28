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
    ChangeLog
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

# Extract a full file list from previous list
set(DIST_FILES)
foreach(path ${DIST_PATHS})
  if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${path})
    file(GLOB_RECURSE files 
         RELATIVE "${CMAKE_SOURCE_DIR}/${path}"
         "${CMAKE_SOURCE_DIR}/${path}/*"
        )
    foreach(file ${files})
      list(APPEND DIST_FILES "${path}/${file}")
    endforeach()
  else()
    list(APPEND DIST_FILES ${path})
  endif()
endforeach()

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

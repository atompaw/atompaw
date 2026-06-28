# ATOMPAW - Execute autotools wipeout script
# Expecte variables : SOURCE_DIR, STAMP_FILE

# Nothing to do if configure already exists
if (NOT EXISTS "${STAMP_FILE}")
  return()
endif()

# Call wipeout script
execute_process(
  COMMAND ./wipeout.sh
  WORKING_DIRECTORY "${SOURCE_DIR}"
  OUTPUT_QUIET
  ERROR_QUIET
)

file(REMOVE "${SOURCE_DIR}/configure~")
file(REMOVE "${STAMP_FILE}")

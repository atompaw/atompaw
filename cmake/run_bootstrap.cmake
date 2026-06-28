# ATOMPAW - Execute autotools bootstrap script
# Expected variables : SOURCE_DIR, STAMP_FILE

# Nothing to do if configure already exists
if(EXISTS "${SOURCE_DIR}/configure")
    return()  # rien à faire, ./configure existe déjà
endif()

# Call bootstrap script
file(TOUCH "${STAMP_FILE}")
execute_process(
    COMMAND ./bootstrap.sh
    WORKING_DIRECTORY "${SOURCE_DIR}"
    OUTPUT_QUIET
    ERROR_QUIET
)
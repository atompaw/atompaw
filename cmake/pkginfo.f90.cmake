MODULE PkgInfo

  IMPLICIT NONE

  PUBLIC

  CHARACTER(len=*), PARAMETER :: atp_package = "@PROJECT_NAME@"
  CHARACTER(len=*), PARAMETER :: atp_version = "@PROJECT_VERSION@"
  CHARACTER(len=*), PARAMETER :: atp_host    = "@CMAKE_HOST_SYSTEM@"
  CHARACTER(len=*), PARAMETER :: atp_target  = "@CMAKE_SYSTEM@"

END MODULE PkgInfo

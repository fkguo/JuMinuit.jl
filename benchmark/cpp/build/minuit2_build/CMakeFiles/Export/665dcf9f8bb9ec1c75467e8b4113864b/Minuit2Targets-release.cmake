#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Minuit2::Minuit2Math" for configuration "Release"
set_property(TARGET Minuit2::Minuit2Math APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Minuit2::Minuit2Math PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libMinuit2Math.a"
  )

list(APPEND _cmake_import_check_targets Minuit2::Minuit2Math )
list(APPEND _cmake_import_check_files_for_Minuit2::Minuit2Math "${_IMPORT_PREFIX}/lib/libMinuit2Math.a" )

# Import target "Minuit2::Minuit2" for configuration "Release"
set_property(TARGET Minuit2::Minuit2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Minuit2::Minuit2 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libMinuit2.a"
  )

list(APPEND _cmake_import_check_targets Minuit2::Minuit2 )
list(APPEND _cmake_import_check_files_for_Minuit2::Minuit2 "${_IMPORT_PREFIX}/lib/libMinuit2.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

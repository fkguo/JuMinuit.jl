# Install script for directory: /Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/src/math

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/cpp/build/minuit2_build/src/math/libMinuit2Math.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2Math.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2Math.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2Math.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Minuit2/Fit" TYPE FILE FILES "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Fit/ParameterSettings.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Minuit2/Math" TYPE FILE FILES
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/Error.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/FitMethodFunction.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/Functor.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/GenAlgoOptions.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/IFunction.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/IFunctionfwd.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/IOptions.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/Minimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/MinimizerOptions.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/Util.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/WrappedFunction.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Math/WrappedParamFunction.h"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/cpp/build/minuit2_build/src/math/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()

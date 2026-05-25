# Install script for directory: /Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/cpp/build/minuit2_build/src/libMinuit2.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMinuit2.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Minuit2/Minuit2" TYPE FILE FILES
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ABObj.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ABProd.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ABSum.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ABTypes.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/AnalyticalGradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BFGSErrorUpdator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicFunctionGradient.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicFunctionMinimum.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicMinimumError.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicMinimumParameters.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicMinimumSeed.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/BasicMinimumState.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/CombinedMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/CombinedMinimumBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ContoursError.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/DavidonErrorUpdator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FCNAdapter.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FCNBase.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FCNGradAdapter.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FCNGradientBase.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliChi2FCN.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliErrorUpdator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliFCNAdapter.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliFCNBase.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliGradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliMaximumLikelihoodFCN.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliStandardChi2FCN.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FumiliStandardMaximumLikelihoodFCN.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FunctionGradient.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FunctionMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/FunctionMinimum.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/GenericFunction.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/GradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/HessianGradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/InitialGradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LASymMatrix.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LAVector.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LaInverse.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LaOuterProduct.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LaProd.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/LaSum.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MPIProcess.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MatrixInverse.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumError.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumErrorUpdator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumParameters.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumSeed.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumSeedGenerator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinimumState.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinosError.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/Minuit2Minimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MinuitParameter.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnApplication.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnConfig.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnContours.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnCovarianceSqueeze.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnCross.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnEigen.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnFcn.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnFumiliMinimize.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnFunctionCross.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnGlobalCorrelationCoeff.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnHesse.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnLineSearch.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnMachinePrecision.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnMatrix.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnMigrad.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnMinimize.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnMinos.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnParabola.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnParabolaFactory.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnParabolaPoint.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnParameterScan.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnPlot.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnPosDef.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnPrint.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnScan.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnSeedGenerator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnSimplex.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnStrategy.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnTiny.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnTraceObject.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnUserCovariance.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnUserFcn.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnUserParameterState.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnUserParameters.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnUserTransformation.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/MnVectorTransform.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ModularFunctionMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/NegativeG2LineSearch.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/Numerical2PGradientCalculator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ParametricFunction.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ScanBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/ScanMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SimplexBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SimplexMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SimplexParameters.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SimplexSeedGenerator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SinParameterTransformation.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SqrtLowParameterTransformation.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/SqrtUpParameterTransformation.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/StackAllocator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/VariableMetricBuilder.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/VariableMetricEDMEstimator.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/VariableMetricMinimizer.h"
    "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/inc/Minuit2/VectorOuterProduct.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/cpp/build/minuit2_build/src/math/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/cpp/build/minuit2_build/src/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()

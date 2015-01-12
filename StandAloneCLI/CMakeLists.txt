PROJECT(ABC_CLI)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS} )
INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})  

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/basicimg
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/bias
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/brainseg
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/common
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/register
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/robust
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/spr
  ${CMAKE_CURRENT_SOURCE_DIR}/../Engine/xmlio
)

SET(ABC_CLI_SRCS
  ../Engine/brainseg/EMSParameters.cxx
  ../Engine/brainseg/EMSegmentationFilter_float+float.cxx
  ../Engine/brainseg/runEMS.cxx
  ../Engine/common/Log.cxx
  ../Engine/common/MersenneTwisterRNG.cxx
  ../Engine/common/Timer.cxx
  ../Engine/common/muFile.cxx
  ../Engine/register/AmoebaOptimizer.cxx
  ../Engine/register/AtlasRegistrationMethod_float+float.cxx
  ../Engine/register/ChainedAffineTransform3D.cxx
  ../Engine/register/GradientDescentOptimizer.cxx
  ../Engine/register/PairRegistrationMethod_float.cxx
  ../Engine/register/PowellOptimizer.cxx
  ../Engine/register/SimulatedAnnealingOptimizer.cxx
  ../Engine/robust/FastMCDSampleFilter.cxx
  ../Engine/robust/KruskalMSTClusteringProcess.cxx
  ../Engine/spr/KMeansEstimator.cxx
  ../Engine/xmlio/EMSParametersXMLFile.cxx
  main_cli.cxx
)

ADD_EXECUTABLE(ABC_MINC ${ABC_CLI_SRCS})

TARGET_LINK_LIBRARIES(ABC_MINC ${ITK_LIBRARIES})

INSTALL(TARGETS
          ABC_MINC
        RUNTIME DESTINATION bin )

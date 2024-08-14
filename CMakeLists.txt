cmake_minimum_required( VERSION 2.8.12 )

project(plotIt LANGUAGES CXX)

include(GNUInstallDirs)
include(ExternalProject)

if(EXISTS "$ENV{CMAKE_PREFIX_PATH}/include/boost")
  message(STATUS "Will use $ENV{CMAKE_PREFIX_PATH} as base path for boost")
  set(BOOST_ROOT $ENV{CMAKE_PREFIX_PATH})
  set(Boost_NO_BOOST_CMAKE ON)
endif()

find_package(ROOT REQUIRED COMPONENTS HistPainter Tree)
find_package(Boost REQUIRED COMPONENTS filesystem regex system)

ExternalProject_Add(
  yaml-cpp-build
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.tar.gz
  CMAKE_ARGS -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/external -DCMAKE_CXX_FLAGS='-Wno-error=maybe-uninitialized'
  )
add_library(yaml-cpp STATIC IMPORTED)
add_dependencies(yaml-cpp yaml-cpp-build)
set_target_properties(yaml-cpp PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/external/lib/libyaml-cpp.a)
ExternalProject_Add(
  tclap
  URL https://github.com/mirror/tclap/archive/refs/tags/v1.2.4.tar.gz
  BUILD_IN_SOURCE 1
  PATCH_COMMAND ./autotools.sh
  CONFIGURE_COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/external
  INSTALL_COMMAND make -ki install
  )

set(SRCS
  src/plotIt.cc
  src/summary.cc
  src/systematics.cc
  src/TH1Plotter.cc
  src/types.cc
  src/utilities.cc
  src/uuid.cc
  )

add_executable(plotIt ${SRCS})
add_dependencies(plotIt tclap)
# workaround, should be inherited from ROOT dependency targets (if present), but is not specified there for versions below 6.18.00
if((${ROOT_VERSION} VERSION_LESS "6.18.00"))
  if(${ROOT_cxx17_FOUND})
    target_compile_features(plotIt PRIVATE cxx_std_17)
  elseif(${ROOT_cxx14_FOUND})
    target_compile_features(plotIt PRIVATE cxx_std_14)
  elseif(${ROOT_cxx11_FOUND})
    target_compile_features(plotIt PRIVATE cxx_std_11)
  endif()
endif()
if(TARGET ROOT::Tree AND TARGET ROOT::HistPainter)
  target_link_libraries(plotIt ROOT::HistPainter ROOT::Tree dl Boost::filesystem Boost::regex Boost::system yaml-cpp)
  target_include_directories(plotIt PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> ${CMAKE_CURRENT_BINARY_DIR}/external/include)
else()
  target_link_libraries(plotIt ${ROOT_LIBRARIES} dl Boost::filesystem Boost::regex Boost::system yaml-cpp)
  target_include_directories(plotIt PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> ${CMAKE_CURRENT_BINARY_DIR}/external/include ${ROOT_INCLUDE_DIRS})
endif()
install(TARGETS plotIt
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

cmake_minimum_required(VERSION 3.10)

# Settings
project(LIE_GROUP_CODE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -O3)



## Packages Import
##########################################################
# CODAC-LIB (and IBEX-LIB) package
# Adding IBEX
find_package(IBEX REQUIRED)
ibex_init_common() # IBEX should have installed this function
message(STATUS "Found IBEX version ${IBEX_VERSION}")

# Adding Codac
find_package(CODAC REQUIRED)
message(STATUS "Found Codac version ${CODAC_VERSION}")

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})

###########################################################
# CAPD package
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_CAPD REQUIRED  capd capd-gui mpcapd mpcapd-gui capdRedHom )  # this looks for .pc file
message(STATUS "Found CAPD")
include_directories(${PKG_CAPD_INCLUDE_DIRS})

###########################################################
# Ipegenerator
find_package(ipegenerator REQUIRED)
message(STATUS "Found ipegenerator")
include_directories(${IPEGENERATOR_INCLUDE_DIRS})


## Source files and headers definition
###########################################################
# src bundles for executables

set (SRC_EXAMPLE1_ARTICLE "src/optimum_speed/article_figures/example_1_article_form.cpp")
set (SRC_EXAMPLE2_ARTICLE "src/optimum_speed/article_figures/example_2_article_form.cpp")
set (SRC_EXAMPLE3_ARTICLE "src/optimum_speed/article_figures/example_3_article_form.cpp")
set (SRC_EXAMPLE4_ARTICLE "src/optimum_speed/article_figures/example_4_article_form.cpp")






set(LIBS ${IPEGENERATOR_LIBRARIES} ${CODAC_LIBRARIES} Ibex::ibex ${CODAC_LIBRARIES} ${PKG_CAPD_LDFLAGS})

## Executables/Libs creation
###########################################################


set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/../lib)
include_directories("src/tools")


set(EXECUTABLE_OUTPUT_PATH bin)

set(LIBS   tools ${LIBS} )


add_subdirectory(src)

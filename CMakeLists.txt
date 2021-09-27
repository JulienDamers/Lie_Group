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
include_directories(${IPEGENERATOR_INCLUDE_DIRS})


## Source files and headers definition
###########################################################
# src bundles for executables
set (SRC_EXAMPLE1_VIBES "src/optimum_speed/example_1.cpp")
set (SRC_EXAMPLE2_VIBES "src/optimum_speed/example_2.cpp")
set (SRC_EXAMPLE3_VIBES "src/optimum_speed/example_3.cpp")
set (SRC_EXAMPLE4_VIBES "src/optimum_speed/example_4.cpp")

set (SRC_EXAMPLE1_ARTICLE "src/optimum_speed/example_1_article_form.cpp")
set (SRC_EXAMPLE2_ARTICLE "src/optimum_speed/example_2_article_form.cpp")
set (SRC_EXAMPLE3_ARTICLE "src/optimum_speed/example_3_article_form.cpp")
set (SRC_EXAMPLE4_ARTICLE "src/optimum_speed/example_4_article_form.cpp")

set (SRC_EXAMPLE1_CN_ARTICLE "src/contractor_network_version/ex_1_cn_method.cpp")
set (SRC_EXAMPLE2_CN_ARTICLE "src/contractor_network_version/ex_2_cn_method.cpp")
set (SRC_EXAMPLE3_CN_ARTICLE "src/contractor_network_version/ex_3_cn_method.cpp")
set (SRC_EXAMPLE4_CN_ARTICLE "src/contractor_network_version/ex_4_cn_method.cpp")

list(APPEND SRC_TOOLS_LIB
        "src/tools/tools.h"
        "src/tools/lie_group_ex3_separator.h"
        "src/tools/lie_group_ex3_separator.cpp"
        "src/tools/lie_group_ex4_separator.h"
        "src/tools/lie_group_ex4_separator.cpp"
        "src/tools/pyibex_SepProj.h"
        "src/tools/pyibex_SepProj.cpp"
        "src/tools/ctc_cn.h"
        "src/tools/ctc_cn.cpp"
        "src/tools/sivia.h"
        "src/tools/sivia.cpp"
        )



set(LIBS ${IPEGENERATOR_LIBRARIES} ${CODAC_LIBRARIES} Ibex::ibex ${CODAC_LIBRARIES} ${PKG_CAPD_LDFLAGS})

## Executables/Libs creation
###########################################################


set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/../lib)
include_directories("src/tools")
add_library(tools ${SRC_TOOLS_LIB})
target_include_directories(tools PUBLIC "src/tools" ${CODAC_INCLUDE_DIRS})
target_link_libraries(tools PUBLIC ${CODAC_LIBRARIES} Ibex::ibex ${CODAC_LIBRARIES} ${IPEGENERATOR_LIBRARIES})


set(EXECUTABLE_OUTPUT_PATH bin)

set(LIBS   tools ${LIBS} )

########################################################################################################################
# Vibes version
########################################################################################################################

# optimum speed mode

add_executable(example_1_vibes ${SRC_EXAMPLE1_VIBES})
target_compile_options(example_1_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_vibes PUBLIC ${LIBS})

add_executable(example_2_vibes ${SRC_EXAMPLE2_VIBES})
target_compile_options(example_2_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_vibes ${LIBS})

add_executable(example_3_vibes ${SRC_EXAMPLE3_VIBES})
target_compile_options(example_3_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_vibes ${LIBS})

add_executable(example_4_vibes ${SRC_EXAMPLE4_VIBES})
target_compile_options(example_4_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_vibes ${LIBS})

# CN mode

add_executable(example_1_cn_article ${SRC_EXAMPLE1_CN_ARTICLE})
target_compile_options(example_1_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_cn_article PUBLIC ${LIBS})

add_executable(example_2_cn_article ${SRC_EXAMPLE2_CN_ARTICLE})
target_compile_options(example_2_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_cn_article PUBLIC ${LIBS})

add_executable(example_3_cn_article ${SRC_EXAMPLE3_CN_ARTICLE})
target_compile_options(example_3_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_cn_article PUBLIC ${LIBS})

add_executable(example_4_cn_article ${SRC_EXAMPLE4_CN_ARTICLE})
target_compile_options(example_4_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_cn_article PUBLIC ${LIBS})

########################################################################################################################
# Ipegenerator version
########################################################################################################################


add_executable(example_1_article ${SRC_EXAMPLE1_ARTICLE})
target_compile_options(example_1_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_article PUBLIC ${LIBS})

add_executable(example_2_article ${SRC_EXAMPLE2_ARTICLE})
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_article PUBLIC ${LIBS})

add_executable(example_3_article ${SRC_EXAMPLE3_ARTICLE})
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_article PUBLIC ${LIBS})

add_executable(example_4_article ${SRC_EXAMPLE4_ARTICLE})
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_article PUBLIC ${LIBS})


########################################################################################################################
# Comparisons for reviews version
########################################################################################################################

add_executable(integration_ex1 src/comparisons_review_2/integration_ex_1.cpp)
target_compile_options(integration_ex1 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex1 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex1 PUBLIC ${LIBS})

add_executable(integration_ex2 src/comparisons_review_2/integration_ex_2.cpp)
target_compile_options(integration_ex2 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex2 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex2 PUBLIC ${LIBS})

add_executable(integration_ex3 src/comparisons_review_2/integration_ex_3.cpp)
target_compile_options(integration_ex3 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex3 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex3 PUBLIC ${LIBS})

add_executable(integration_ex4 src/comparisons_review_2/integration_ex_4.cpp)
target_compile_options(integration_ex4 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex4 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex4 PUBLIC ${LIBS})

add_executable(sivia_capd_ex2 src/comparisons_review_2/sivia_capd_integration.cpp)
target_compile_options(sivia_capd_ex2 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(sivia_capd_ex2 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(sivia_capd_ex2 PUBLIC ${LIBS})

add_executable(particles src/comparisons_review_2/particles_comparison.cpp)
target_compile_options(particles PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(particles PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(particles PUBLIC ${LIBS})

add_executable(reachability src/comparisons_review_2/reachability.cpp)
target_compile_options(reachability PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(reachability PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(reachability PUBLIC ${LIBS})
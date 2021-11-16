add_executable(integration_ex1 integration_ex_1.cpp)
target_compile_options(integration_ex1 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex1 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex1 PUBLIC ${LIBS})

add_executable(integration_ex2 integration_ex_2.cpp)
target_compile_options(integration_ex2 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex2 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex2 PUBLIC ${LIBS})

add_executable(integration_ex3 integration_ex_3.cpp)
target_compile_options(integration_ex3 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex3 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex3 PUBLIC ${LIBS})

add_executable(integration_ex4 integration_ex_4.cpp)
target_compile_options(integration_ex4 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(integration_ex4 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(integration_ex4 PUBLIC ${LIBS})

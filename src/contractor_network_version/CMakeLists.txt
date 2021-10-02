add_executable(example_1_cn_article ex_1_cn_method.cpp)
target_compile_options(example_1_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_cn_article PUBLIC ${LIBS})

add_executable(example_2_cn_article ex_2_cn_method.cpp)
target_compile_options(example_2_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_cn_article PUBLIC ${LIBS})

add_executable(example_3_cn_article ex_3_cn_method.cpp)
target_compile_options(example_3_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_cn_article PUBLIC ${LIBS})

add_executable(example_4_cn_article ex_4_cn_method.cpp)
target_compile_options(example_4_cn_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_cn_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_cn_article PUBLIC ${LIBS})
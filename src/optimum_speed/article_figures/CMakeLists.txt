add_executable(example_1_article example_1_article_form.cpp)
target_compile_options(example_1_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_article PUBLIC ${LIBS})

add_executable(example_2_article example_2_article_form.cpp)
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_article PUBLIC ${LIBS})

add_executable(example_3_article example_3_article_form.cpp)
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_article PUBLIC ${LIBS})

add_executable(example_4_article example_4_article_form.cpp)
target_compile_options(example_2_article PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_article PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_article PUBLIC ${LIBS})

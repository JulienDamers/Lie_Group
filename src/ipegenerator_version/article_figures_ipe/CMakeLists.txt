add_executable(example_1_ipe example_1_article_form.cpp)
target_compile_options(example_1_ipe PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_ipe PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_ipe PUBLIC ${LIBS})

add_executable(example_2_ipe example_2_article_form.cpp)
target_compile_options(example_2_ipe PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_ipe PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_ipe PUBLIC ${LIBS})

#add_executable(example_3_ipe example_3_article_form.cpp)
#target_compile_options(example_2_ipe PUBLIC ${CODAC_CXX_FLAGS})
#target_include_directories(example_3_ipe PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
#target_link_libraries(example_3_ipe PUBLIC ${LIBS})

add_executable(example_4_ipe example_4_article_form.cpp)
target_compile_options(example_2_ipe PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_ipe PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_ipe PUBLIC ${LIBS})

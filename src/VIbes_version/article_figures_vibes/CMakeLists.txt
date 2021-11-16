add_executable(example_1_vibes example_1.cpp)
target_compile_options(example_1_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_1_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_1_vibes PUBLIC ${LIBS})

add_executable(example_2_vibes example_2.cpp)
target_compile_options(example_2_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_2_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_2_vibes ${LIBS})

add_executable(example_3_vibes example_3.cpp)
target_compile_options(example_3_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_3_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_3_vibes ${LIBS})

add_executable(example_4_vibes example_4.cpp)
target_compile_options(example_4_vibes PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(example_4_vibes PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(example_4_vibes ${LIBS})
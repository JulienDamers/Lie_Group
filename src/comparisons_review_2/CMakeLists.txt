add_executable(sivia_capd_ex2 sivia_capd_integration.cpp)
target_compile_options(sivia_capd_ex2 PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(sivia_capd_ex2 PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(sivia_capd_ex2 PUBLIC ${LIBS})

add_executable(particles particles_comparison.cpp)
target_compile_options(particles PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(particles PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(particles PUBLIC ${LIBS})

add_executable(reachability reachability.cpp)
target_compile_options(reachability PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(reachability PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(reachability PUBLIC ${LIBS})


#add_executable(collins collins.cpp)
#target_compile_options(collins PUBLIC ${CODAC_CXX_FLAGS})
#target_include_directories(collins PUBLIC ${CODAC_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
#target_link_libraries(collins PUBLIC ${LIBS})




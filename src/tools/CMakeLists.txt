list(APPEND SRC_TOOLS_LIB
        "tools.h"
        "lie_group_ex3_separator.h"
        "lie_group_ex3_separator.cpp"
        "lie_group_ex4_separator.h"
        "lie_group_ex4_separator.cpp"
        "pyibex_SepProj.h"
        "pyibex_SepProj.cpp"
        "ctc_cn.h"
        "ctc_cn.cpp"
        "sivia.h"
        "sivia.cpp"
        )


add_library(tools ${SRC_TOOLS_LIB})
target_include_directories(tools PUBLIC "src/tools" ${CODAC_INCLUDE_DIRS})
target_link_libraries(tools PUBLIC ${CODAC_LIBRARIES} Ibex::ibex ${CODAC_LIBRARIES} ${IPEGENERATOR_LIBRARIES})
FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master) # or whatever tag you want
FetchContent_MakeAvailable(matplotplusplus)

file(GLOB
  SOURCE_FILES
  *.cpp
)

foreach(ITEM ${SOURCE_FILES})
    get_filename_component(EXECUTABLE ${ITEM} NAME_WE)
    add_executable(${EXECUTABLE}
        ${ITEM}
    )
    target_include_directories(${EXECUTABLE}
        PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/test/include"
    )
    target_link_libraries(${EXECUTABLE}
        PRIVATE
        matplot
        eigen
    )
endforeach()
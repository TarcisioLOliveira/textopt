set(PROJECT_NAME textopt)

set(SOURCES
    main.cpp
    MMASolver.cpp
    render.cpp
    opt_func.cpp
    texture.cpp
    analysis.cpp
    slp.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

set(REQ_LIBS sfml-graphics)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} ${REQ_LIBS} OpenMP::OpenMP_CXX)
else()
    target_link_libraries(${PROJECT_NAME} ${REQ_LIBS})
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)

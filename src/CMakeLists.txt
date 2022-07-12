set(PROJECT_NAME textopt)

set(SOURCES
    main.cpp
    MMASolver.cpp
    render.cpp
    opt_func.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} sfml-graphics)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)

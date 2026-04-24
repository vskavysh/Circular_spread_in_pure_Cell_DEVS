set(EXECUTABLE_NAME brain_wave_celldevs)

add_executable(${EXECUTABLE_NAME}
    main.cpp
)

if(DEFINED ENV{CADMIUM})
    target_include_directories(${EXECUTABLE_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        $ENV{CADMIUM}/include
    )
else()
    target_include_directories(${EXECUTABLE_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
endif()

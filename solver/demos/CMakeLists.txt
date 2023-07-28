set(demos arenstorf brownian brusselator brusselator_dirk exp heat pendulum lorenz lotka_volterra curtiss_hirschfelder transport)
set(demos_eigen brusselator_dirk)
set(demos_samurai transport)

foreach(demo IN LISTS demos)
    message(STATUS "demo=${demo}")
    add_executable(${demo} ${demo}.cpp)

    if(${demo} IN_LIST demos_eigen)
        find_package (Eigen3 3.4 REQUIRED NO_MODULE)
        target_link_libraries(${demo} Eigen3::Eigen)
    endif()
    if(${demo} IN_LIST demos_samurai)
        find_package (samurai REQUIRED)
        target_link_libraries(${demo} samurai::samurai)
    endif()

    target_link_libraries(${demo} ponio)
endforeach()

foreach(demo IN LISTS demos)
    add_custom_target(visu_${demo} python ${CMAKE_SOURCE_DIR}/solver/demos/${demo}.py)
    add_custom_target(${demo}_visu make visu_${demo})
endforeach()

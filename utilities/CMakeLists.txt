#file(GLOB children RELATIVE ${current} ${current}/*)
#foreach(subdir ${children})
#    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child} AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt)
#	    add_subdirectory(${subdir})
#    endif()
#endforeach()

################################################################################
# Utilities
################################################################################
if(VEGAFEM_BUILD_UTILITIES)
    set(VegaFEM_utilities
        computeDistanceField
        finiteDifferenceTest
        isosurfaceMesher
        objMergeFiles
        tetMesher
    )

    foreach(exe ${VegaFEM_utilities})
      file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.cpp")
      add_executable(${exe} ${sources})
      target_link_libraries(${exe} PRIVATE ${PROJECT_NAME})
      if($<TARGET_RUNTIME_DLLS:${exe}>)
        add_custom_command(TARGET ${exe} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${exe}> $<TARGET_FILE_DIR:${exe}> 
          COMMAND_EXPAND_LISTS
        )
        install(TARGETS $<IF:$<PLATFORM_ID:Windows>,${exe},>
            RUNTIME_DEPENDENCIES
            PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
            POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
            DIRECTORIES "$<TARGET_FILE_DIR:${exe}>"
        )
      endif()
      install(TARGETS ${exe} RUNTIME DESTINATION bin CONFIGURATIONS Release)
    endforeach()



    set(VegaFEM_utilities
        generateInterpolant
        generateInterpolationMatrix
        generateMassMatrix
        generateSurfaceMesh
    )
    foreach(target_name ${VegaFEM_utilities})
        add_executable(${target_name} "volumetricMeshUtilities/${target_name}.cpp")
        target_link_libraries(${target_name} PRIVATE ${PROJECT_NAME})
        if($<TARGET_RUNTIME_DLLS:${target_name}>)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${target_name}> $<TARGET_FILE_DIR:${target_name}> 
                COMMAND_EXPAND_LISTS
            )
            install(TARGETS ${target_name}
                RUNTIME_DEPENDENCIES
                PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
                POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
                DIRECTORIES "$<TARGET_FILE_DIR:${target_name}> "
            )
        endif()
        install(TARGETS ${target_name} RUNTIME DESTINATION bin CONFIGURATIONS Release)
    endforeach()


    if(VEGAFEM_BUILD_COPYLEFT)
    set(VegaFEM_utilities_gpl
        immersionMesher
        virtualTetsDriver
    )
    foreach(exe ${VegaFEM_utilities_gpl})
      file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.cpp")
      add_executable(${exe} ${sources})
      target_link_libraries(${exe} PRIVATE vegafem)
      if($<TARGET_RUNTIME_DLLS:${exe}>)
        add_custom_command(TARGET ${exe} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${exe}> $<TARGET_FILE_DIR:${exe}> 
          COMMAND_EXPAND_LISTS
        )
        install(TARGETS ${exe}
            RUNTIME_DEPENDENCIES
            PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
            POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
            DIRECTORIES "$<TARGET_FILE_DIR:${exe}> "
        )
      endif()
      install(TARGETS ${exe} RUNTIME DESTINATION bin CONFIGURATIONS Release)
    endforeach()
    endif()
endif()

################################################################################
# Utilities graphics
################################################################################
if(VEGAFEM_BUILD_UTILITIES_GRAPHICS)
    find_package(X11 REQUIRED) # freeglut on Linux
    find_package(glui CONFIG REQUIRED) # clothBW-rt displayObj editShapeARAP interactiveDeformableSimulator reducedDynamicSolver

    set(VegaFEM_utilities
        clothBW-rt
        displayObj
        editShapeARAP
        interactiveDeformableSimulator
        # largeModalDeformationFactory
        reducedDynamicSolver-rt
    )
    foreach(target_name ${VegaFEM_utilities})
      file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS "${target_name}/*.c" "${target_name}/*.cpp")
      add_executable(${target_name} ${sources})
      target_link_libraries(${target_name} PRIVATE ${PROJECT_NAME} 
          OpenGL::GL 
          GLUT::GLUT
          $<IF:$<TARGET_EXISTS:glui::glui>,glui::glui,glui::glui_static> #reducedDynamicSolver-rt interactiveDeformableSimulator editShapeARAP displayObj clothBW-rt
          $<$<PLATFORM_ID:Linux,Darwin>:X11::X11>
      )
      if($<TARGET_RUNTIME_DLLS:${target_name}>)
        add_custom_command(TARGET ${target_name} POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${target_name}> $<TARGET_FILE_DIR:${target_name}> 
          COMMAND_EXPAND_LISTS
        )
        install(TARGETS ${target_name}
           RUNTIME_DEPENDENCIES
           PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
           POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
           DIRECTORIES "$<TARGET_FILE_DIR:${target_name}>"
       )
      endif()
      install(TARGETS ${target_name} RUNTIME DESTINATION bin CONFIGURATIONS Release)
    endforeach()
    target_link_libraries(displayObj PRIVATE GLEW::GLEW)



    set(exe largeModalDeformationFactory) #requires wxwidget
    file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/${exe}/*.cpp")
    add_executable(${exe} WIN32 ${sources})
    target_link_libraries(${exe} PRIVATE ${PROJECT_NAME} 
        OpenGL::GL 
        GLUT::GLUT
        $<$<PLATFORM_ID:Linux,Darwin>:X11::X11>
    )
    find_package(wxWidgets CONFIG REQUIRED)
    target_link_libraries(largeModalDeformationFactory PRIVATE wx::core wx::base wx::gl)

    if($<TARGET_RUNTIME_DLLS:${exe}>)
        add_custom_command(TARGET ${exe} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${exe}> $<TARGET_FILE_DIR:${exe}> 
            COMMAND_EXPAND_LISTS
        )
        install(TARGETS ${exe}
            RUNTIME_DEPENDENCIES
            PRE_EXCLUDE_REGEXES "api-ms-" "ext-ms-"
            POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
            DIRECTORIES "$<TARGET_FILE_DIR:${exe}>"
        )
    endif()
    install(TARGETS ${exe} RUNTIME DESTINATION bin CONFIGURATIONS Release)
endif()


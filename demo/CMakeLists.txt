file(GLOB SRC_DEMO "*.cpp") # get all demo source files
foreach (FILE ${SRC_DEMO}) 
    get_filename_component(NAME_WITH_EXT ${FILE} NAME)
    string(REPLACE ".cpp" "" TARGET_NAME "${NAME_WITH_EXT}") # set demo application name as the filename without the extension
    add_executable (${TARGET_NAME} ${FILE})
    target_link_libraries(${TARGET_NAME} PUBLIC ctrModelLib) # link the ctrModelLib
    if(MSVC) # enable optimization
      target_compile_options(${TARGET_NAME} PRIVATE /O2 /DEIGEN_NO_DEBUG)
    else()
      target_compile_options(${TARGET_NAME} PRIVATE -O3 -DEIGEN_NO_DEBUG)
    endif()
endforeach ()
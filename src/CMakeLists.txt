file(GLOB SRC "*.cpp")

add_library(ctrModelLib STATIC ${SRC})

if(MSVC)
  target_compile_options(ctrModelLib PRIVATE /O2 /DEIGEN_NO_DEBUG)
else()
  target_compile_options(ctrModelLib PRIVATE -O3 -DEIGEN_NO_DEBUG)
endif()
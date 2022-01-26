SET(PROJECT_NAME Singular_Value_Decomposition)
SET(TEST_NAMES "UnitTest;SIMDTest;StreamTest;ThreadTest")

add_definitions(-DENABLE_AVX_INSTRUCTION_SET)
add_definitions(-DENABLE_MIC_INSTRUCTION_SET)

foreach(TEST_NAME ${TEST_NAMES})
  message("creating target for ${PROJECT_NAME}_${TEST_NAME}")
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.cpp)
    add_executable(${PROJECT_NAME}_${TEST_NAME}
      ${TEST_NAME}.cpp
      ../../Kernels/${PROJECT_NAME}/${PROJECT_NAME}.cpp
      ../../References/${PROJECT_NAME}/${PROJECT_NAME}_Reference.cpp

      ../../TestDeps/PTHREAD_QUEUE.cpp
      )

    target_include_directories(${PROJECT_NAME}_${TEST_NAME}
      PUBLIC ../../Kernels/${PROJECT_NAME}
      PUBLIC ../../References/${PROJECT_NAME}
      )
  else()
    message("${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.cpp does not exit")
  endif()
endforeach()

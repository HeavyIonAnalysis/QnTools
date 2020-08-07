set(SOURCES TestCuts.cpp)
add_executable(QnToolsCorrectionTests ${SOURCES})
target_link_libraries(QnToolsCorrectionTests  PRIVATE gtest_main gmock QnToolsCorrection QnToolsDataFrame PRIVATE ROOT::Core ROOT::Hist)
gtest_add_tests(TARGET QnToolsCorrectionTests  )

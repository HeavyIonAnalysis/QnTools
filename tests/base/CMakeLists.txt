
add_executable(TestQnToolsBase TestAxis.cpp)
target_link_libraries(TestQnToolsBase PRIVATE gtest_main QnToolsBase ROOT::Core)
gtest_add_tests(TARGET TestQnToolsBase)

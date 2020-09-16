project(QnToolsBase)
set(SOURCES QVector.cpp DataContainerHelper.cpp StatCollect.cpp StatCalculate.cpp Statistics.cpp BootStrap.cpp)
set(HEADERS DataContainer.hpp DataContainerHelper.hpp Axis.hpp QVector.hpp
             Cuts.hpp Statistics.hpp StatCalculate.hpp StatCollect.hpp BootStrap.hpp Stat.hpp)
set(HEADERSQNTOOLS ${HEADERS})
list(TRANSFORM HEADERSQNTOOLS PREPEND QnTools/)

include_directories(${QnTools_SOURCE_DIR}/include/base)
ROOT_GENERATE_DICTIONARY(G__${PROJECT_NAME} ${HEADERSQNTOOLS} LINKDEF ${PROJECT_NAME}LinkDef.h)

add_library(${PROJECT_NAME} SHARED ${SOURCES} G__${PROJECT_NAME})
add_library(QnTools::Base ALIAS ${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY EXPORT_NAME Base)

target_compile_definitions(${PROJECT_NAME} PUBLIC "-DUSE_ROOT")
target_include_directories(${PROJECT_NAME}
        PRIVATE
            ${QnTools_SOURCE_DIR}/include/base/QnTools
        PUBLIC
            $<BUILD_INTERFACE:${QnTools_SOURCE_DIR}/include/base/QnTools>
            $<BUILD_INTERFACE:${QnTools_SOURCE_DIR}/include/base>
            $<INSTALL_INTERFACE:include/QnTools>
        )
target_link_libraries(${PROJECT_NAME} PRIVATE ROOT::Core ROOT::Gpad)

# link to the gtest library if test are enabled
IF (QnTools_BUILD_TESTS)
    list(APPEND SOURCES QVector.test.cpp
                        DataContainer.test.cpp
                        Statistics.test.cpp
                        StatCollect.test.cpp
                        StatCalculate.test.cpp)
    add_executable(${PROJECT_NAME}UnitTests ${SOURCES})
    target_link_libraries(${PROJECT_NAME}UnitTests PRIVATE gtest_main ${PROJECT_NAME} ROOT::Core ROOT::Gpad ROOT::Hist)
    gtest_add_tests(TARGET ${PROJECT_NAME}UnitTests)
ENDIF(QnTools_BUILD_TESTS)

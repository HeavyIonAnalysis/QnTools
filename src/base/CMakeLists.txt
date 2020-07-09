project(QnToolsBase)
set(SOURCES QVector.cpp DataContainerHelper.cpp Stats.cpp ReSamples.cpp)
set(HEADERS DataContainer.hpp DataContainerHelper.hpp Axis.hpp QVector.hpp Stats.hpp Cuts.hpp Statistic.hpp ReSamples.hpp)

ROOT_GENERATE_DICTIONARY(G__${PROJECT_NAME} ${HEADERS} LINKDEF ${PROJECT_NAME}LinkDef.h)

add_library(${PROJECT_NAME} SHARED ${SOURCES} G__${PROJECT_NAME})
add_library(QnTools::Base ALIAS ${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY EXPORT_NAME Base)

target_compile_definitions(${PROJECT_NAME} PUBLIC "-DUSE_ROOT")
target_include_directories(${PROJECT_NAME}
        PRIVATE
        ${ROOT_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
        PUBLIC
        $<INSTALL_INTERFACE:include/QnTools>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )
target_link_libraries(${PROJECT_NAME} PUBLIC ${ROOT_LIBRARIES})

# link to the gtest library if test are enabled
IF (CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND SOURCES QVector.test.cpp
                        DataContainer.test.cpp)
    add_executable(${PROJECT_NAME}UnitTests ${SOURCES})
    target_link_libraries(${PROJECT_NAME}UnitTests PRIVATE gtest_main PUBLIC ${PROJECT_NAME})
    gtest_add_tests(TARGET ${PROJECT_NAME}UnitTests)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

project(QnToolsDataFrame)

add_library(${PROJECT_NAME} INTERFACE)

set(dataframe_headers common/AverageHelper.hpp
        common/AxesConfiguration.hpp
        common/EqualBinningHelper.hpp
        common/TemplateHelpers.hpp
        correction/RecenterAction.hpp
        correlation/CorrelationAction.hpp
        correlation/CorrelationAction.hpp
        correlation/ReSampleHelper.hpp)

list(APPEND DataFrameIncludePath
        "${ROOT_INCLUDE_DIRS}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/common"
        "${CMAKE_CURRENT_SOURCE_DIR}/correction"
        "${CMAKE_CURRENT_SOURCE_DIR}/correlation")
target_include_directories(${PROJECT_NAME}
        INTERFACE
        $<INSTALL_INTERFACE:include/QnTools>
        "$<BUILD_INTERFACE:${DataFrameIncludePath}>"
        )
target_link_libraries(${PROJECT_NAME} INTERFACE QnTools::Base ${ROOT_LIBRARIES})
add_library(QnTools::DataFrame ALIAS ${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY EXPORT_NAME DataFrame)

# link to the gtest library if test are enabled
IF (CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND SOURCES correction/RecenterAction.test.cpp)
    add_executable(${PROJECT_NAME}UnitTests ${SOURCES})
    target_link_libraries(${PROJECT_NAME}UnitTests PRIVATE gtest_main PUBLIC ${PROJECT_NAME})
    gtest_add_tests(TARGET ${PROJECT_NAME}UnitTests)
ENDIF (CMAKE_BUILD_TYPE MATCHES Debug)

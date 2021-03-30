project(QnToolsDataFrame)

add_library(${PROJECT_NAME} INTERFACE)

set(dataframe_headers common/AverageHelper.hpp
        common/AxesConfiguration.hpp
        common/BinningEqualizer.hpp
        common/TemplateFunctions.hpp
        correction/CorrectionBuilder.hpp
        correction/CorrectionAction.hpp
        correction/QVectorFunctor.hpp
        correlation/CorrelationAction.hpp
        correlation/CorrelationBuilder.hpp
        correlation/ReSampleFunctor.hpp
        toymc/ParticleGeneratorFunctor.hpp
        toymc/TrackingDetectorFunctor.hpp
        toymc/ChannelDetectorFunctor.hpp)

list(APPEND DataFrameIncludePath
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/common"
        "${CMAKE_CURRENT_SOURCE_DIR}/correction"
        "${CMAKE_CURRENT_SOURCE_DIR}/correlation"
        "${CMAKE_CURRENT_SOURCE_DIR}/toymc")
target_include_directories(${PROJECT_NAME}
        INTERFACE
        $<INSTALL_INTERFACE:include/QnTools>
        "$<BUILD_INTERFACE:${DataFrameIncludePath}>"
        )
target_link_libraries(${PROJECT_NAME} INTERFACE QnTools::Base ROOT::ROOTDataFrame)
add_library(QnTools::DataFrame ALIAS ${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY EXPORT_NAME DataFrame)

# link to the gtest library if test are enabled
IF (QnTools_BUILD_TESTS)
    list(APPEND SOURCES
            common/AverageHelper.test.cpp
            common/AxesConfiguration.test.cpp
            common/QVectorFunctor.test.cpp
            common/QVectorNormalizationFunctor.test.cpp
            correction/CorrectionAction.test.cpp
            correction/CorrectionBuilder.test.cpp
            correlation/CorrelationAction.test.cpp
            correlation/CorrelationBuilder.test.cpp
            toymc/ParticleGeneratorFunctor.test.cpp
            toymc/TrackingDetectorFunctor.test.cpp
            toymc/ChannelDetectorFunctor.test.cpp
            )
    add_executable(${PROJECT_NAME}UnitTests ${SOURCES})
    target_link_libraries(${PROJECT_NAME}UnitTests PRIVATE QnTools::Base gtest_main ${PROJECT_NAME})
    gtest_add_tests(TARGET ${PROJECT_NAME}UnitTests)
ENDIF (QnTools_BUILD_TESTS)

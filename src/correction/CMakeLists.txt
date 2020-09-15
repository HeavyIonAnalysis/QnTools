project(QnToolsCorrection)
set(SOURCES
        CorrectionHistogramBase.cpp
        CorrectionHistogramChannelizedSparse.cpp
        CorrectionHistogramSparse.cpp
        GainEqualization.cpp
        CorrectionProfile3DCorrelations.cpp
        CorrectionProfileChannelized.cpp
        CorrectionProfileChannelizedIngress.cpp
        CorrectionProfileComponents.cpp
        CorrectionProfileCorrelationComponents.cpp
        SubEvent.cpp
        SubEventChannels.cpp
        Alignment.cpp
        SubEventTracks.cpp
        Recentering.cpp
        TwistAndRescale.cpp
        CorrectionManager.cpp
        QAHistogram.cpp
        Detector.cpp)

set(HEADERS
        CorrectionManager.hpp
        CorrectionOnInputData.hpp
        CorrectionOnQnVector.hpp
        CorrectionAxis.hpp
        CorrectionAxisSet.hpp
        CorrectionBase.hpp
        CorrectionHistogramBase.hpp
        CorrectionsSet.hpp
        CorrectionHistogramChannelizedSparse.hpp
        CorrectionHistogramSparse.hpp
        CorrectionProfile3DCorrelations.hpp
        CorrectionProfileChannelized.hpp
        CorrectionProfileChannelizedIngress.hpp
        CorrectionDataVector.hpp
        CorrectionProfileComponents.hpp
        CorrectionProfileCorrelationComponents.hpp
        SubEvent.hpp
        SubEventChannels.hpp
        SubEventTracks.hpp
        GainEqualization.hpp
        Alignment.hpp
        Recentering.hpp
        TwistAndRescale.hpp
        RunList.hpp
        DetectorList.hpp
        InputVariable.hpp
        QAHistogram.hpp
        )
set(HEADERSQNTOOLS ${HEADERS})
list(TRANSFORM HEADERSQNTOOLS PREPEND QnTools/)

include_directories(${QnTools_SOURCE_DIR}/include/correction)
ROOT_GENERATE_DICTIONARY(G__${PROJECT_NAME} ${HEADERSQNTOOLS} LINKDEF ${PROJECT_NAME}LinkDef.h
        OPTIONS -I ${QnTools_SOURCE_DIR}/include/base/QnTools)

add_library(${PROJECT_NAME} SHARED ${SOURCES} G__${PROJECT_NAME})
add_library(QnTools::Correction ALIAS ${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY EXPORT_NAME Correction)
target_compile_definitions(${PROJECT_NAME} PUBLIC "-DUSE_ROOT")
target_include_directories(${PROJECT_NAME}
        PRIVATE
        ${QnTools_SOURCE_DIR}/include/correction/QnTools
        PUBLIC
        $<INSTALL_INTERFACE:include/QnTools>
        $<BUILD_INTERFACE:${QnTools_SOURCE_DIR}/include/correction/QnTools>
        )
target_link_libraries(${PROJECT_NAME} PUBLIC QnTools::Base PRIVATE ROOT::Core ROOT::Hist)
# link to the gtest library if test are enabled
IF (CMAKE_BUILD_TYPE MATCHES Debug)
    list(APPEND SOURCES )
    add_executable(${PROJECT_NAME}UnitTests ${SOURCES})
    target_link_libraries(${PROJECT_NAME}UnitTests PRIVATE gtest_main ${PROJECT_NAME} ROOT::Core ROOT::Hist)
    gtest_add_tests(TARGET ${PROJECT_NAME}UnitTests)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

//
// Created by eugene on 23/07/2020.
//

#include <vector>

#include <CorrectionManager.hpp>
#include <TRandom.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using namespace ::testing;
using namespace Qn;

struct ICutFHolder {
  virtual bool staticCutF1(const double &v) = 0;
  virtual bool staticCutF2(double v) = 0;
  virtual bool staticCutF3(int v) = 0;
  virtual bool dynamicCutF(const std::vector<double> &vargs) = 0;
};

struct CutMocksHolder : public ICutFHolder {
  MOCK_METHOD(bool, staticCutF1, (const double &v), (override));
  MOCK_METHOD(bool, staticCutF2, (double v), (override));
  MOCK_METHOD(bool, staticCutF3, (int v), (override));
  MOCK_METHOD(bool, dynamicCutF, (const std::vector<double> &vargs),
              (override));
};

class CutsTest : public ::testing::Test {
 public:
 protected:
  void SetUp() override {
    manager = new CorrectionManager;
    manager->SetFillValidationQA(false);
    manager->SetFillOutputTree(false);
    manager->SetFillCalibrationQA(false);

    manager->AddVariable("phi", 0, 1);
    manager->AddVariable("pt", 1, 1);
    manager->AddDetector("DetPhi", DetectorType::TRACK, "phi", "Ones", {},
                        {1, 2, 3});
  }
  void TearDown() override {
    delete manager;
  }

  void Fill() {
    for (int i = 0; i < 1000; ++i) {
      manager->Reset();
      if (manager->ProcessEvent()) {
        manager->GetVariableContainer()[0] =
            gRandom->Uniform(-TMath::Pi(), TMath::Pi());
        manager->GetVariableContainer()[1] = gRandom->Exp(1.);
        manager->FillTrackingDetectors();
      }

      manager->ProcessCorrections();
    }
  }

 protected:
  CorrectionManager* manager;

};

TEST_F(CutsTest, StaticCut) {
  auto mocks = new CutMocksHolder;
  EXPECT_CALL(*mocks, staticCutF1).Times(1000).WillRepeatedly(Return(true));
  EXPECT_CALL(*mocks, staticCutF2).Times(1000).WillRepeatedly(Return(true));
  EXPECT_CALL(*mocks, staticCutF3).Times(1000).WillRepeatedly(Return(true));

  manager->AddCutOnDetector(
      "DetPhi", {"pt"},
      [mocks](const double &v) { return mocks->staticCutF1(v); }, "static_cut_1");

  manager->AddCutOnDetector(
      "DetPhi", {"pt"},
      [mocks](double v) { return mocks->staticCutF2(v); }, "static_cut_2");

  manager->AddCutOnDetector(
      "DetPhi", {"pt"},
      [mocks](int v) { return mocks->staticCutF3(v); }, "static_cut_2");

  manager->InitializeOnNode();

  Fill();

  delete mocks;
}

//TEST_F(CutsTest, DynamicCut) {
//  auto mocks = new CutMocksHolder;
//  EXPECT_CALL(*mocks, dynamicCutF).Times(1000).WillRepeatedly(Return(true));
//
//  manager->AddCutOnDetector(
//      "DetPhi", {"pt"},
//      [mocks](const std::vector<double>& args) { return mocks->dynamicCutF(args); }, "static_cut");
//
//  manager->InitializeOnNode();
//
//  Fill();
//
//  delete mocks;
//}

}  // namespace

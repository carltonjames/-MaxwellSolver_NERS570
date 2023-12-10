#include <gtest/gtest.h>
#include "FDTD.hpp"

class FDTDTest : public ::testing::Test {
protected:
    FDTD* fdtd;

    void SetUp() override {
        fdtd = new FDTD(/* parameters */);
    }

    void TearDown() override {
        delete fdtd;
    }
};

TEST_F(FDTDTest, ConstructorTest) {
    // Test if FDTD is constructed correctly
    // For example, check if dimensions are set correctly
    EXPECT_EQ(fdtd->getXDim(), expectedXDim);
    EXPECT_EQ(fdtd->getYDim(), expectedYDim);
    EXPECT_EQ(fdtd->getZDim(), expectedZDim);
}

TEST_F(FDTDTest, SetSourcesTest) {
    
    auto source1 = std::make_unique<Source>();
    auto source2 = std::make_unique<Source>();
    std::vector<std::unique_ptr<Source>> sources;
    sources.push_back(std::move(source1));
    sources.push_back(std::move(source2));

    fdtd->setSources(std::move(sources));
    
    EXPECT_EQ(fdtd->getSourceCount(), 2);
}

TEST_F(FDTDTest, RunSimulationTest) {
    int result = fdtd->runSimulation();
    EXPECT_EQ(result, 0);
}

TEST_F(FDTDTest, FieldUpdateTest) {
    fdtd->runSimulationStep();
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include "Point.hpp"
#include <gtest/gtest.h>

class PointTest : public ::testing::Test {
protected:
    Point point;
    CurrentLine currentLine;
    void SetUp() override {
        point.setLocation({ 1.0f, 2.0f, 3.0f });
        point.setCharge(5.0f);
    }
};

TEST_F(PointTest, ChargeAtLocation) {
    EXPECT_FLOAT_EQ(point.getCharge({ 1.0f, 2.0f, 3.0f }), 5.0f);
}

TEST_F(PointTest, ChargeOutsideLocation) {
    EXPECT_FLOAT_EQ(point.getCharge({ 4.0f, 5.0f, 6.0f }), 0.0f);
}

TEST_F(PointTest, GetName) {
    EXPECT_EQ(point.getName(), "point");
}

TEST_F(CurrentLineTest, ChargeDistribution) {
    EXPECT_FLOAT_EQ(currentLine.getCharge({ 50, 50, 10 }), 1.0f);
    EXPECT_FLOAT_EQ(currentLine.getCharge({ 51, 50, 10 }), 0.0f);
}

TEST_F(CurrentLineTest, GetName) {
    EXPECT_EQ(currentLine.getName(), "CurrentLine");
}
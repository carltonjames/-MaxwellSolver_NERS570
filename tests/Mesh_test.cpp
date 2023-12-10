#include <gtest/gtest.h>
#include "Mesh.hpp"

class MeshTest : public ::testing::Test {
protected:
    Mesh* mesh;

    void SetUp() override {
        mesh = new Mesh(/* parameters */);
    }

    void TearDown() override {
        delete mesh;
    }
};

TEST_F(MeshTest, ConstructorTest) {
    // Test if Mesh is constructed correctly
    // For example, check if dimensions are set correctly
    EXPECT_EQ(mesh->getXDim(), expectedXDim);
    EXPECT_EQ(mesh->getYDim(), expectedYDim);
    EXPECT_EQ(mesh->getZDim(), expectedZDim);
}

TEST_F(MeshTest, FieldInitializationTest) {
    // Test if fields are initialized to zero
    auto field = mesh->getField({ 0, 0, 0 });
    EXPECT_FLOAT_EQ(field.E[0], 0.0);
    EXPECT_FLOAT_EQ(field.E[1], 0.0);
    EXPECT_FLOAT_EQ(field.E[2], 0.0);
    
    EXPECT_FLOAT_EQ(field.B[0], 0.0);
    EXPECT_FLOAT_EQ(field.B[1], 0.0);
    EXPECT_FLOAT_EQ(field.B[2], 0.0);

    EXPECT_FLOAT_EQ(field.D[0], 0.0);
    EXPECT_FLOAT_EQ(field.D[1], 0.0);
    EXPECT_FLOAT_EQ(field.D[2], 0.0);

    EXPECT_FLOAT_EQ(field.H[0], 0.0);
    EXPECT_FLOAT_EQ(field.H[1], 0.0);
    EXPECT_FLOAT_EQ(field.H[2], 0.0);
}

TEST_F(MeshTest, AddSourceTest) {
    // Add a source and test if it's added correctly
    auto source = std::make_unique<Source>();
    mesh->addSource(std::move(source));

    // Check if the source count is correct
    EXPECT_EQ(mesh->getSourceCount(), 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
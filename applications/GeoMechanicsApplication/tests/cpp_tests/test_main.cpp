#include <iostream>

#include <gtest/gtest.h>

#include "includes/kernel.h"
#include "geo_mechanics_application.h"

class GeoEnvironment : public ::testing::Environment
{
public:
    ~GeoEnvironment() override {}

    void SetUp() override
    {
        mpGeoApp = Kratos::make_shared<Kratos::KratosGeoMechanicsApplication>();
        mKernel.ImportApplication(mpGeoApp);
    }
    
    void TearDown() override {}

private:
    Kratos::Kernel mKernel;
    Kratos::KratosGeoMechanicsApplication::Pointer mpGeoApp;
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    testing::AddGlobalTestEnvironment(new GeoEnvironment);
    return RUN_ALL_TESTS();
}

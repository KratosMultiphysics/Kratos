#include <iostream>

#include "includes/kernel.h"
#include "geo_mechanics_application.h"
#include "testing/tester.h"

int main(int argc, char* argv[]) {
    std::cout << "Starting KratosGeoMechanicsFastSuite tests" << std::endl;

    int return_value = EXIT_SUCCESS;

    try {
        Kratos::Kernel kernel;
        Kratos::KratosGeoMechanicsApplication::Pointer geo_app = Kratos::make_shared<Kratos::KratosGeoMechanicsApplication>();
        kernel.ImportApplication(geo_app);

        auto& tester = Kratos::Testing::Tester::GetInstance();
        tester.SetVerbosity(Kratos::Testing::Tester::Verbosity::FAILED_TESTS_OUTPUTS);

        return_value = tester.RunTestSuite("KratosGeoMechanicsFastSuite");
        // return_value = tester.RunTestCases("*CalculateInclinedNormalFlux*");

    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout  << std::endl << "Finished KratosGeoMechanicsFastSuite tests" << std::endl;

    return return_value;
}

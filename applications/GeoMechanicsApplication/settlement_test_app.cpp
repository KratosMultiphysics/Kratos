#include <iostream>
#include <string>

#include "custom_workflows/dgeosettlement.h"

using namespace Kratos;


int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Too few arguments: provide at least the working directory and a project parameters file" << std::endl;
        return 1;
    }

    auto dummy_log_callback = [](const char*){};
    auto dummy_report_progress = [](double){};
    auto dummy_report_textual_progress = [](const char*){};
    auto never_cancel = [](){ return false; };

    const std::string working_directory{argv[1]};
    KratosGeoSettlement geo_settlement;
    for (auto i = 2; i < argc; ++i) {
        const std::string project_parameters_file_name{argv[i]};
        std::cout << "About to run stage using '" << project_parameters_file_name << "'" << std::endl;
        try {
            geo_settlement.RunStage(working_directory,
                                    project_parameters_file_name,
                                    dummy_log_callback,
                                    dummy_report_progress,
                                    dummy_report_textual_progress,
                                    never_cancel);
        }
        catch (const std::exception& e) {
            std::cerr << "Exception caught: " << e.what() << std::endl;
            return 2;
        }
        catch (...) {
            std::cerr << "Non-standard exception caught" << std::endl;
            return 3;
        }
    }

    return 0;
}

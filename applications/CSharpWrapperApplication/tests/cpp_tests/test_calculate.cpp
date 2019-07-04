//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "external_interface.h"
#include "utilities/os_utilities.h"

namespace Kratos
{
    namespace Testing
    {
        void CreateMDPAFile()
        {
            std::filebuf buffer;
            buffer.open(OSUtilities::GetCurrentWorkingDir() + "/file.mdpa",std::ios::out);
            std::ostream os(&buffer);
            os << "Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties  0\n    DENSITY 2700.000000\n    YOUNG_MODULUS 7000000.000000\n    POISSON_RATIO 0.300000\n    BODY_FORCE [3] (0.000000,0.000000,0.000000)\n    THICKNESS 1.000000\nEnd Properties\n\nBegin Nodes\n        1        0.0        0.0         0.0                               //node number, coord x, cord y, coord z\n        2        1.0        0.0         0.0                               //node number, coord x, cord y, coord z\n        3        1.0        1.0         0.0                               //node number, coord x, cord y, coord z\n        4        0.0        1.0         1.0                               //node number, coord x, cord y, coord z\nEnd Nodes\n\nBegin Elements SmallDisplacementElement3D4N\n    1 0 1 2 3 4\nEnd Elements\n\nBegin SubModelPart BasePart // Note that this would be a sub sub modelpart\n    Begin SubModelPartNodes\n        1\n        2\n    End SubModelPartNodes\n    Begin SubModelPart inner_part\n        Begin SubModelPartNodes\n            1\n        End SubModelPartNodes\n    End SubModelPart\nEnd SubModelPart";
            buffer.close();
        }

        /**
        * Checks the correct work of update node position
        */
        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperUpdateNodePosition, KratosCSharpWrapperApplicationFastSuite)
        {
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<Element>::Has("SmallDisplacementElement3D4N"))
                return void();

            CreateMDPAFile();
            const std::string file_name = OSUtilities::GetCurrentWorkingDir() + "/file.mdpa";
            CSharpKratosWrapper::CSharpInterface::init(file_name.c_str());
            float *x = CSharpKratosWrapper::CSharpInterface::getXCoordinates();
            float *y = CSharpKratosWrapper::CSharpInterface::getYCoordinates();
            float *z = CSharpKratosWrapper::CSharpInterface::getZCoordinates();
            int n = CSharpKratosWrapper::CSharpInterface::getNodesCount();

            for (int i = 0; i < n; i++) {
                std::cout << x[i] << " " << y[i] << " " << z[i] << std::endl;
            }

            std::cout << "----------------------------------------" << std::endl;
            CSharpKratosWrapper::CSharpInterface::updateNodePos(2, x[2]-2, y[2], z[2]);
            CSharpKratosWrapper::CSharpInterface::updateNodePos(10, x[10], y[10]-2.5, z[10]);
            CSharpKratosWrapper::CSharpInterface::calculate();
            x = CSharpKratosWrapper::CSharpInterface::getXCoordinates();
            y = CSharpKratosWrapper::CSharpInterface::getYCoordinates();
            z = CSharpKratosWrapper::CSharpInterface::getZCoordinates();

            for (int i = 0; i < n; i++) {
                std::cout << x[i] << " " << y[i] << " " << z[i] << std::endl;
            }

            remove((OSUtilities::GetCurrentWorkingDir() + "/file.mdpa").c_str());
        }

    } // namespace Testing
}  // namespace Kratos.

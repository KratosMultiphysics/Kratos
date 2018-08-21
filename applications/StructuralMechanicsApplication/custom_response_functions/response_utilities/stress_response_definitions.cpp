//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


// System includes

// External includes

// Project includes
#include "stress_response_definitions.h"

namespace Kratos
{
    namespace StressResponseDefinitions
    {
        TracedStressType ConvertStringToTracedStressType(const std::string& Str)
        {
            if(Str == "FX")
                return TracedStressType::FX;
            else if(Str == "FY")
                return TracedStressType::FY;
            else if(Str == "FZ")
                return TracedStressType::FZ;
            else if(Str == "MX")
                return TracedStressType::MX;
            else if(Str == "MY")
                return TracedStressType::MY;
            else if(Str == "MZ")
                return TracedStressType::MZ;
            else if(Str == "FXX")
                return TracedStressType::FXX;
            else if(Str == "FXY")
                return TracedStressType::FXY;
            else if(Str == "FXZ")
                return TracedStressType::FXZ;
            else if(Str == "FYX")
                return TracedStressType::FYX;
            else if(Str == "FYY")
                return TracedStressType::FYY;
            else if(Str == "FYZ")
                return TracedStressType::FYZ;
            else if(Str == "FZX")
                return TracedStressType::FZX;
            else if(Str == "FZY")
                return TracedStressType::FZY;
            else if(Str == "FZZ")
                return TracedStressType::FZZ;
            else if(Str == "MXX")
                return TracedStressType::MXX;
            else if(Str == "MXY")
                return TracedStressType::MXY;
            else if(Str == "MXZ")
                return TracedStressType::MXZ;
            else if(Str == "MYX")
                return TracedStressType::MYX;
            else if(Str == "MYY")
                return TracedStressType::MYY;
            else if(Str == "MYZ")
                return TracedStressType::MYZ;
            else if(Str == "MZX")
                return TracedStressType::MZX;
            else if(Str == "MZY")
                return TracedStressType::MZY;
            else if(Str == "MZZ")
                return TracedStressType::MZZ;
            else if(Str == "PK2")
                return TracedStressType::PK2;
            else
                KRATOS_ERROR << "Chosen stress type \n" <<Str<<"\" is not available!" << std::endl;
        }

        StressTreatment ConvertStringToStressTreatment(const std::string& Str)
        {
            if(Str == "mean")
                return StressTreatment::Mean;
            else if(Str == "node")
                return StressTreatment::Node;
            else if(Str == "GP")
                return StressTreatment::GaussPoint;
            else
                KRATOS_ERROR << "Chosen stress treatment \n" <<Str<<"\" is not available!" << std::endl;
        }

    }  // namespace StressResponseDefinitions.

}  // namespace Kratos.


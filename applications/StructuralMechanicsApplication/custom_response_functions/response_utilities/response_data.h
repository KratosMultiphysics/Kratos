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

#if !defined(KRATOS_RESPONSE_DATA_H_INCLUDED )
#define  KRATOS_RESPONSE_DATA_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"


namespace Kratos
{
    enum TracedStressType
    {
        FX,
        FY,
        FZ,
        MX,
        MY,
        MZ,
        FXX,
        FXY,
        FXZ,
        FYX,
        FYY,
        FYZ,
        FZX,
        FZY,
        FZZ,
        MXX,
        MXY,
        MXZ,
        MYX,
        MYY,
        MYZ,
        MZX,
        MZY,
        MZZ,
        StressTypeNotAvailible
    };

    enum StressTreatment
    {
        mean,
        node,
        GP,
        StressTreatmentNotAvailible
    };


class ResponseData
{
public:

      ///@name Type Definitions
      ///@{

      /// Pointer definition of ResponseData
      KRATOS_CLASS_POINTER_DEFINITION(ResponseData);

      ///@}
      ///@name Life Cycle
      ///@{

     /// Default constructor.
      ResponseData()
      {
      }

      /// Destructor.
      virtual ~ResponseData()
      {
      }

      ///@}
      ///@name Operators
      ///@{

      ///@}
      ///@name Operations
      ///@{


    TracedStressType ConvertStressType(const std::string& Str)
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
        else
            return TracedStressType::StressTypeNotAvailible;	
    }

    StressTreatment ConvertStressTreatment(const std::string& Str)
    {	
        if(Str == "mean") 
            return StressTreatment::mean;
        else if(Str == "node")
            return StressTreatment::node;
        else if(Str == "GP")
            return StressTreatment::GP;
        else
            return StressTreatment::StressTreatmentNotAvailible;
    }
    
private:
};// class ResponseData

    
}  // namespace Kratos.

#endif // KRATOS_RESPONSE_DATA_H_INCLUDED  defined



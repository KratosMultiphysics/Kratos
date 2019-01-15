//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#ifndef KRATOS_COMPUTE_GRADIENT_NUMERICAL_PROCESS_H
#define KRATOS_COMPUTE_GRADIENT_NUMERICAL_PROCESS_H


#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

#include <string>
#include <iostream>
#include <sstream>


namespace Kratos
{

class GetEquationId: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(GetEquationId);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for GetEquationId Process
//     GetEquationId(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for GetEquationId Process
    GetEquationId(ModelPart& rModelPart,
                Vector& rResult,
                int rElementId                             
                    ):
        Process(),
        mrModelPart(rModelPart),
        mrResult(rResult),
        mrElementId(rElementId)
    {
    }

    /// Destructor.
    ~GetEquationId() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;

        auto it = mrModelPart.pGetElement(mrElementId);
        KRATOS_WATCH(it->Id())
        std::vector<std::size_t> result;
        it -> EquationIdVector(result,mrModelPart.GetProcessInfo());
        mrResult.resize(result.size(), false);
        KRATOS_WATCH(result.size())
        KRATOS_WATCH(mrResult.size())       

        for (std::size_t i_dof  = 0;i_dof<result.size();i_dof++){
            mrResult(i_dof)=result[i_dof];
        }
        KRATOS_WATCH(result)
                 
        KRATOS_CATCH("");
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GetEquationId";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GetEquationId";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Vector& mrResult;
    int mrElementId;
    Flags mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GetEquationId& operator=(GetEquationId const& rOther);

    /// Copy constructor.
    GetEquationId(GetEquationId const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GetEquationId& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GetEquationId& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_ComputeLift_LEVEL_SET_PROCESS_H

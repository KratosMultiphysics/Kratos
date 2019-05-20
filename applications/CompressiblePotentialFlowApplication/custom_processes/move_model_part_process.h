//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nuñez
//

#ifndef KRATOS_MOVE_MODEL_PART_PROCESS_H
#define KRATOS_MOVE_MODEL_PART_PROCESS_H


#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

class MoveModelPartProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(MoveModelPartProcess);

    // Constructor for MoveModelPartProcess Process
    MoveModelPartProcess(ModelPart& rModelPart,
                     Paramaters& ThisParameters
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "origin"                        : [0.0,0.0,0.0],
            "rotation_angle"                : 0.0,
            "sizing_multiplier"             : 1.0

        })" );
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mOrigin = ThisParameters["origin"].GetVector();
        mRotationAngle = ThisParameters["rotation_angle"].GetDouble();
        mSizingMultiplier = ThisParameters["sizing_multiplier"].GetDouble();
    }

    /// Destructor.
    ~MoveModelPartProcess() override {}

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;

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
        return "MoveModelPartProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MoveModelPartProcess";
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
    Vector mOrigin;
    double mRotationAngle;
    double mSizingMultiplier;
    Parameters mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MoveModelPartProcess& operator=(MoveModelPartProcess const& rOther);

    /// Copy constructor.
    MoveModelPartProcess(MoveModelPartProcess const& rOther);


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
                                  MoveModelPartProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MoveModelPartProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_KUTTA_CONDITION_PROCESS_H

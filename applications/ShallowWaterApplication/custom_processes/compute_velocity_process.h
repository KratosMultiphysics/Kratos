//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_COMPUTE_VELOCITY_PROCESS_H_INCLUDED
#define KRATOS_COMPUTE_VELOCITY_PROCESS_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

namespace Kratos {
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** 
 * @ingroup ShallowWaterApplication
 * @class ComputeVelocityProcess
 * @brief This process computes the velocity from the conserved variables using the mass matrix to avoid local instabilities
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) ComputeVelocityProcess : public Process {
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ComputeVelocityProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeVelocityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ComputeVelocityProcess(ModelPart& rModelPart, double Epsilon);

    /// Destructor.
    virtual ~ComputeVelocityProcess() = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()() {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() override;

    virtual void Clear();

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
    virtual std::string Info() const override {
        return "ComputeVelocityProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "ComputeVelocityProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    double mEpsilon;
    Matrix mMassMatrix;
    unsigned int mNumNodes;
    bool mInitializeWasPerformed = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void Initialize();

    void ComputeMassMatrix();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator
    // ComputeVelocityProcess& operator=(ComputeVelocityProcess const& rOther);

    /// Copy constructor
    // ComputeVelocityProcess(ComputeVelocityProcess const& rOther);

    ///@}

}; // Class ComputeVelocityProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ComputeVelocityProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ComputeVelocityProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COMPUTE_VELOCITY_PROCESS_H_INCLUDED  defined

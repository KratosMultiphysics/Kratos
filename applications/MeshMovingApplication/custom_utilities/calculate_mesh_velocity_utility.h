//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Andreas Winterstein
//


#if !defined(KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED )
#define  KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Short class definition.
/** Detail class definition.
 */
class CalculateMeshVelocityUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateMeshVelocityUtility
    KRATOS_CLASS_POINTER_DEFINITION(CalculateMeshVelocityUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CalculateMeshVelocityUtility(ModelPart& rModelPart,
                                 Parameters Settings);

    /// Destructor.
    virtual ~CalculateMeshVelocityUtility();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CalculateMeshVelocities();


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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    CalculateMeshVelocityUtility& operator=(CalculateMeshVelocityUtility const& rOther);

    /// Copy constructor.
    CalculateMeshVelocityUtility(CalculateMeshVelocityUtility const& rOther);


    ///@}

}; // Class CalculateMeshVelocityUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CalculateMeshVelocityUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CalculateMeshVelocityUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_MESH_VELOCITY_UTILITY_H_INCLUDED  defined

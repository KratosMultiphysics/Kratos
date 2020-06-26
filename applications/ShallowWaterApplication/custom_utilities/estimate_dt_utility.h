//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_ESTIMATE_DT_SHALLOW_H_INCLUDED
#define KRATOS_ESTIMATE_DT_SHALLOW_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
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

/// Utility to estimate the time step in terms of the courant number.
/** The velocity can be the sum of the convective velocity and the wave speed
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) EstimateDtShallow
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EstimateDtShallow
    KRATOS_CLASS_POINTER_DEFINITION(EstimateDtShallow);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EstimateDtShallow(ModelPart& rThisModelPart, Parameters ThisParameters);

    /// Destructor.
    ~EstimateDtShallow(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double EstimateDt() const;

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;


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

    const ModelPart& mrModelPart;
    bool mEstimateDt;
    double mConstantDt;
    double mCourant;
    bool mConsiderFroude;
    double mMinDt;
    double mMaxDt;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    double EstimateTimeStep() const;

    double NodalCharacteristicTime(const Node<3>& rNode, double gravity) const;

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
    EstimateDtShallow& operator=(EstimateDtShallow const& rOther);

    /// Copy constructor.
    EstimateDtShallow(EstimateDtShallow const& rOther);


    ///@}

}; // Class EstimateDtShallow

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                EstimateDtShallow& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const EstimateDtShallow& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ESTIMATE_DT_SHALLOW_H_INCLUDED  defined

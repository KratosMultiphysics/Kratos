//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//

#if !defined(KRATOS_MC_YIELD_CRITERION_H_INCLUDED)
#define      KRATOS_MC_YIELD_CRITERION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions

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
class MCYieldCriterion
    : public YieldCriterion
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MisesHuberYieldCriterion
    KRATOS_CLASS_POINTER_DEFINITION( MCYieldCriterion );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MCYieldCriterion();

    /// Initialization constructor.
    MCYieldCriterion(HardeningLawPointer pHardeningLaw);

    /// Copy constructor.
    MCYieldCriterion(MCYieldCriterion const& rOther);

    /// Assignment operator.
    MCYieldCriterion& operator=(MCYieldCriterion const& rOther);


    /// Destructor.
    ~MCYieldCriterion() override;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double& CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rCohesion, const double& rFrictionAngle) override;

    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

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

    double GetSmoothingLodeAngle();

    double GetPI();

    double GetSmoothingHiperbolic();
    
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MisesHuberYieldCriterion& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MisesHuberYieldCriterion& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRESCA_YIELD_CRITERION_H_INCLUDED  defined 



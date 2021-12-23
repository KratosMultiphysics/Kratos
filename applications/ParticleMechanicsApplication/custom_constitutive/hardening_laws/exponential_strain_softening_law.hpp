//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_EXPONENTIAL_STRAIN_SOFTENING_LAW_H_INCLUDED )
#define  KRATOS_EXPONENTIAL_STRAIN_SOFTENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hardening_laws/particle_hardening_law.hpp"

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
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) ExponentialStrainSofteningLaw
        : public ParticleHardeningLaw
{
public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExponentialStrainSofteningLaw
    KRATOS_CLASS_POINTER_DEFINITION( ExponentialStrainSofteningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExponentialStrainSofteningLaw();


    /// Copy constructor.
    ExponentialStrainSofteningLaw(ExponentialStrainSofteningLaw const& rOther);

    /// Assignment operator.
    ExponentialStrainSofteningLaw& operator=(ExponentialStrainSofteningLaw const& rOther);

    /// Destructor.
    ~ExponentialStrainSofteningLaw();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /*
    * @brief This function return the softening rule for internal friction angle, cohesion, and internal dilatancy angle:
    * @param[in/out] rHardening Exponential softening rate of change
    * @param[in] rAlpha Plastic deviatoric strain
    * @param[in] rThisVariable Identifier variable: INTERNAL_FRICTION_ANGLE, COHESION, INTERNAL_DILATANCY_ANGLE
    * @return Softening rate of change parameter
    */
    double& CalculateHardening(double &rHardening, const double &rAlpha, const Variable<double>& rThisVariable) override;

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
    // std::string Info() const;

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const;


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

}; // Class ExponentialStrainSoftening

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   LinearIsotropicHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const LinearIsotropicHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXPONENTIAL_STRAIN_SOFTENING_LAW_H_INCLUDED defined


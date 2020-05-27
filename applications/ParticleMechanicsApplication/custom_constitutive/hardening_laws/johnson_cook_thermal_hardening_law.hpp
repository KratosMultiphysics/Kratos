//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

#if !defined(KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED

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
class JohnsonCookThermalHardeningLaw
	: public ParticleHardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of JohnsonCookThermalHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalHardeningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    JohnsonCookThermalHardeningLaw();


    /// Copy constructor.
    JohnsonCookThermalHardeningLaw(JohnsonCookThermalHardeningLaw const& rOther);

    /// Assignment operator.
    JohnsonCookThermalHardeningLaw& operator=(JohnsonCookThermalHardeningLaw const& rOther);

    /// Destructor.
    ~JohnsonCookThermalHardeningLaw() override;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this hardening law
     */
    ParticleHardeningLaw::Pointer Clone() const override;


    ///@}
    ///@name Operations
    ///@{

    //double& CalculateHardening(double &rHardening, const Parameters& rValues) override;

    double& CalculateHardening(double& rHardening,
        const double rPlasticStrainRate, const double rEquivalentPlasticStrain, const double rTemperature) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
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

}; // Class JohnsonCookThermalHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// // input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   JohnsonCookThermalHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const JohnsonCookThermalHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED  defined

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Daniel Diez
//

#if !defined (KRATOS_NEWTONIAN_TWO_FLUID_2D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_TWO_FLUID_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "newtonian_2d_law.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

/**
 * Defines a Newtonian constitutive law in 2D for two-fluid formulations.
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY (read from the nodes!!)
 * 2) C_SMAGORINSKY
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NewtonianTwoFluid2DLaw : public Newtonian2DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef std::size_t             SizeType;
    
    /**
     * Counted pointer of NewtonianTwoFluid2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewtonianTwoFluid2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewtonianTwoFluid2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewtonianTwoFluid2DLaw (const NewtonianTwoFluid2DLaw& rOther);

    /**
     * Destructor.
     */
    ~NewtonianTwoFluid2DLaw() override;

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override;

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
    
    double ComputeEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const;

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

    void EvaluateInPoint(
        double& rResult,
        const Variable<double>& rVariable,
        ConstitutiveLaw::Parameters& rParameters) const;

    double EquivalentStrainRate(ConstitutiveLaw::Parameters& rParameters) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class NewtonianTwoFluid2DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_TWO_FLUID_2D_H_INCLUDED  defined 

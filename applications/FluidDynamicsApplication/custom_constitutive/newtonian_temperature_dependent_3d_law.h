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
//

#if !defined (KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_3D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "newtonian_3d_law.h"

namespace Kratos
{

/**
 * Defines a Newtonian constitutive law in 3D.
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NewtonianTemperatureDependent3DLaw : public Newtonian3DLaw

{
public:
    /**
     * Type Definitions
     */

    typedef std::size_t             SizeType;

    /**
     * Counted pointer of NewtonianTemperatureDependent3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewtonianTemperatureDependent3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewtonianTemperatureDependent3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewtonianTemperatureDependent3DLaw (const NewtonianTemperatureDependent3DLaw& rOther);


    /**
     * Destructor.
     */
    ~NewtonianTemperatureDependent3DLaw() override;

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;

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

    /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
    double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override;

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class NewtonianTemperatureDependent3DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_3D_H_INCLUDED  defined

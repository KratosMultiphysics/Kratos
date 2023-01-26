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

#if !defined (KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_2D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "newtonian_2d_law.h"

namespace Kratos
{
/**
 * Defines a Newtonian constitutive law for 2D
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NewtonianTemperatureDependent2DLaw : public Newtonian2DLaw
{
public:
    /**
     * Type Definitions
     */

    typedef std::size_t             SizeType;

    /**
     * Counted pointer of Newtonian3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewtonianTemperatureDependent2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewtonianTemperatureDependent2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewtonianTemperatureDependent2DLaw (const NewtonianTemperatureDependent2DLaw& rOther);

    /**
     * Destructor.
     */
    ~NewtonianTemperatureDependent2DLaw() override;

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

    double GetEffectiveViscosity(ConstitutiveLaw::Parameters &rParameters) const override;

    ///@}
    ///@name Protected Operations
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    ///@}
}; // Class NewtonianTemperatureDependent2DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_TEMPERATURE_DEPENDENT_LAW_2D_H_INCLUDED  defined

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_FROZEN_TURBULENCE_NEWTONIAN_LAW_H_INCLUDED)
#define KRATOS_RANS_FROZEN_TURBULENCE_NEWTONIAN_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/constitutive_law.h"

// Application includes


namespace Kratos
{
/**
 * @brief This class is extending Newtonian2DLaw in FluidDynamicsApplication
 *
 * This class is used to extend Newtonian2DLaw in FluidDynamicsApplication
 * to include turbulent viscosity from Bossinesq Hypothesis.
 *
 */
template<class TUnfrozenTurbulenceConsititutiveLaw>
class KRATOS_API(RANS_APPLICATION) RansFrozenTurbulenceNewtonianLaw : public TUnfrozenTurbulenceConsititutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TUnfrozenTurbulenceConsititutiveLaw;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Counted pointer of Newtonian3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(RansFrozenTurbulenceNewtonianLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    RansFrozenTurbulenceNewtonianLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    RansFrozenTurbulenceNewtonianLaw(const RansFrozenTurbulenceNewtonianLaw& rOther);

    /**
     * Destructor.
     */
    ~RansFrozenTurbulenceNewtonianLaw() = default;

    /**
     * Operations needed by the base class:
     */

    ///@}
    ///@name Operations
    ///@{

    void CalculateDerivative(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rFunctionVariable,
        const Variable<double>& rDerivativeVariable,
        double& rOutput) override;


    ///@}
    ///@name Access
    ///@{

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    std::string Info() const override;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
}; // Class RansFrozenTurbulenceNewtonianLaw
} // namespace Kratos.
#endif // KRATOS_RANS_FROZEN_TURBULENCE_NEWTONIAN_LAW_H_INCLUDED  defined

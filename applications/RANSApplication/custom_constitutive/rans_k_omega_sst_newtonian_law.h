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

#pragma once

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
template<unsigned int TDim, class TPrimalBaseType>
class KRATOS_API(RANS_APPLICATION) RansKOmegaSSTNewtonianLaw : public TPrimalBaseType
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TPrimalBaseType;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Counted pointer of Newtonian3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(RansKOmegaSSTNewtonianLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    RansKOmegaSSTNewtonianLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    RansKOmegaSSTNewtonianLaw(const RansKOmegaSSTNewtonianLaw& rOther);

    /**
     * Destructor.
     */
    ~RansKOmegaSSTNewtonianLaw();

    /**
     * Operations needed by the base class:
     */

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to perform all the checks
     * needed on the input provided. Checks can be "expensive" as the function
     * is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;


    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameters,
        const Variable<double>& rThisVariable,
        double& rValue) override;


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

protected:
    ///@name Protected Operations
    ///@{

    /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
    double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override;

    ///@}
private:
    ///@name Private Operations
    ///@{

    double CalculateTurbulentViscosity(ConstitutiveLaw::Parameters& rParameters) const;

    double GetDynamicViscosity(ConstitutiveLaw::Parameters& rParameters) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
}; // Class RansKOmegaSSTNewtonianLaw
} // namespace Kratos.
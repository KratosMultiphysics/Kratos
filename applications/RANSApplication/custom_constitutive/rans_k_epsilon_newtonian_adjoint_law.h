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

#if !defined (KRATOS_RANS_K_EPSILONNEWTONIAN_ADJOINT_LAW_H_INCLUDED)
#define  KRATOS_RANS_K_EPSILONNEWTONIAN_ADJOINT_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

// Application includes


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// This class contains the common infrastructure for fluid constitutive laws.
template<class TBaseClassType>
class KRATOS_API(RANS_APPLICATION) RansKEpsilonNewtonianAdjointLaw : public TBaseClassType
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(RansKEpsilonNewtonianAdjointLaw);

    using BaseType = TBaseClassType;

    using GeometryType = typename BaseType::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    RansKEpsilonNewtonianAdjointLaw(ConstitutiveLaw& rConstitutiveLaw);

    /// Copy constructor.
    RansKEpsilonNewtonianAdjointLaw (const RansKEpsilonNewtonianAdjointLaw& rOther);

    /// Destructor
    virtual ~RansKEpsilonNewtonianAdjointLaw() override;

    ///@}
    ///@name Operations
    ///@{

    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) override;

    double CalculateEffectiveViscosityDerivative(
        ConstitutiveLaw::Parameters& rValuesDerivative,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType NodeIndex,
        const Variable<double>& rDerivativeVariable) override;

    ///@}
    ///@name Input and output
    ///@{

    /// @return A short string identifying this constitutive law instance.
    std::string Info() const override;

    /// Print basic information about this constitutive law instance.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print detailed information about this constitutive law instance and its managed data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

}; // Class RansKEpsilonNewtonianAdjointLaw

}  // namespace Kratos.

#endif // KRATOS_RANS_K_EPSILONNEWTONIAN_ADJOINT_LAW_H_INCLUDED  defined

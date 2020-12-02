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

#if !defined (KRATOS_FLUID_ADJOINT_CONSTITUTIVE_LAW)
#define  KRATOS_FLUID_ADJOINT_CONSTITUTIVE_LAW

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// This class contains the common infrastructure for fluid constitutive laws.
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FluidAdjointConstitutiveLaw);

    using GeometryType = ConstitutiveLaw::GeometryType;

    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FluidAdjointConstitutiveLaw(ConstitutiveLaw& rConstitutiveLaw);

    /// Copy constructor.
    FluidAdjointConstitutiveLaw (const FluidAdjointConstitutiveLaw& rOther);

    /// Destructor
    virtual ~FluidAdjointConstitutiveLaw();

    ///@}
    ///@name Operations
    ///@{

    /// Initialize a new instance of this type of law
    virtual FluidAdjointConstitutiveLaw::Pointer Clone() const;

    virtual void CalculateMaterialResponseCauchyDerivative(
        ConstitutiveLaw::Parameters& rValuesDerivative,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType NodeIndex,
        const Variable<double>& rDerivativeVariable,
        const double EffectiveViscosity,
        const double EffectiveViscosityDerivative);

    virtual int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo);

    virtual double CalculateEffectiveViscosityDerivative(
        ConstitutiveLaw::Parameters& rValuesDerivative,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType NodeIndex,
        const Variable<double>& rDerivativeVariable);

    ///@}
    ///@name Inquiry
    ///@{

    ConstitutiveLaw& GetPrimalConstitutiveLaw();

    ///@}
    ///@name Input and output
    ///@{

    /// @return A short string identifying this constitutive law instance.
    virtual std::string Info() const;

    /// Print basic information about this constitutive law instance.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print detailed information about this constitutive law instance and its managed data.
    virtual void PrintData(std::ostream& rOStream) const;

    ///@}

protected:

    ///@name Protected Operations
    ///@{

    void NewtonianConstitutiveMatrixDerivatives2D(
        const double EffectiveViscosityDerivative,
        Matrix& rCDerivative) const;

    void NewtonianConstitutiveMatrixDerivatives3D(
        const double EffectiveViscosityDerivative,
        Matrix& rCDerivative) const;

    ///@}

private:

    ///@name Member Variables
    ///@{

    ConstitutiveLaw& mrPrimalConstitutiveLaw;

    ///@}

}; // Class FluidAdjointConstitutiveLaw

}  // namespace Kratos.

#endif // KRATOS_FLUID_ADJOINT_CONSTITUTIVE_LAW  defined

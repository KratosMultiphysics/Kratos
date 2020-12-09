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

#if !defined(KRATOS_RANS_CONSTITUTIVE_LAW_EXTENSION_H_INCLUDED)
#define KRATOS_RANS_CONSTITUTIVE_LAW_EXTENSION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/serializer.h"

namespace Kratos
{
///@name Classes
///@{

/**
 * @brief This class is extending FluidConstitutiveLaws in FluidDynamicsApplication
 *
 * This class is used to extend FluidConstitutiveLaws in FluidDynamicsApplication
 * to include turbulent viscosity from Bossinesq Hypothesis.
 *
 * @tparam TFluidConstitutiveLawType        FluidConstitutiveLaw type
 */
template <class TFluidConstitutiveLawType>
class RansConstitutiveLawExtension : public TFluidConstitutiveLawType
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TFluidConstitutiveLawType;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(RansConstitutiveLawExtension);

    ///@}
    ///@name Life Cycle
    ///@{

    RansConstitutiveLawExtension();

    RansConstitutiveLawExtension(const RansConstitutiveLawExtension& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    ~RansConstitutiveLawExtension() override;

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
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class RansConstitutiveLawExtension

///@}

} // namespace Kratos.
#endif // KRATOS_RANS_CONSTITUTIVE_LAW_EXTENSION_H_INCLUDED  defined

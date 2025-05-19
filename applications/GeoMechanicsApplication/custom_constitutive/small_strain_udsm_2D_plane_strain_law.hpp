// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// System includes
#include "includes/define.h"

// Project includes
#include "custom_constitutive/small_strain_udsm_3D_law.hpp"

namespace Kratos
{
///@addtogroup ConstitutiveModelsApplication
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
class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM2DPlaneStrainLaw : public SmallStrainUDSM3DLaw
{
public:
    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    /// The size type definition
    using SizeType = std::size_t;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_2D;

    /// Pointer definition of SmallStrainUDSM2DPlaneStrainLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM2DPlaneStrainLaw);

    //@}
    //@name Life Cycle
    //@{

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    using SmallStrainUDSM3DLaw::SetValue;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override { return Dimension; }

    /**
     * @brief Voigt tensor size:
     */
    [[nodiscard]] SizeType GetStrainSize() const override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override { return "SmallStrainUDSM2DPlaneStrainLaw"; }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "SmallStrainUDSM2DPlaneStrainLaw Data";
    }

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
    void SetExternalStressVector(Vector& rStressVector) override;

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallStrainUDSM3DLaw

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos
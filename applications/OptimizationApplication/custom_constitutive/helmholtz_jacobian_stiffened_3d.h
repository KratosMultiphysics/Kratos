//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//


#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
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

/**
 * @class HelmholtzJacobianStiffened3D
 * @ingroup OptimizationApplication
 * @brief This class defines constitutive model for shape filtering of solid 3D cases.
 * This current implementation is basically elastic isotropic materical which computes
 * young modules based on element's jacobian and helmholtz filter radius.
 * @details This class derives from the base constitutive law
 * @author Reza Najian Asl
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HelmholtzJacobianStiffened3D
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The process info type definition
    typedef ProcessInfo      ProcessInfoType;

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw         BaseType;

    /// The size type definition
    typedef std::size_t             SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    // Adding the respective using to avoid overload conflicts
    using BaseType::Has;
    using BaseType::GetValue;

    /// Counted pointer of HelmholtzJacobianStiffened3D
    KRATOS_CLASS_POINTER_DEFINITION( HelmholtzJacobianStiffened3D );

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HelmholtzJacobianStiffened3D();

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HelmholtzJacobianStiffened3D (const HelmholtzJacobianStiffened3D& rOther);

    /**
     * @brief Destructor.
     */
    ~HelmholtzJacobianStiffened3D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Dimension of the law:
    */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return VoigtSize;
    };

    /**
     * @brief It calculates the value of a specified variable (Matrix case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

private:

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class HelmholtzJacobianStiffened3D
}  // namespace Kratos.

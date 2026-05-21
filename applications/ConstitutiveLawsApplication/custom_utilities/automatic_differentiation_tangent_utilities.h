// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"


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
 * @class AdvancedConstitutiveLawUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the tangent tensor of the constitutive law
 * @details The methods are static, so it can be called without constructing the class
 * @details The python sympy files that generate the c++ code can be found in applications/ConstitutiveLawsApplication/python_scripts/symbolic_generation/iso_damage_tangent_tensor
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @tparam TSofteningType exponential or linear softening type of damage
 * @author Alejandro Cornejo
 */
template <class TYieldSurfaceType, SizeType TSofteningType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) AutomaticDifferentiationTangentUtilities
{
  public:
    ///@name Type definitions
    ///@{
    typedef std::size_t SizeType;

    typedef TYieldSurfaceType YieldSurfaceType;

    /// We define the dimension
    static constexpr SizeType Dimension =  YieldSurfaceType::Dimension;

    /// We define the Voigt size
    static constexpr SizeType VoigtSize = YieldSurfaceType::VoigtSize;

    /// The matrix type definition
    typedef Matrix MatrixType;

    /// the vector type definition
    typedef Vector VectorType;

    /// The definition of the bounded vector type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, VoigtSize> BoundedMatrixVoigtType;

    /// Node type definition
    typedef Node NodeType;

    /// Geometry definitions
    typedef Geometry<NodeType> GeometryType;


    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the second invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI2 The second invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateTangentTensorIsotropicDamage(ConstitutiveLaw::Parameters rValues);

};
}
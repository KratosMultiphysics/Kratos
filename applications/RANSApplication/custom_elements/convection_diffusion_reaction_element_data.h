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

#if !defined(KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_SCALAR_CONVECTION_DIFFUSION_REACTION_ELEMENT_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Base class to hold Convection-Diffusion-Reaction transport equation data
 *
 * This class is used as the base class for Convection-Diffusion-Reaction
 * transport equation data. The derrived classes should implement the following methods
 *
 *      1. CalculateConstants      -> This method is used to calculate constants. Called only once before
 *                                    the gauss point iteration loop. (This is a public member method)
 *      2. CalculateGaussPointData -> This method is used to calculate gauss point quantities. Usually this method
 *                                    is used to update member variables in the protected region in this class.
 *                                    This method is called for each gauss point in gauss point iteration loop.
 *                                    (This is a public member method)
 *      3. Check                   -> This is a public static method. Used for checking only.
 *      4. GetScalarVariable       -> This is a public static method. Returns the variable being solved in this
 *                                    Convection-Diffusion-Reaction transport equation.
 *      5. GetName                 -> This is a public static method. Returns just the name.
 *
 * These data containers are used in the following element formulations.
 *      @see ConvectionDiffusionReactionCrossWindStabilizedElement
 *      @see ConvectionDiffusionReactionElement
 *      @see ConvectionDiffusionReactionResidualBasedFluxCorrectedElement
 *
 * Followings are few example use cases
 *      @see KEpsilonElementData::KElementData
 *      @see KEpsilonElementData::EpsilonElementData
 *
 * @tparam TDim
 */
template<unsigned int TDim>
class ConvectionDiffusionReactionElementData
{
public:
    ///@name Type Definitions
    ///@{

    using GeometryType = Geometry<Node>;

    using ArrayD = array_1d<double, TDim>;

    ///@}
    ///@name Life Cycle
    ///@{

    ConvectionDiffusionReactionElementData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo)
        : mrGeometry(rGeometry),
          mrProperties(rProperties),
          mrConstitutiveLaw(*(rProperties.GetValue(CONSTITUTIVE_LAW)))
    {
        mConstitutiveLawParameters =
            ConstitutiveLaw::Parameters(rGeometry, rProperties, rProcessInfo);
    }

    virtual ~ConvectionDiffusionReactionElementData() = default;

    ///@}
    ///@name Operations
    ///@{

    virtual void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
    }

    ///@}
    ///@name Access
    ///@{

    ConstitutiveLaw::Parameters& GetConstitutiveLawParameters() { return mConstitutiveLawParameters; }

    ConstitutiveLaw& GetConstitutiveLaw() { return mrConstitutiveLaw; }

    const GeometryType& GetGeometry() const { return mrGeometry; }

    const Properties& GetProperties() const { return mrProperties; }

    ArrayD GetEffectiveVelocity() const { return mEffectiveVelocity; }

    double GetEffectiveKinematicViscosity() const { return mEffectiveKinematicViscosity; }

    double GetReactionTerm() const { return mReactionTerm; }

    double GetSourceTerm() const { return mSourceTerm; }

    ///@}

protected:
    ///@name Protected Members
    ///@{

    /**
     * These variables are supposed to be modified by derrived class's CalculateGaussPointData
     * method. "CalculateGaussPointData" is not available in the base class in order to avoid
     * runtime virtual table search because these classes are used templated in the elements.
     */

    ArrayD mEffectiveVelocity;
    double mEffectiveKinematicViscosity;
    double mReactionTerm;
    double mSourceTerm;

    ///@}

private:
    ///@name Private Members
    ///@{

    const GeometryType& mrGeometry;
    const Properties& mrProperties;
    ConstitutiveLaw& mrConstitutiveLaw;
    ConstitutiveLaw::Parameters mConstitutiveLawParameters;

    ///@}
};

///@}
} // namespace Kratos

#endif
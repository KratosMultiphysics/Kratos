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

#if !defined(KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVES_H)
#define KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/condition.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <
    unsigned int TDim,
    unsigned int TNumNodes>
class VMSMonolithicKBasedWallConditionDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Condition::PropertiesType;

    static constexpr IndexType TBlockSize = TDim + 1;

    static constexpr IndexType TElementLocalSize = TBlockSize * TNumNodes;

    using VectorF = BoundedVector<double, TElementLocalSize>;

    using ArrayD = array_1d<double, 3>;

    using VectorN = BoundedVector<double, TNumNodes>;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    static void InitializeCondition(
        Condition& rCondition,
        const ProcessInfo& rProcessInfo);

    ///@}
    ///@name Class Declarations
    ///@{

    class Data;

    ///@}
    ///@name Classes
    ///@{

    class ResidualsContributions
    {
    public:
        ///@name Life Cycle
        ///@{

        ResidualsContributions(Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void AddGaussPointResidualsContributions(
            VectorF& rResidual,
            const double W,
            const Vector& rN);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;

        ///@}
    };

    template<class TDerivativesType>
    class VariableDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        static constexpr IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        static constexpr IndexType TDerivativesSize = TDerivativeDimension * TNumNodes;

        static constexpr double VelocityDerivativeFactor = TDerivativesType::VelocityDerivativeFactor;

        static constexpr double TurbulentKineticEnergyDerivativeFactor = TDerivativesType::TurbulentKineticEnergyDerivativeFactor;

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Calculates derivatives of residual w.r.t. condition nodes
         *
         * This method is used to calculate gauss point residual derivatives w.r.t.
         * condition nodes and its directions.
         *
         * @param rResidualDerivative           Output residual derivative
         * @param ConditionNodeIndex            Condition node index
         * @param ParentElementNodeIndex        Parent element node index corresponding to condition node given by ConditionNodeIndex
         * @param DirectionIndex                Direction index
         * @param W                             Gauss point weight
         * @param rN                            Shape function values
         * @param WDerivative                   Gauss point weight derivative
         * @param DetJDerivative                Jacobian determinant derivative
         */
        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
            const IndexType ConditionNodeIndex,
            const IndexType ParentElementNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative);

        /**
         * @brief Calculates derivatives of residual w.r.t. parent only nodes
         *
         * This method is used to calculate gauss point residual derivatives w.r.t.
         * parent element only nodes for the given condition. (common nodes' derivatives should not be calculated here)
         *
         * All the conditions have their parent element which the condition belongs to.
         * This can be identified through NEIGHBOUR_ELEMENTS of the condition
         *
         * @param rResidualDerivative           Output residual derivative
         * @param ParentElementNodeIndex        Parent element node index
         * @param DirectionIndex                Direction index
         * @param W                             Gauss point weight    (calculated from condition geometry)
         * @param rN                            Shape function values (calculated from condition geometry)
         */
        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
            const IndexType ParentElementNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        ResidualsContributions mResidualWeightDerivativeContributions;

        ///@}
    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            Condition& rCondition,
            Element& rParentElement,
            const ProcessInfo& rProcessInfo);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const int Step = 0);

        ///@}

    private:
        ///@name Private Members
        ///@{

        Condition& mrCondition;
        Element& mrParentElement;

        IndexType mGaussIntegrationPointIndex;

        /// Primal data

        double mWallVelocityMagnitude;
        double mDensity;
        double mUTau;
        double mYPlus;
        double mWallHeight;
        double mKappa;
        double mBeta;
        double mYPlusLimit;
        double mTkeBasedUTau;
        double mUBasedUTau;
        double mInvKappa;
        double mNormalMagnitude;
        double mKinematicViscosity;
        double mCmu;
        double mCmu25;
        double mTurbulentKineticEnergy;

        ArrayD mWallVelocity;
        ArrayD mUnitNormal;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType>
        friend class VariableDerivatives;
        friend class ResidualsContributions;

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVES_H

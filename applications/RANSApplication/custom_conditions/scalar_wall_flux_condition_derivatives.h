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

#if !defined(KRATOS_SCALAR_WALL_FLUX_CONDITION_DERIVATIVES_H)
#define KRATOS_SCALAR_WALL_FLUX_CONDITION_DERIVATIVES_H

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
#include "custom_conditions/scalar_wall_flux_condition_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <
    unsigned int TDim,
    unsigned int TNumNodes,
    class TConditionDataType>
class ScalarWallFluxConditionDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Condition::PropertiesType;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

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
            VectorN& rResidual,
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

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorN& rResidualDerivative,
            const IndexType ConditionNodeIndex,
            const IndexType ParentElementNodeIndex,
            const IndexType DirectionIndex,
            const double W,
            const Vector& rN,
            const double WDerivative,
            const double DetJDerivative);

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorN& rResidualDerivative,
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
            const Condition& rCondition,
            const Element& rParentElement,
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

        const Condition& mrCondition;
        const Element& mrParentElement;
        TConditionDataType mrConditionData;

        ScalarWallFluxConditionData::Parameters mParameters;

        IndexType mGaussPointIndex;

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

#endif // KRATOS_SCALAR_WALL_FLUX_CONDITION_DERIVATIVES_H

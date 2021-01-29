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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/constitutive_law.h"
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
    unsigned int TNumNodes,
    class TElementDataType>
class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Element::PropertiesType;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    using MatrixNN = BoundedMatrix<double, TNumNodes, TNumNodes>;

    static constexpr IndexType TEquationOffset = TElementDataType::TEquationOffset;

    ///@}
    ///@name Static Operations
    ///@{

    static int Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    ///@}
    ///@name Class Declarations
    ///@{

    class Data;

    ///@}
    ///@name Classes
    ///@{

    class ResidualContributions
    {
    public:
        ///@name Life Cycle
        ///@{

        ResidualContributions(Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Vector& rResidual,
            const ProcessInfo& rProcessInfo);

        void AddResidualContribution(
            Vector& rResidual,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        void Finalize(
            Vector& rResidual,
            const ProcessInfo& rProcessInfo);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;

        ///@}
    };

    template<class TDerivativesType>
    class FirstDerivative
    {
    public:
        ///@name Type Definitions
        ///@{

        static constexpr double TSelfWeight = TDerivativesType::TSelfWeight;

        static constexpr IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        static constexpr IndexType TDerivativesSize = TDerivativeDimension * TNumNodes;

        ///@}
        ///name@ Life Cycle
        ///@{

        FirstDerivative(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Vector& rResidualDerivative,
            const ProcessInfo& rProcessInfo);

        void CalculateResidualDerivative(
            Vector& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        void Finalize(
            Vector& rResidualDerivative,
            const ProcessInfo& rProcessInfo);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;
        ResidualContributions mResidualWeightDerivativeContributions;
        BoundedVector<MatrixNN, TDerivativesSize> mPrimalMatrixDerivatives;

        ///@}
    };

    class SecondDerivative
    {
    public:
        ///@name Life Cycle
        ///@{

        SecondDerivative(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Vector& rResidualDerivative,
            const ProcessInfo& rProcessInfo);

        void AddResidualDerivativeContributions(
            Vector& rResidualDerivative,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        void Finalize(
            Vector& rResidualDerivative,
            const ProcessInfo& rProcessInfo);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;

        ///@}
    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            const Element& rElement,
            const ProcessInfo& rProcessInfo);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(const ProcessInfo& rProcessInfo);

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        TElementDataType mrElementData;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType>
        friend class FirstDerivative;
        friend class SecondDerivative;
        friend class ResidualContributions;

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H

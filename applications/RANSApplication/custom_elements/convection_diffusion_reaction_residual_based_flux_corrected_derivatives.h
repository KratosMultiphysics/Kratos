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
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_constitutive/fluid_constitutive_law.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Element::PropertiesType;

    static constexpr IndexType TBlockSize = 1;

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

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
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo);

        void AddResidualContribution(
            BoundedVector<double, TElementLocalSize>& rResidual,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}

    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
    };

    template<class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
    class VariableDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        constexpr static IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        constexpr static double TSelfWeight = (TEquationOffset == TDerivativeOffset);

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo);

        void CalculateResidualDerivative(
            BoundedVector<double, TElementLocalSize>& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative);

        ///@}

    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;
        ResidualContributions mResidualWeightDerivativeContributions;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
    };

    class SecondDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        ///@}
        ///@name Life Cycle
        ///@{

        SecondDerivatives(
            const Element& rElement,
            FluidConstitutiveLaw& rFluidConstitutiveLaw);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo);

        void AddResidualDerivativeContributions(
            Matrix& rOutput,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        FluidConstitutiveLaw& mrFluidConstitutiveLaw;

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
            FluidConstitutiveLaw& rFluidConstitutiveLaw);

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
        FluidConstitutiveLaw& mrFluidConstitutiveLaw;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
        friend class VariableDerivatives;
        friend class ResidualContributions;

        ///@}
    };

    ///@}
    ///@name Static Operations
    ///@{

    ///@}
};

} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_DERIVATIVES_H

//  KratosFluidDynamicsApplication
//
//  License:		 BSD License
//					 license: FluidDynamicsApplication/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_RANS_EVM_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_RANS_EVM_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>

// Project includes

// Application includes
#include "custom_elements/vms_adjoint_element.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief An adjoint element for discrete shape sensitivity of VMS fluid element.
 *
 * @see VMS monolithic fluid element
 */
template <unsigned int TDim, class TRANSEvmVMSAdjointData, unsigned int TNumNodes = TDim + 1>
class KRATOS_API(RANS_APPLICATION) RANSEvmVMSAdjoint : public VMSAdjointElement<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(RANSEvmVMSAdjoint);

    constexpr static unsigned int TBlockSize = TDim + 1;

    constexpr static unsigned int TFluidLocalSize = TBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    typedef Element::NodeType NodeType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::VectorType VectorType;

    typedef Element::MatrixType MatrixType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef Element::EquationIdVectorType EquationIdVectorType;

    typedef BoundedMatrix<double, TNumNodes, TDim> ShapeFunctionDerivativesType;

    typedef VMSAdjointElement<TDim> BaseType;

    using BaseType::CalculateFirstDerivativesLHS;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit RANSEvmVMSAdjoint(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    RANSEvmVMSAdjoint(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    RANSEvmVMSAdjoint(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    ~RANSEvmVMSAdjoint() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element of this type.
     *
     * @return pointer to the newly created element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief Checks for proper element geometry, nodal variables and dofs.
     *
     * @return 0 after successful completion.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        BaseType::Check(rCurrentProcessInfo);

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            NodeType& r_node = this->GetGeometry()[iNode];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RELAXED_ACCELERATION, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        }
        return 0;

        KRATOS_CATCH("")
    }

    void CalculateFirstDerivativesLHS(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        this->CalculatePrimalGradientOfVMSSteadyTerm(rLeftHandSideMatrix, rCurrentProcessInfo);
        this->AddPrimalGradientOfVMSMassTerm(
            rLeftHandSideMatrix, RELAXED_ACCELERATION, -1.0, rCurrentProcessInfo);

        // Adding additional terms coming from the eddy viscosity hypothesis
        AddTurbulentViscosityPartialDerivativePrimalGradientOfVMSSteadyTerm(
            rLeftHandSideMatrix, rCurrentProcessInfo);
        AddTurbulentViscosityPartialDerivativePrimalGradientOfVMSMassTerm(
            rLeftHandSideMatrix, RELAXED_ACCELERATION, -1.0, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(rLeftHandSideMatrix); // transpose
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize> local_matrix;
            this->CalculateSensitivityMatrix(rSensitivityVariable, local_matrix,
                                             rCurrentProcessInfo);

            if (rOutput.size1() != local_matrix.size1() ||
                rOutput.size2() != local_matrix.size2())
                rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);
            noalias(rOutput) = local_matrix;
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            this->CalculateShapeGradientOfVMSSteadyTerm(rOutput, rCurrentProcessInfo);
            this->AddShapeGradientOfVMSMassTerm(rOutput, RELAXED_ACCELERATION,
                                                -1.0, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    void CalculateResidualScalarDerivatives(const Variable<double>& rVariable,
                                            BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        this->CalculateElementTotalSteadyResidualScalarDerivatives(
            rOutput, rVariable, rCurrentProcessInfo);
        this->AddElementTotalMassResidualScalarDerivatives(
            rOutput, RELAXED_ACCELERATION, rVariable, -1.0, rCurrentProcessInfo);
    }

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "RANSEvmVMSAdjoint" << this->GetGeometry().WorkingSpaceDimension()
               << "D #" << this->Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RANSEvmVMSAdjoint" << this->GetGeometry().WorkingSpaceDimension()
                 << "D #" << this->Id() << std::endl;
        rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber() << std::endl;
    }

    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
        rOStream << "Geometry Data: " << std::endl;
        this->GetGeometry().PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    virtual void CalculateElementData(TRANSEvmVMSAdjointData& rData,
                                      const Vector& rShapeFunctions,
                                      const Matrix& rShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "RANSEvmVMSAdjoint "
                        "CalculateElementData method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual void CalculateTurbulentKinematicViscosityVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const TRANSEvmVMSAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR
            << "Calling base "
               "CalculateTurbulentKinematicViscosityVelocityDerivatives "
               "method in "
               "RANSEvmVMSAdjoint "
               "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    virtual void CalculateTurbulentKinematicViscosityScalarDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TRANSEvmVMSAdjointData& rCurrentData,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base "
                        "CalculateTurbulentKinematicViscosityScalarDerivatives "
                        "method in "
                        "RANSEvmVMSAdjoint "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    void AddElementTotalMassResidualScalarDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutputMatrix,
        const Variable<array_1d<double, 3>>& rVariable,
        const Variable<double>& rDerivativeVariable,
        double alpha,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d<double, TNumNodes> N;
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

        TRANSEvmVMSAdjointData current_data;
        this->CalculateElementData(current_data, N, DN_DX, rCurrentProcessInfo);

        // Density
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity, VISCOSITY, N);
        Viscosity *= Density;

        // u
        array_1d<double, TDim> Velocity;
        this->EvaluateInPoint(Velocity, VELOCITY, N);

        // u * Grad(N)
        array_1d<double, TNumNodes> DensityVelGradN;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            DensityVelGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
                DensityVelGradN[i] += Density * DN_DX(i, d) * Velocity[d];
        }

        // Stabilization parameters
        double VelNorm = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
            VelNorm += Velocity[d] * Velocity[d];
        VelNorm = std::sqrt(VelNorm);
        const double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize, Density,
                                               Viscosity, rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be consistent
        // with the fluid element used to solve for VELOCITY and PRESSURE.
        BoundedVector<double, TNumNodes> TauOneDeriv;
        // BoundedVector<double, TNumNodes> TauTwoDeriv;
        TauOneDeriv.clear();

        BoundedVector<double, TNumNodes> NuTDerivative;
        this->CalculateTurbulentKinematicViscosityScalarDerivatives(
            NuTDerivative, rDerivativeVariable, current_data, rCurrentProcessInfo);

        const double CoefOne = -4.0 * Density * TauOne * TauOne / (ElemSize * ElemSize);

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            TauOneDeriv[i] = CoefOne * NuTDerivative[i];
        }

        // rVariable (x)
        array_1d<double, TDim> X;
        this->EvaluateInPoint(X, rVariable, N);

        // x * Grad(N)
        array_1d<double, TNumNodes> DensityXGradN;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            DensityXGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
                DensityXGradN[i] += Density * DN_DX(i, d) * X[d];
        }

        // Primal gradient of (lumped) VMS mass matrix multiplied with vector
        IndexType FirstRow(0);
        // Loop over nodes
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                for (IndexType m = 0; m < TDim; ++m)
                {
                    double valmn = 0.0;

                    valmn += DensityVelGradN[i] * TauOneDeriv[j] * Density * X[m];

                    // Adding it in a transposed manner
                    rOutputMatrix(j, FirstRow + m) += alpha * Volume * valmn;
                }

                rOutputMatrix(j, FirstRow + TDim) +=
                    alpha * Volume * DensityXGradN[i] * TauOneDeriv[j];
            } // Node block columns

            FirstRow += TBlockSize;
        } // Node block rows

        KRATOS_CATCH("")
    }

    void CalculateElementTotalSteadyResidualScalarDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rAdjointMatrix,
        const Variable<double>& rDerivativeVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rAdjointMatrix.clear();

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d<double, TNumNodes> N;
        double Volume;

        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

        TRANSEvmVMSAdjointData current_data;
        this->CalculateElementData(current_data, N, DN_DX, rCurrentProcessInfo);

        // Density
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity, VISCOSITY, N);
        Viscosity *= Density;

        // u
        array_1d<double, TDim> Velocity;
        this->EvaluateInPoint(Velocity, VELOCITY, N);

        // u * Grad(N)
        array_1d<double, TNumNodes> DensityVelGradN;
        noalias(DensityVelGradN) = Density * prod(DN_DX, Velocity);

        // Grad(u)
        BoundedMatrix<double, TDim, TDim> DensityGradVel;
        this->CalculateVelocityGradient(DensityGradVel, DN_DX);

        // Div(u)
        double DivVel = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
            DivVel += DensityGradVel(d, d);

        DensityGradVel *= Density;

        // Grad(p)
        array_1d<double, TDim> GradP;
        this->CalculatePressureGradient(GradP, DN_DX);

        // ( Grad(u) * Grad(N) )^T
        BoundedMatrix<double, TNumNodes, TDim> DN_DX_DensityGradVel;
        noalias(DN_DX_DensityGradVel) = prod(DN_DX, DensityGradVel);

        // ( u * Grad(u) * Grad(N) )^T
        array_1d<double, TNumNodes> DN_DX_DensityGradVel_Vel;
        noalias(DN_DX_DensityGradVel_Vel) = prod(DN_DX_DensityGradVel, Velocity);

        // u * Grad(u)
        array_1d<double, TDim> DensityGradVel_Vel;
        noalias(DensityGradVel_Vel) = prod(DensityGradVel, Velocity);

        // Grad(N)^T * Grad(p)
        array_1d<double, TNumNodes> DN_DX_GradP;
        noalias(DN_DX_GradP) = prod(DN_DX, GradP);

        // Grad(N)^T * BodyForce
        array_1d<double, TDim> BodyForce;
        array_1d<double, TNumNodes> DN_DX_BodyForce;
        this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
        BodyForce *= Density;
        noalias(DN_DX_BodyForce) = prod(DN_DX, BodyForce);

        // Stabilization parameters TauOne, TauTwo
        double VelNorm = norm_2(Velocity);
        double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize, Density,
                                               Viscosity, rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be
        // consistent with the fluid element used to solve for VELOCITY and
        // PRESSURE.
        BoundedVector<double, TNumNodes> TauOneDeriv;
        BoundedVector<double, TNumNodes> TauTwoDeriv;

        BoundedVector<double, TNumNodes> NuTDerivative;
        this->CalculateTurbulentKinematicViscosityScalarDerivatives(
            NuTDerivative, rDerivativeVariable, current_data, rCurrentProcessInfo);

        const double CoefOne = -4.0 * Density * TauOne * TauOne / (ElemSize * ElemSize);

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            TauOneDeriv[i] = CoefOne * NuTDerivative[i];
            TauTwoDeriv[i] = Density * NuTDerivative[i];
        }

        // Here, -(\partial R / \partial W) is calculated. This is the discrete
        // derivative of the fluid residual w.r.t the fluid variables and therefore
        // includes many of the terms defined in the fluid element. Neglecting the
        // transient terms of the fluid element, this matrix is identical to the
        // Jacobian of the fluid residual used for Newton-Raphson iterations. The
        // matrix is transposed at the end to get the adjoint system matrix.

        const double coeff = DivVel * Density * 2.0 / 3.0;

        IndexType FirstRow(0);
        // Loop over nodes
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                for (IndexType m = 0; m < TDim; ++m)
                {
                    const double dNi_dUm =
                        inner_prod(row(DN_DX, i), row(DensityGradVel, m));

                    double valmn = 0.0;

                    // Stabilization, lsq convection
                    // (u * Grad(v)) * TauOne * (u * Grad(u))
                    valmn += DensityVelGradN[i] * TauOneDeriv[j] * DensityGradVel_Vel[m];
                    // Stabilization, lsq divergence
                    // Div(v) * TauTwo * Div(u)
                    valmn += DN_DX(i, m) * TauTwoDeriv[j] * DivVel;

                    // Stabilization, convection-pressure
                    // (u * Grad(v)) * TauOne * Grad(p)
                    valmn += TauOneDeriv[j] * DensityVelGradN[i] * GradP[m];

                    // Stabilization, convection-BodyForce
                    // (u * Grad(v)) * TauOne * f
                    valmn -= DensityVelGradN[i] * TauOneDeriv[j] * BodyForce[m];

                    valmn += NuTDerivative[j] * DN_DX_DensityGradVel(i, m);
                    valmn += NuTDerivative[j] * dNi_dUm;
                    valmn -= NuTDerivative[j] * DN_DX(i, m) * coeff;

                    rAdjointMatrix(j, FirstRow + m) -= Volume * valmn;
                }

                double valpn = 0.0;

                // Stabilization, lsq pressure
                // TauOne * Grad(q) * Grad(p)
                valpn += DN_DX_GradP[i] * TauOneDeriv[j];

                // Stabilization, pressure-convection
                // Grad(q) * TauOne * (u * Grad(u))
                valpn += DN_DX_DensityGradVel_Vel[i] * TauOneDeriv[j];

                // Stabilization, pressure-BodyForce
                // Grad(q) * TauOne * f
                valpn -= DN_DX_BodyForce[i] * TauOneDeriv[j];

                rAdjointMatrix(j, FirstRow + TDim) -= Volume * valpn;
            } // Node block columns

            FirstRow += TBlockSize;
        } // Node block rows

        KRATOS_CATCH("")
    }

    /**
     * @brief Adds primal gradient of the VMS mass matrix multiplied by a vector.
     *
     * Calculates \f$ \partial_{\mathbf{w}^n} (\mathbf{M}^n\mathbf{x}) \f$.
     * \f$\mathbf{w}^n\f$ is the vector of primal variables at the current
     * adjoint step. \f$\mathbf{M}^n\f$ is the VMS mass matrix. \f$\mathbf{x}\f$
     * is a constant vector with zeros for pressure dofs. The variable
     * determines the values for velocity dofs.
     */

    void AddTurbulentViscosityPartialDerivativePrimalGradientOfVMSMassTerm(
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rOutputMatrix,
        const Variable<array_1d<double, 3>>& rVariable,
        double alpha,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d<double, TNumNodes> N;
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

        TRANSEvmVMSAdjointData current_data;
        this->CalculateElementData(current_data, N, DN_DX, rCurrentProcessInfo);

        // Density
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity, VISCOSITY, N);
        Viscosity *= Density;

        // u
        array_1d<double, TDim> Velocity;
        this->EvaluateInPoint(Velocity, VELOCITY, N);

        // u * Grad(N)
        array_1d<double, TNumNodes> DensityVelGradN;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            DensityVelGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
                DensityVelGradN[i] += Density * DN_DX(i, d) * Velocity[d];
        }

        // Stabilization parameters
        double VelNorm = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
            VelNorm += Velocity[d] * Velocity[d];
        VelNorm = std::sqrt(VelNorm);
        const double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize, Density,
                                               Viscosity, rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be consistent
        // with the fluid element used to solve for VELOCITY and PRESSURE.
        BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
        BoundedMatrix<double, TNumNodes, TDim> NuTDerivative;

        this->CalculateTurbulentKinematicViscosityVelocityDerivatives(
            NuTDerivative, current_data, rCurrentProcessInfo);

        const double CoefOne = -4.0 * Density * TauOne * TauOne / (ElemSize * ElemSize);

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                TauOneDeriv(i, d) = CoefOne * NuTDerivative(i, d);
            }
        }

        // rVariable (x)
        array_1d<double, TDim> X;
        this->EvaluateInPoint(X, rVariable, N);

        // x * Grad(N)
        array_1d<double, TNumNodes> DensityXGradN;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            DensityXGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; ++d)
                DensityXGradN[i] += Density * DN_DX(i, d) * X[d];
        }

        // Primal gradient of (lumped) VMS mass matrix multiplied with vector
        IndexType FirstRow(0), FirstCol(0);
        // Loop over nodes
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                for (IndexType m = 0; m < TDim; ++m)
                {
                    for (IndexType n = 0; n < TDim; ++n)
                    {
                        double valmn = 0.0;

                        valmn += DensityVelGradN[i] * TauOneDeriv(j, n) * Density * X[m];

                        // Adding it in a transposed manner
                        rOutputMatrix(FirstRow + m, FirstCol + n) += alpha * Volume * valmn;
                    }

                    rOutputMatrix(FirstRow + TDim, FirstCol + m) +=
                        alpha * Volume * DensityXGradN[i] * TauOneDeriv(j, m);
                }

                FirstCol += TBlockSize;
            } // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        } // Node block rows

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the elemental contribution to the steady adjoint system matrix.
     *
     * This function returns elemental contributions for:
     *
     * \f[
     * \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)
     * \f]
     *
     * where the current adjoint step is the \f$n^{th}\f$ time step.
     */
    void AddTurbulentViscosityPartialDerivativePrimalGradientOfVMSSteadyTerm(
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rAdjointMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d<double, TNumNodes> N;
        double Volume;

        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

        TRANSEvmVMSAdjointData current_data;
        this->CalculateElementData(current_data, N, DN_DX, rCurrentProcessInfo);

        // Density
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity, VISCOSITY, N);
        Viscosity *= Density;

        // u
        array_1d<double, TDim> Velocity;
        this->EvaluateInPoint(Velocity, VELOCITY, N);

        // u * Grad(N)
        array_1d<double, TNumNodes> DensityVelGradN;
        noalias(DensityVelGradN) = Density * prod(DN_DX, Velocity);

        // Grad(u)
        BoundedMatrix<double, TDim, TDim> DensityGradVel;
        this->CalculateVelocityGradient(DensityGradVel, DN_DX);

        // Div(u)
        double DivVel = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
            DivVel += DensityGradVel(d, d);

        DensityGradVel *= Density;

        // Grad(p)
        array_1d<double, TDim> GradP;
        this->CalculatePressureGradient(GradP, DN_DX);

        // ( Grad(u) * Grad(N) )^T
        BoundedMatrix<double, TNumNodes, TDim> DN_DX_DensityGradVel;
        noalias(DN_DX_DensityGradVel) = prod(DN_DX, DensityGradVel);

        // ( u * Grad(u) * Grad(N) )^T
        array_1d<double, TNumNodes> DN_DX_DensityGradVel_Vel;
        noalias(DN_DX_DensityGradVel_Vel) = prod(DN_DX_DensityGradVel, Velocity);

        // u * Grad(u)
        array_1d<double, TDim> DensityGradVel_Vel;
        noalias(DensityGradVel_Vel) = prod(DensityGradVel, Velocity);

        // Grad(N)^T * Grad(p)
        array_1d<double, TNumNodes> DN_DX_GradP;
        noalias(DN_DX_GradP) = prod(DN_DX, GradP);

        // Grad(N)^T * BodyForce
        array_1d<double, TDim> BodyForce;
        array_1d<double, TNumNodes> DN_DX_BodyForce;
        this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
        BodyForce *= Density;
        noalias(DN_DX_BodyForce) = prod(DN_DX, BodyForce);

        // Stabilization parameters TauOne, TauTwo
        double VelNorm = norm_2(Velocity);
        double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne, TauTwo, VelNorm, ElemSize, Density,
                                               Viscosity, rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be
        // consistent with the fluid element used to solve for VELOCITY and
        // PRESSURE.
        BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
        BoundedMatrix<double, TNumNodes, TDim> TauTwoDeriv;

        BoundedMatrix<double, TNumNodes, TDim> NuTDerivative;
        this->CalculateTurbulentKinematicViscosityVelocityDerivatives(
            NuTDerivative, current_data, rCurrentProcessInfo);

        const double CoefOne = -4.0 * Density * TauOne * TauOne / (ElemSize * ElemSize);

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                TauOneDeriv(i, d) = CoefOne * NuTDerivative(i, d);
                TauTwoDeriv(i, d) = Density * NuTDerivative(i, d);
            }
        }

        // Here, -(\partial R / \partial W) is calculated. This is the discrete
        // derivative of the fluid residual w.r.t the fluid variables and therefore
        // includes many of the terms defined in the fluid element. Neglecting the
        // transient terms of the fluid element, this matrix is identical to the
        // Jacobian of the fluid residual used for Newton-Raphson iterations. The
        // matrix is transposed at the end to get the adjoint system matrix.

        const double coeff = DivVel * Density * 2.0 / 3.0;

        IndexType FirstRow(0), FirstCol(0);
        // Loop over nodes
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                for (IndexType m = 0; m < TDim; ++m)
                {
                    const double dNi_dUm =
                        inner_prod(row(DN_DX, i), row(DensityGradVel, m));

                    for (IndexType n = 0; n < TDim; ++n)
                    {
                        double valmn = 0.0;

                        // Stabilization, lsq convection
                        // (u * Grad(v)) * TauOne * (u * Grad(u))
                        valmn += DensityVelGradN[i] * TauOneDeriv(j, n) *
                                 DensityGradVel_Vel[m];
                        // Stabilization, lsq divergence
                        // Div(v) * TauTwo * Div(u)
                        valmn += DN_DX(i, m) * TauTwoDeriv(j, n) * DivVel;

                        // Stabilization, convection-pressure
                        // (u * Grad(v)) * TauOne * Grad(p)
                        valmn += TauOneDeriv(j, n) * DensityVelGradN[i] * GradP[m];

                        // Stabilization, convection-BodyForce
                        // (u * Grad(v)) * TauOne * f
                        valmn -= DensityVelGradN[i] * TauOneDeriv(j, n) * BodyForce[m];

                        valmn += NuTDerivative(j, n) * DN_DX_DensityGradVel(i, m);
                        valmn += NuTDerivative(j, n) * dNi_dUm;
                        valmn -= NuTDerivative(j, n) * DN_DX(i, m) * coeff;

                        rAdjointMatrix(FirstRow + m, FirstCol + n) -= Volume * valmn;
                    }

                    double valpn = 0.0;

                    // Stabilization, lsq pressure
                    // TauOne * Grad(q) * Grad(p)
                    valpn += DN_DX_GradP[i] * TauOneDeriv(j, m);

                    // Stabilization, pressure-convection
                    // Grad(q) * TauOne * (u * Grad(u))
                    valpn += DN_DX_DensityGradVel_Vel[i] * TauOneDeriv(j, m);

                    // Stabilization, pressure-BodyForce
                    // Grad(q) * TauOne * f
                    valpn -= DN_DX_BodyForce[i] * TauOneDeriv(j, m);

                    rAdjointMatrix(FirstRow + TDim, FirstCol + m) -= Volume * valpn;
                }

                FirstCol += TBlockSize;
            } // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        } // Node block rows

        KRATOS_CATCH("")
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY;

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY;

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Unaccessible methods
    ///@{

    RANSEvmVMSAdjoint& operator=(RANSEvmVMSAdjoint const& rOther);

    RANSEvmVMSAdjoint(RANSEvmVMSAdjoint const& rOther);

    ///@}

}; // class RANSEvmVMSAdjoint

///@} // Kratos classes

///@name Input and output
///@{

/// Defines an input stream operator that does nothing.
template <unsigned int TDim, class TRANSEvmVMSAdjointData>
inline std::istream& operator>>(std::istream& rIStream,
                                RANSEvmVMSAdjoint<TDim, TRANSEvmVMSAdjointData>& rThis)
{
    return rIStream;
}

/// Defines an output stream operator that prints element info.
template <unsigned int TDim, class TRANSEvmVMSAdjointData>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RANSEvmVMSAdjoint<TDim, TRANSEvmVMSAdjointData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // FluidDynamicsApplication group

} // namespace Kratos

#endif // KRATOS_RANS_EVM_VMS_ADJOINT_ELEMENT_H_INCLUDED defined

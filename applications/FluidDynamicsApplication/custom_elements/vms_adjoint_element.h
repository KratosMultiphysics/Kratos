//  KratosFluidDynamicsApplication
//
//  License:		 BSD License
//					 license: FluidDynamicsApplication/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <array>

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/cfd_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/indirect_scalar.h"
#include "utilities/adjoint_extensions.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief An adjoint element for discrete shape sensitivity of VMS fluid element.
 *
 * @see VMS monolithic fluid element
 */
template< unsigned int TDim >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) VMSAdjointElement: public Element {

    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement) : mpElement{pElement}
        {
        }

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
            std::size_t index = 0;
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
            if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Z, Step);
            }
            rVector[index] = IndirectScalar<double>{}; // pressure
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
            std::size_t index = 0;
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
            if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Z, Step);
            }
            rVector[index] = IndirectScalar<double>{}; // pressure
        }

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension() + 1);
            std::size_t index = 0;
            rVector[index++] =
                MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
            rVector[index++] =
                MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
            if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Z, Step);
            }
            rVector[index] = IndirectScalar<double>{}; // pressure
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
        }

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
        }
    };

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VMSAdjointElement);

    constexpr static unsigned int TNumNodes = TDim + 1;

    constexpr static unsigned int TBlockSize = TDim + 1;

    constexpr static unsigned int TFluidLocalSize = TBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::VectorType VectorType;

    typedef std::array<double, TFluidLocalSize> ArrayType;

    typedef Element::MatrixType MatrixType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef std::array<Dof<double>::Pointer, TFluidLocalSize> DofsArrayType;

    typedef Element::EquationIdVectorType EquationIdVectorType;

    typedef std::array<std::size_t, TFluidLocalSize> EquationIdArrayType;

    typedef BoundedMatrix<double, TNumNodes, TDim>
    ShapeFunctionDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    VMSAdjointElement(IndexType NewId = 0) : Element(NewId)
    {
    }

    VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    VMSAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    ~VMSAdjointElement() override
    {}

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    }

    /**
     * @brief Creates a new element of this type.
     *
     * @return pointer to the newly created element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY

        return Element::Pointer(new VMSAdjointElement<TDim>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties));

        KRATOS_CATCH("")
    }

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY

        return Element::Pointer(new VMSAdjointElement<TDim>(NewId, pGeom, pProperties));

        KRATOS_CATCH("")
    }

    /**
     * @brief Checks for proper element geometry, nodal variables and dofs.
     *
     * @return 0 after successful completion.
     */
    int Check(const ProcessInfo &/*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        // Check the element id and geometry.
        ProcessInfo UnusedProcessInfo;
        int ReturnValue = Element::Check(UnusedProcessInfo);

        // Check if adjoint and fluid variables are defined.
        if (ADJOINT_FLUID_VECTOR_1.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_FLUID_VECTOR_1 Key is 0. "
                    "Check if the application was correctly registered.","");
        if (ADJOINT_FLUID_VECTOR_2.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_FLUID_VECTOR_2 Key is 0. "
                    "Check if the application was correctly registered.","");
        if (ADJOINT_FLUID_VECTOR_3.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_FLUID_VECTOR_3 Key is 0. "
                    "Check if the application was correctly registered.","");
        if (ADJOINT_FLUID_SCALAR_1.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_FLUID_SCALAR_1 Key is 0. "
                    "Check if the application was correctly registered.","");
        if (VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "VELOCITY Key is 0. "
                    "Check if the application was correctly registered.","");
        if (ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ACCELERATION Key is 0. "
                    "Check if the application was correctly registered.","");
        if (PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "PRESSURE Key is 0. "
                    "Check if the application was correctly registered.","");

        // Check if the nodes have adjoint and fluid variables and adjoint dofs.
        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_1) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_VECTOR_1 variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_2) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_VECTOR_2 variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_3) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_VECTOR_3 variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(ADJOINT_FLUID_SCALAR_1) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_SCALAR_1 variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing VELOCITY variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ACCELERATION variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing PRESSURE variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_X) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_Y) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_VECTOR_1_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_VECTOR_1 component degree of freedom on node ",
                        this->GetGeometry()[iNode].Id());
            if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_FLUID_SCALAR_1) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_FLUID_SCALAR_1 component degree of freedom on node ",
                        this->GetGeometry()[iNode].Id());
        }

        return ReturnValue;

        KRATOS_CATCH("")
    }

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        ArrayType values;
        this->GetValuesArray(values, Step);
        if (rValues.size() != TFluidLocalSize)
            rValues.resize(TFluidLocalSize, false);
        std::copy(values.begin(), values.end(), rValues.begin());
    }

    void GetValuesArray(ArrayType& rValues, int Step = 0)
    {
        GeometryType& rGeom = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            const array_1d<double, 3>& rVel =
                rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
            for (IndexType d = 0; d < TDim; d++)
                rValues[LocalIndex++] = rVel[d];
            rValues[LocalIndex++] =
                rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
        }
    }

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TFluidLocalSize)
            rValues.resize(TFluidLocalSize, false);

        rValues.clear();
    }

    void GetFirstDerivativesArray(ArrayType& rValues, int Step = 0)
    {
        std::fill(rValues.begin(), rValues.end(), 0.0);
    }

    /// Returns the adjoint acceleration values stored in this element's nodes.
    void GetSecondDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        ArrayType values;
        this->GetSecondDerivativesArray(values, Step);
        if (rValues.size() != TFluidLocalSize)
            rValues.resize(TFluidLocalSize, false);
        std::copy(values.begin(), values.end(), rValues.begin());
    }

    void GetSecondDerivativesArray(ArrayType& rValues, int Step = 0)
    {
        GeometryType& rGeom = this->GetGeometry();
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
        {
            const array_1d<double, 3>& rAccel =
                rGeom[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
            for (IndexType d = 0; d < TDim; d++)
                rValues[LocalIndex++] = rAccel[d];
            rValues[LocalIndex++] = 0.0; // pressure dof
        }
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

        KRATOS_CATCH("")
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        if (rLeftHandSideMatrix.size1() != TFluidLocalSize ||
            rLeftHandSideMatrix.size2() != TFluidLocalSize)
            rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);

        rLeftHandSideMatrix.clear();
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the adjoint matrix for velocity and pressure.
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * velocity and pressure transposed:
     *
     * \f[
     *    \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{w}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     *
     * where \f$\mathbf{w}^n\f$ is the vector of nodal velocities and pressures
     * stored at the current step. For steady problems, the ACCELERATION
     * (\f$\dot{\mathbf{w}}^n\f$) must be set to zero on the nodes. For
     * the Bossak method, \f$\dot{\mathbf{w}}^{n-\alpha}\f$ must be stored in
     * ACCELERATION.
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
        this->CalculateFirstDerivativesLHS(LHS, rCurrentProcessInfo);
        rLeftHandSideMatrix.resize(LHS.size1(), LHS.size2());
        noalias(rLeftHandSideMatrix) = LHS;
    }

    virtual void CalculateFirstDerivativesLHS(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLeftHandSideMatrix,
                                              ProcessInfo& rCurrentProcessInfo)
    {
        this->CalculatePrimalGradientOfVMSSteadyTerm(rLeftHandSideMatrix, rCurrentProcessInfo);
        this->AddPrimalGradientOfVMSMassTerm(rLeftHandSideMatrix, ACCELERATION,
                                             -1.0, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(rLeftHandSideMatrix); // transpose
    }

    /**
     * @brief Calculates the adjoint matrix for acceleration.
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * acceleration:
     *
     * \f[
     *    \partial_{\dot{\mathbf{w}}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\dot{\mathbf{w}}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
        this->CalculateSecondDerivativesLHS(LHS, rCurrentProcessInfo);
        rLeftHandSideMatrix.resize(LHS.size1(), LHS.size2(), false);
        noalias(rLeftHandSideMatrix) = LHS;
    }

    void CalculateSecondDerivativesLHS(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        this->CalculateVMSMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        rLeftHandSideMatrix = -trans(rLeftHandSideMatrix); // transpose
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

        KRATOS_CATCH("")
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::runtime_error,
                           "this function is not implemented.", "")

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the sensitivity matrix.
     *
     * \f[
     *    \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{s}}(\mathbf{M}^n \dot{\mathbf{w}}^{n-\alpha})^T
     * \f]
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize> local_matrix;
            this->CalculateShapeGradientOfVMSSteadyTerm(local_matrix, rCurrentProcessInfo);
            this->AddShapeGradientOfVMSMassTerm(local_matrix, ACCELERATION, -1.0, rCurrentProcessInfo);
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

    void GetDofList(DofsVectorType& rElementalDofList,
                    ProcessInfo& rCurrentProcessInfo) override
    {
        DofsArrayType dofs;
        this->GetDofArray(dofs, rCurrentProcessInfo);
        if (rElementalDofList.size() != dofs.size())
            rElementalDofList.resize(dofs.size());
        std::copy(dofs.begin(), dofs.end(), rElementalDofList.begin());
    }

    void GetDofArray(DofsArrayType& rElementalDofList,
                    ProcessInfo& /*rCurrentProcessInfo*/);

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override
    {
        EquationIdArrayType ids;
        this->EquationIdArray(ids, rCurrentProcessInfo);
        if (rResult.size() != ids.size())
            rResult.resize(ids.size());
        std::copy(ids.begin(), ids.end(), rResult.begin());
    }

    void EquationIdArray(EquationIdArrayType& rResult, ProcessInfo& /*rCurrentProcessInfo*/);

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "VMSAdjointElement" << this->GetGeometry().WorkingSpaceDimension()
        << "D #" << this->Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "VMSAdjointElement"
        << this->GetGeometry().WorkingSpaceDimension() << "D #"
        << this->Id() << std::endl;
        rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber()
        << std::endl;
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

    /// Calculate VMS-stabilized (lumped) mass matrix.
    void CalculateVMSMassMatrix(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rMassMatrix,
                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rMassMatrix.clear();

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

        // Density
        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
        Viscosity *= Density;

        // u
        array_1d< double, TDim > Velocity;
        this->EvaluateInPoint(Velocity,VELOCITY,N);

        // u * Grad(N)
        array_1d< double, TNumNodes > DensityVelGradN;
        for (IndexType i = 0; i < TNumNodes; i++)
        {
            DensityVelGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; d++)
                DensityVelGradN[i] += Density * DN_DX(i,d) * Velocity[d];
        }

        // Stabilization parameters
        double VelNorm = 0.0;
        for (IndexType d = 0; d < TDim; d++)
            VelNorm += Velocity[d] * Velocity[d];
        VelNorm = std::sqrt(VelNorm);
        const double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                Density,Viscosity,rCurrentProcessInfo);

        // Lumped mass
        const double LumpedMass = Density * Volume / static_cast<double>(TNumNodes);
        IndexType DofIndex = 0;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType d = 0; d < TDim; ++d)
            {
                rMassMatrix(DofIndex, DofIndex) += LumpedMass;
                ++DofIndex;
            }
            ++DofIndex; // Skip pressure Dof
        }

        // Stabilization, convection-acceleration
        IndexType FirstRow(0), FirstCol(0);
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                const double diag = DensityVelGradN[i] * TauOne * Density * N[j];

                for (IndexType d = 0; d < TDim; ++d)
                {
                    rMassMatrix(FirstRow+d,FirstCol+d) += Volume * diag;
                    rMassMatrix(FirstRow+TDim,FirstCol+d) +=
                            Volume * DN_DX(i,d) * TauOne * Density * N[j];
                }

                FirstCol += TBlockSize;
            }  // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        }  // Node block rows

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
    void AddPrimalGradientOfVMSMassTerm(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rOutputMatrix,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        double alpha,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

        // Density
        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
        Viscosity *= Density;

        // u
        array_1d< double, TDim > Velocity;
        this->EvaluateInPoint(Velocity,VELOCITY,N);

        // u * Grad(N)
        array_1d< double, TNumNodes > DensityVelGradN;
        for (IndexType i = 0; i < TNumNodes; i++)
        {
            DensityVelGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; d++)
                DensityVelGradN[i] += Density * DN_DX(i,d) * Velocity[d];
        }

        // Stabilization parameters
        double VelNorm = 0.0;
        for (IndexType d = 0; d < TDim; d++)
            VelNorm += Velocity[d] * Velocity[d];
        VelNorm = std::sqrt(VelNorm);
        const double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                Density,Viscosity,rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be consistent
        // with the fluid element used to solve for VELOCITY and PRESSURE.
        BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
        BoundedMatrix<double, TNumNodes, TDim> TauTwoDeriv;

        if (VelNorm > 0.0)
        {
            const double CoefOne =-2.0 * Density * TauOne * TauOne / (ElemSize * VelNorm);
            const double CoefTwo = 0.5 * Density * ElemSize / VelNorm;

            for (IndexType i = 0; i < TNumNodes; ++i)
            {
                for (IndexType d = 0; d < TDim; ++d)
                {
                    TauOneDeriv(i,d) = CoefOne * N[i] * Velocity[d];
                    TauTwoDeriv(i,d) = CoefTwo * N[i] * Velocity[d];
                }
            }
        }

        // rVariable (x)
        array_1d< double, TDim > X;
        this->EvaluateInPoint(X,rVariable,N);

        // x * Grad(N)
        array_1d< double, TNumNodes > DensityXGradN;
        for (IndexType i = 0; i < TNumNodes; i++)
        {
            DensityXGradN[i] = 0.0;
            for (IndexType d = 0; d < TDim; d++)
                DensityXGradN[i] += Density * DN_DX(i,d) * X[d];
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

                        valmn += DensityVelGradN[i] * TauOneDeriv(j,n) * Density * X[m];

                        valmn += Density * N[j] * DN_DX(i,n) * TauOne * Density * X[m];

                        rOutputMatrix(FirstRow+m,FirstCol+n) += alpha * Volume * valmn;
                    }

                    rOutputMatrix(FirstRow+TDim,FirstCol+m) +=
                    alpha * Volume * DensityXGradN[i] * TauOneDeriv(j,m);
                }

                FirstCol += TBlockSize;
            }  // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        }  // Node block rows

        KRATOS_CATCH("")
    }

    /**
     * @brief Adds shape gradient of the VMS mass matrix multiplied by a vector.
     *
     * Calculates \f$ \partial_{\mathbf{s}} (\mathbf{M}^n\mathbf{x})^T \f$.
     * \f$\mathbf{s}\f$ is the vector of nodal coordinates. \f$\mathbf{M}^n\f$.
     * is the VMS mass matrix at the current adjoint step. \f$\mathbf{x}\f$ is
     * a constant vector with zeros for pressure dofs. The variable
     * determines the values for velocity dofs.
     */
    void AddShapeGradientOfVMSMassTerm(BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutputMatrix,
                                       const Variable<array_1d<double, 3>>& rVariable,
                                       double alpha,
                                       const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

        // Density
        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
        Viscosity *= Density;

        // u
        array_1d< double, TDim > Velocity;
        this->EvaluateInPoint(Velocity,VELOCITY,N);

        // u * Grad(N)
        array_1d< double, TNumNodes > DensityVelGradN;
        noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

        // Det(J)
        const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
        array_1d< double, TCoordLocalSize > DetJDerivatives;
        this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

        // Stabilization parameters TauOne, TauTwo
        double VelNorm = norm_2(Velocity);
        double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                Density,Viscosity,rCurrentProcessInfo);

        // Vector values
        array_1d< double, TFluidLocalSize > X;
        IndexType DofIndex = 0;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            array_1d< double, 3 >& rValue = this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
            for (IndexType d = 0; d < TDim; ++d)
            {
                X[DofIndex++] = rValue[d];
            }
            X[DofIndex++] = 0.0; // pressure dof
        }

        array_1d< double, TFluidLocalSize > Derivative;

        // We compute the derivative w.r.t each coordinate of each node and
        // assign it to the corresponding row of the shape derivatives matrix.
        for (IndexType iCoord = 0; iCoord < TCoordLocalSize; ++iCoord)
        {
            // Det(J)'
            double DetJDeriv = DetJDerivatives[iCoord];

            // DN_DX'
            BoundedMatrix< double, TNumNodes, TDim > DN_DX_Deriv;
            for (IndexType i = 0; i < TNumNodes; ++i)
            {
                for (IndexType d = 0; d < TDim; ++d)
                {
                    DN_DX_Deriv(i,d) = -DN_DX(iCoord / TDim,d) * DN_DX(i,iCoord % TDim);
                }
            }

            // Volume'
            double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

            // u * Grad(N)'
            array_1d< double, TNumNodes > DensityVelGradNDeriv;
            noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv,Velocity);

            // TauOne', TauTwo'
            double TauOneDeriv, TauTwoDeriv;
            this->CalculateStabilizationParametersDerivative(
                    TauOneDeriv,TauTwoDeriv,TauOne,TauTwo,VelNorm,ElemSize,Density,
                    Viscosity,DetJDeriv);

            BoundedMatrix< double, TFluidLocalSize, TFluidLocalSize > LHS;
            array_1d< double, TFluidLocalSize > RHS;
            for (IndexType i=0; i < TFluidLocalSize; i++)
            {
                RHS[i] = 0.0;
                for (IndexType j=0; j < TFluidLocalSize; j++)
                    LHS(i,j) = 0.0;
            }

            // The usual lumped mass matrix
            double LumpedMassDeriv = Density * VolumeDeriv / static_cast<double>(TNumNodes);
            IndexType DofIndex = 0;
            for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
            {
                for (IndexType d = 0; d < TDim; ++d)
                {
                    LHS(DofIndex, DofIndex) += LumpedMassDeriv;
                    ++DofIndex;
                }
                ++DofIndex; // Skip pressure Dof
            }

            // Stabilization, convection-acceleration
            IndexType FirstRow(0), FirstCol(0);
            for (IndexType i = 0; i < TNumNodes; ++i)
            {
                for (IndexType j = 0; j < TNumNodes; ++j)
                {
                    double diag = 0.0;
                    double ddiag = 0.0;

                    diag += DensityVelGradN[i] * TauOne * Density * N[j];
                    ddiag += DensityVelGradNDeriv[i] * TauOne * Density * N[j];
                    ddiag += DensityVelGradN[i] * TauOneDeriv * Density * N[j];

                    for (IndexType n = 0; n < TDim; ++n)
                    {
                        double valn = DN_DX(i,n) * TauOne * Density * N[j];
                        double dvaln = 0.0;
                        dvaln += DN_DX_Deriv(i,n) * TauOne * Density * N[j];
                        dvaln += DN_DX(i,n) * TauOneDeriv * Density * N[j];

                        LHS(FirstRow+n,FirstCol+n) += VolumeDeriv * diag
                        + Volume * ddiag;

                        LHS(FirstRow+TDim,FirstCol+n) += VolumeDeriv * valn
                        + Volume * dvaln;
                    }

                    FirstCol += TBlockSize;
                } // Node block columns

                FirstRow += TBlockSize;
                FirstCol = 0;
            } // Node block rows

            // Assign the derivative w.r.t this coordinate to the
            // shape derivative mass matrix.
            noalias(Derivative) = prod(LHS,X);
            for (IndexType k = 0; k < TFluidLocalSize; ++k)
            {
                rOutputMatrix(iCoord,k) += alpha * Derivative[k];
            }
        }
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
    void CalculatePrimalGradientOfVMSSteadyTerm(
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rAdjointMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rAdjointMatrix.clear();

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

        // Density
        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
        Viscosity *= Density;

        // u
        array_1d< double, TDim > Velocity;
        this->EvaluateInPoint(Velocity,VELOCITY,N);

        // u * Grad(N)
        array_1d< double, TNumNodes > DensityVelGradN;
        noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

        // Grad(u)
        BoundedMatrix< double, TDim, TDim > DensityGradVel;
        this->CalculateVelocityGradient(DensityGradVel,DN_DX);

        // Div(u)
        double DivVel = 0.0;
        for (IndexType d = 0; d < TDim; ++d)
            DivVel += DensityGradVel(d,d);

        DensityGradVel *= Density;

        // Grad(p)
        array_1d< double, TDim > GradP;
        this->CalculatePressureGradient(GradP,DN_DX);

        // ( Grad(u) * Grad(N) )^T
        BoundedMatrix< double, TNumNodes, TDim > DN_DX_DensityGradVel;
        noalias(DN_DX_DensityGradVel) = prod(DN_DX,DensityGradVel);

        // ( u * Grad(u) * Grad(N) )^T
        array_1d< double, TNumNodes > DN_DX_DensityGradVel_Vel;
        noalias(DN_DX_DensityGradVel_Vel) = prod(DN_DX_DensityGradVel,Velocity);

        // u * Grad(u)
        array_1d< double, TDim > DensityGradVel_Vel;
        noalias(DensityGradVel_Vel) = prod(DensityGradVel,Velocity);

        // Grad(N)^T * Grad(p)
        array_1d< double, TNumNodes > DN_DX_GradP;
        noalias(DN_DX_GradP) = prod(DN_DX,GradP);

        // Grad(N)^T * BodyForce
        array_1d< double, TDim > BodyForce;
        array_1d< double, TNumNodes > DN_DX_BodyForce;
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
        BodyForce *= Density;
        noalias(DN_DX_BodyForce) = prod(DN_DX,BodyForce);

        // Stabilization parameters TauOne, TauTwo
        double VelNorm = norm_2(Velocity);
        double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                Density,Viscosity,rCurrentProcessInfo);

        // Derivatives of TauOne, TauTwo w.r.t velocity. These definitions
        // depend on the definitions of TauOne and TauTwo and should be consistent
        // with the fluid element used to solve for VELOCITY and
        // PRESSURE.
        BoundedMatrix<double, TNumNodes, TDim> TauOneDeriv;
        BoundedMatrix<double, TNumNodes, TDim> TauTwoDeriv;

        if (VelNorm > 0.0)
        {
            double CoefOne =-2.0 * Density * TauOne * TauOne / (ElemSize * VelNorm);
            double CoefTwo = 0.5 * Density * ElemSize / VelNorm;

            for (IndexType i = 0; i < TNumNodes; ++i)
            {
                for (IndexType d = 0; d < TDim; ++d)
                {
                    TauOneDeriv(i,d) = CoefOne * N[i] * Velocity[d];
                    TauTwoDeriv(i,d) = CoefTwo * N[i] * Velocity[d];
                }
            }
        }

        // Here, -(\partial R / \partial W) is calculated. This is the discrete
        // derivative of the fluid residual w.r.t the fluid variables and therefore
        // includes many of the terms defined in the fluid element. Neglecting the
        // transient terms of the fluid element, this matrix is identical to the
        // Jacobian of the fluid residual used for Newton-Raphson iterations. The
        // matrix is transposed at the end to get the adjoint system matrix.

        IndexType FirstRow(0), FirstCol(0);
        // Loop over nodes
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            for (IndexType j = 0; j < TNumNodes; ++j)
            {
                double diag = 0.0;

                // Convective term, v * (u * Grad(u))
                diag += N[i] * DensityVelGradN[j];

                // Stabilization, lsq convection
                // (u * Grad(v)) * TauOne * (u * Grad(u))
                diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];

                for (IndexType m = 0; m < TDim; ++m)
                {
                    for (IndexType n = 0; n < TDim; ++n)
                    {
                        double valmn = 0.0;

                        // Convective term, v * (u * Grad(u))
                        valmn += N[i] * N[j] * DensityGradVel(m,n);

                        // Stabilization, lsq convection
                        // (u * Grad(v)) * TauOne * (u * Grad(u))
                        valmn += DensityVelGradN[i] * TauOne * N[j] * DensityGradVel(m,n);
                        valmn += DensityVelGradN[i] * TauOneDeriv(j,n)
                        * DensityGradVel_Vel[m];
                        valmn += Density * N[j] * DN_DX(i,n) * TauOne
                        * DensityGradVel_Vel[m];

                        // Stabilization, lsq divergence
                        // Div(v) * TauTwo * Div(u)
                        valmn += DN_DX(i,m) * TauTwo * DN_DX(j,n);
                        valmn += DN_DX(i,m) * TauTwoDeriv(j,n) * DivVel;

                        // Stabilization, convection-pressure
                        // (u * Grad(v)) * TauOne * Grad(p)
                        valmn += TauOneDeriv(j,n) * DensityVelGradN[i] * GradP[m];
                        valmn += Density * TauOne * N[j] * DN_DX(i,n) * GradP[m];

                        // Stabilization, convection-BodyForce
                        // (u * Grad(v)) * TauOne * f
                        valmn -= N[j] * DN_DX(i,n) * TauOne * Density * BodyForce[m];
                        valmn -= DensityVelGradN[i] * TauOneDeriv(j,n) * BodyForce[m];

                        rAdjointMatrix(FirstRow+m,FirstCol+n) += Volume * valmn;
                    }

                    rAdjointMatrix(FirstRow+m,FirstCol+m) += Volume * diag;

                    double valmp = 0.0;
                    double valpn = 0.0;

                    // Pressure term
                    // Div(v) * p
                    valmp -= DN_DX(i,m) * N[j];

                    // Stabilization, convection-pressure
                    // (u * Grad(v)) * TauOne * Grad(p)
                    valmp += TauOne * DensityVelGradN[i] * DN_DX(j,m);

                    // Divergence term
                    // q * Div(u)
                    valpn += N[i] * DN_DX(j,m);

                    // Stabilization, lsq pressure
                    // TauOne * Grad(q) * Grad(p)
                    valpn += DN_DX_GradP[i] * TauOneDeriv(j,m);

                    // Stabilization, pressure-convection
                    // Grad(q) * TauOne * (u * Grad(u))
                    valpn += DN_DX(i,m) * TauOne * DensityVelGradN[j];
                    valpn += DN_DX_DensityGradVel(i,m) * TauOne * N[j];
                    valpn += DN_DX_DensityGradVel_Vel[i] * TauOneDeriv(j,m);

                    // Stabilization, pressure-BodyForce
                    // Grad(q) * TauOne * f
                    valpn -= DN_DX_BodyForce[i] * TauOneDeriv(j,m);

                    rAdjointMatrix(FirstRow+m,FirstCol+TDim) += Volume * valmp;
                    rAdjointMatrix(FirstRow+TDim,FirstCol+m) += Volume * valpn;
                }

                // Stabilization, lsq pressure
                // TauOne * Grad(q) * Grad(p)
                double valpp = 0.0;
                for (IndexType d = 0; d < TDim; ++d)
                {
                    valpp += DN_DX(i,d) * DN_DX(j,d);
                }
                valpp *= TauOne;

                rAdjointMatrix(FirstRow+TDim,FirstCol+TDim) += Volume * valpp;

                FirstCol += TBlockSize;
            }  // Node block columns

            FirstRow += TBlockSize;
            FirstCol = 0;
        }  // Node block rows

        // Viscous term
        this->AddViscousTerm(rAdjointMatrix,DN_DX,Viscosity * Volume);

        // change the sign for consistency with definition
        noalias(rAdjointMatrix) = -rAdjointMatrix;

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculate the partial derivatives of damping term w.r.t. shape parameters.
     *
     * This function returns elemental contributions for:
     *
     * \f[
     * \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     * \f]
     *
     * \f$\mathbf{s}\f$ are the coordinates of the element's nodes.
     *
     * This function is only valid when the determinant of the Jacobian is constant
     * over the element.
     */
    void CalculateShapeGradientOfVMSSteadyTerm(
        BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rShapeDerivativesMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get shape functions, shape function gradients and element volume (area in
        // 2D). Only one integration point is used so the volume is its weight.
        ShapeFunctionDerivativesType DN_DX;
        array_1d< double, TNumNodes > N;
        double Volume;

        GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

        // Density
        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        // Dynamic viscosity
        double Viscosity;
        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
        Viscosity *= Density;

        // u
        array_1d< double, TDim > Velocity;
        this->EvaluateInPoint(Velocity,VELOCITY,N);

        // u * Grad(N)
        array_1d< double, TNumNodes > DensityVelGradN;
        noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

        // Det(J)
        const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
        array_1d< double, TCoordLocalSize > DetJDerivatives;
        this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

        // Stabilization parameters TauOne, TauTwo
        double VelNorm = norm_2(Velocity);
        double ElemSize = this->CalculateElementSize(Volume);
        double TauOne, TauTwo;
        this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
                Density,Viscosity,rCurrentProcessInfo);

        // External body force
        array_1d< double, TDim > BodyForce;
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
        BodyForce *= Density;

        array_1d< double, TFluidLocalSize > FluidValues;

        IndexType DofIndex = 0;
        for (IndexType iNode = 0; iNode < TNumNodes; iNode++)
        {
            array_1d< double, 3 >& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            for (IndexType d = 0; d < TDim; d++)
                FluidValues[DofIndex++] = rVelocity[d];
            FluidValues[DofIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE);
        }

        // We compute the derivative of the residual w.r.t each coordinate of each
        // node and assign it to the corresponding row of the shape derivatives
        // matrix.
        for (IndexType iCoord = 0; iCoord < TCoordLocalSize; ++iCoord)
        {
            // Det(J)'
            double DetJDeriv = DetJDerivatives[iCoord];

            // DN_DX'
            BoundedMatrix< double, TNumNodes, TDim > DN_DX_Deriv;
            for (IndexType i = 0; i < TNumNodes; ++i)
                for (IndexType d = 0; d < TDim; ++d)
                    DN_DX_Deriv(i,d) = -DN_DX(iCoord / TDim,d) * DN_DX(i,iCoord % TDim);

            // Volume'
            double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

            // u * Grad(N)'
            array_1d< double, TNumNodes > DensityVelGradNDeriv;
            noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv,Velocity);

            // TauOne', TauTwo'
            double TauOneDeriv, TauTwoDeriv;
            this->CalculateStabilizationParametersDerivative(
                    TauOneDeriv,TauTwoDeriv,TauOne,TauTwo,VelNorm,ElemSize,Density,
                    Viscosity,DetJDeriv);

            BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
            array_1d< double, TFluidLocalSize > RHS;
            for (IndexType i = 0; i < TFluidLocalSize; i++)
            {
                RHS[i] = 0.0;
                for (IndexType j = 0; j < TFluidLocalSize; j++)
                    LHS(i,j) = 0.0;
            }

            for (IndexType i = 0; i < TNumNodes; ++i)
            {
                for (IndexType j = 0; j < TNumNodes; ++j)
                {
                    // Left-hand side matrix
                    double diag = 0.0;
                    double ddiag = 0.0;

                    // Convective term, v * (u * Grad(u))
                    diag += N[i] * DensityVelGradN[j];
                    ddiag += N[i] * DensityVelGradNDeriv[j];

                    // Stabilization, lsq convection
                    // (u * Grad(v)) * TauOne * (u * Grad(u))
                    diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];
                    ddiag += DensityVelGradNDeriv[i] * TauOne * DensityVelGradN[j]
                    + DensityVelGradN[i] * TauOneDeriv * DensityVelGradN[j]
                    + DensityVelGradN[i] * TauOne * DensityVelGradNDeriv[j];

                    for (IndexType m = 0; m < TDim; ++m)
                    {
                        for (IndexType n = 0; n < TDim; ++n)
                        {
                            // Stabilization, lsq divergence
                            // Div(v) * TauTwo * Div(u)
                            double valmn = DN_DX(i,m) * TauTwo * DN_DX(j,n);
                            double dvalmn = DN_DX_Deriv(i,m) * TauTwo * DN_DX(j,n)
                            + DN_DX(i,m) * TauTwoDeriv * DN_DX(j,n)
                            + DN_DX(i,m) * TauTwo * DN_DX_Deriv(j,n);

                            LHS(i*TBlockSize+m,j*TBlockSize+n) += VolumeDeriv*valmn
                            + Volume*dvalmn;
                        }
                        LHS(i*TBlockSize+m,j*TBlockSize+m) += VolumeDeriv * diag
                        + Volume * ddiag;

                        double valmp = 0.0;
                        double dvalmp = 0.0;
                        // Pressure term
                        // Div(v) * p
                        valmp -= DN_DX(i,m) * N[j];
                        dvalmp -= DN_DX_Deriv(i,m) * N[j];

                        // Stabilization, convection-pressure
                        // (u * Grad(v)) * TauOne * Grad(p)
                        valmp += TauOne * DensityVelGradN[i] * DN_DX(j,m);
                        dvalmp += TauOneDeriv * DensityVelGradN[i] * DN_DX(j,m)
                        + TauOne * DensityVelGradNDeriv[i] * DN_DX(j,m)
                        + TauOne * DensityVelGradN[i] * DN_DX_Deriv(j,m);

                        double valpn = 0.0;
                        double dvalpn = 0.0;
                        // Divergence term
                        // q * Div(u)
                        valpn += N[i] * DN_DX(j,m);
                        dvalpn += N[i] * DN_DX_Deriv(j,m);

                        // Stabilization, pressure-convection
                        // Grad(q) * TauOne * (u * Grad(u))
                        valpn += TauOne * DensityVelGradN[j] * DN_DX(i,m);
                        dvalpn += TauOneDeriv * DensityVelGradN[j] * DN_DX(i,m)
                        + TauOne * DensityVelGradNDeriv[j] * DN_DX(i,m)
                        + TauOne * DensityVelGradN[j] * DN_DX_Deriv(i,m);

                        LHS(i*TBlockSize+m,j*TBlockSize+TDim) += VolumeDeriv * valmp
                        + Volume * dvalmp;
                        LHS(i*TBlockSize+TDim,j*TBlockSize+m) += VolumeDeriv * valpn
                        + Volume * dvalpn;
                    }

                    double valpp = 0.0;
                    double dvalpp = 0.0;
                    // Stabilization, lsq pressure
                    // TauOne * Grad(q) * Grad(p)
                    for (IndexType d = 0; d < TDim; ++d)
                    {
                        valpp += DN_DX(i,d) * DN_DX(j,d) * TauOne;
                        dvalpp += DN_DX_Deriv(i,d) * DN_DX(j,d) * TauOne
                        + DN_DX(i,d) * DN_DX_Deriv(j,d) * TauOne
                        + DN_DX(i,d) * DN_DX(j,d) * TauOneDeriv;
                    }

                    LHS(i*TBlockSize+TDim,j*TBlockSize+TDim) += VolumeDeriv * valpp
                    + Volume * dvalpp;
                } // Node block columns

                // Right-hand side vector
                double DN_DX_BodyForce = 0.0;
                double DN_DX_BodyForceDeriv = 0.0;
                for (IndexType d = 0; d < TDim; ++d)
                {
                    DN_DX_BodyForce += DN_DX(i,d) * BodyForce[d];
                    DN_DX_BodyForceDeriv += DN_DX_Deriv(i,d) * BodyForce[d];
                }

                for (IndexType m = 0; m < TDim; ++m)
                {
                    double valm = 0.0;
                    double dvalm = 0.0;

                    // External body force
                    valm += N[i] * BodyForce[m];

                    // Stabilization, convection-BodyForce
                    // (u * Grad(v)) * TauOne * f
                    valm += TauOne * DensityVelGradN[i] * BodyForce[m];
                    dvalm += TauOneDeriv * DensityVelGradN[i] * BodyForce[m]
                    + TauOne * DensityVelGradNDeriv[i] * BodyForce[m];

                    RHS[i*TBlockSize+m] += VolumeDeriv * valm + Volume * dvalm;
                }

                double valp = TauOne * DN_DX_BodyForce;
                double dvalp = TauOneDeriv * DN_DX_BodyForce
                + TauOne * DN_DX_BodyForceDeriv;

                RHS[i*TBlockSize+TDim] += VolumeDeriv * valp + Volume * dvalp;
            } // Node block rows

            this->AddViscousTermDerivative(LHS,DN_DX,DN_DX_Deriv,Viscosity * Volume,
                    Viscosity * VolumeDeriv);

            // Assign the derivative of the residual w.r.t this coordinate to the
            // shape derivatives matrix.
            array_1d< double, TFluidLocalSize > ResidualDerivative;
            noalias(ResidualDerivative) = RHS - prod(LHS,FluidValues);
            for (IndexType k = 0; k < TFluidLocalSize; ++k)
                rShapeDerivativesMatrix(iCoord,k) = ResidualDerivative[k];
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the gradient matrix of the velocity.
     *
     * The row index corresponds to the velocity component and the column index to
     * the derivative.
     *
     * @param rGradVel velocity gradient matrix
     * @param rDN_DX shape functions' gradients
     */
    void CalculateVelocityGradient(BoundedMatrix<double, TDim, TDim>& rGradVel,
                                   const ShapeFunctionDerivativesType& rDN_DX)
    {
        GeometryType& rGeom = this->GetGeometry();
        // node 0
        const array_1d< double, 3 >& rVel = rGeom[0].FastGetSolutionStepValue(VELOCITY,0);
        for (IndexType m = 0; m < TDim; m++)
            for (IndexType n = 0; n < TDim; n++)
                rGradVel(m,n) = rDN_DX(0,n) * rVel[m];
        // node 1,2,...
        for (IndexType iNode = 1; iNode < TNumNodes; iNode++)
        {
            const array_1d< double, 3 >& rVel = rGeom[iNode].FastGetSolutionStepValue(VELOCITY,0);
            for (IndexType m = 0; m < TDim; m++)
                for (IndexType n = 0; n < TDim; n++)
                    rGradVel(m,n) += rDN_DX(iNode,n) * rVel[m];
        }
    }

    /**
     * @brief Returns the pressure gradient.
     *
     * @param rGradP pressure gradient
     * @param rDN_DX shape functions' gradients
     */
    void CalculatePressureGradient(array_1d<double, TDim>& rGradP,
                                   const ShapeFunctionDerivativesType& rDN_DX)
    {
        GeometryType& rGeom = this->GetGeometry();
        // node 0
        for (IndexType d = 0; d < TDim; d++)
            rGradP[d] = rDN_DX(0,d) * rGeom[0].FastGetSolutionStepValue(PRESSURE);
        // node 1,2,...
        for (IndexType iNode = 1; iNode < TNumNodes; iNode++)
            for (IndexType d = 0; d < TDim; ++d)
                rGradP[d] += rDN_DX(iNode,d) * rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
    }

    /**
     * @brief Returns the element's size.
     *
     * @param Volume the volume (area in 2D) of the element
     */
    double CalculateElementSize(const double Volume);

    /**
     * @brief Returns derivatives of determinant of Jacobian w.r.t coordinates.
     *
     * The derivative of the determinant of the Jacobian w.r.t the jth coordinate
     * of the ith node is stored at the index (i * TDim + j).
     *
     * This function is only valid when the determinant of the Jacobian is constant
     * over the element.
     *
     * @see Triangle2D3
     * @see Tetrahedra3D4
     */
    void CalculateDeterminantOfJacobianDerivatives(array_1d< double, TCoordLocalSize >& rDetJDerivatives);

    /**
     * @brief Returns the VMS stabilization parameters.
     *
     * @param rTauOne momentum stabilization parameter
     * @param rTauTwo divergence stabilization parameter
     * @param VelNorm Euclidean norm of the velocity
     * @param ElemSize size of this element
     * @param Density density of the fluid
     * @param Viscosity dynamic viscosity of the fluid
     */
    void CalculateStabilizationParameters(double& rTauOne,
                                          double& rTauTwo,
                                          double VelNorm,
                                          double ElemSize,
                                          double Density,
                                          double Viscosity,
                                          const ProcessInfo& rCurrentProcessInfo)
    {
        // assume DELTA_TIME < 0 !!!
        double tmp =-rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME];
        tmp += 2.0 * VelNorm / ElemSize;
        tmp *= Density;
        tmp += 4.0 * Viscosity / (ElemSize * ElemSize);
        rTauOne = 1.0 / tmp;
        rTauTwo = Viscosity + 0.5 * Density * ElemSize * VelNorm;
    }

    /**
     * @brief Returns stabilization parameters derived w.r.t a node's coordinate.
     *
     * @param rTauOneDeriv derivative of momentum stabilization parameter
     * @param rTauTwoDeriv derivative of divergence stabilization parameter
     * @param TauOne momentum stabilization parameter
     * @param TauTwo divergence stabilization parameter
     * @param VelNorm Euclidean norm of the velocity
     * @param ElemSize size of this element
     * @param Density density of the fluid
     * @param Viscosity dynamic viscosity of the fluid
     * @param DetJDeriv derivative of the determinant of the Jacobian
     */
    void CalculateStabilizationParametersDerivative(double& rTauOneDeriv,
                                                    double& rTauTwoDeriv,
                                                    double TauOne,
                                                    double TauTwo,
                                                    double VelNorm,
                                                    double ElemSize,
                                                    double Density,
                                                    double Viscosity,
                                                    double DetJDeriv);

    /**
     * @brief Returns a scalar variable at this integration point.
     *
     * @param rResult the value of the scalar variable at this integration point
     * @param rVariable the variable to be evaluated
     * @param rShapeFunc array of shape function values at this integration point
     */
    void EvaluateInPoint(double& rResult,
                         const Variable<double>& rVariable,
                         const array_1d<double, TNumNodes>& rShapeFunc,
                         IndexType step = 0)
    {
        GeometryType& rGeom = this->GetGeometry();
        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(rVariable, step);
        for (IndexType iNode = 1; iNode < TNumNodes; ++iNode)
        {
            rResult += rShapeFunc[iNode]
                      * rGeom[iNode].FastGetSolutionStepValue(rVariable, step);
        }
    }

    /**
     * @brief Returns a vector variable at this integration point.
     *
     * @param rResult the value of the vector variable at this integration point
     * @param rVariable the variable to be evaluated
     * @param rShapeFunc array of shape function values at this integration point
     */
    void EvaluateInPoint(array_1d<double, TDim>& rResult,
                         const Variable<array_1d<double, 3>>& rVariable,
                         const array_1d<double, TNumNodes>& rN,
                         IndexType step = 0)
    {
        GeometryType& rGeom = this->GetGeometry();
        const array_1d< double, 3 >& rNodalValue = rGeom[0].FastGetSolutionStepValue(rVariable, step);
        for (IndexType d=0; d < TDim; d++)
            rResult[d] = rN[0] * rNodalValue[d];
        for (IndexType iNode = 1; iNode < TNumNodes; ++iNode)
        {
            const array_1d< double, 3 >& rNodalValue = rGeom[iNode].FastGetSolutionStepValue(rVariable, step);
            for (IndexType d=0; d < TDim; d++)
                rResult[d] += rN[iNode] * rNodalValue[d];
        }
    }

    /**
     * @brief Adds viscous contributions to adjoint system matrix.
     *
     * @param rResult matrix to add viscous contributions to
     * @param rDN_DX shape functions' gradients
     * @param Weight integration weight including dynamic viscosity
     */
    void AddViscousTerm(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rResult,
                        const ShapeFunctionDerivativesType& rDN_DX,
                        const double Weight);

    /**
     * @brief Adds derivative of viscous term w.r.t a node's coordinate.
     *
     * @param rResult matrix to add viscous contributions to
     * @param rDN_DX shape functions' gradients
     * @param rDN_DX_Deriv shape functions' gradients derived w.r.t the coordinate
     * @param Weight integration weight including dynamic viscosity
     * @param WeightDeriv integration weight derived w.r.t the coordinate
     *
     * @see AddViscousTerm
     */
    void AddViscousTermDerivative(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rResult,
                                  const ShapeFunctionDerivativesType& rDN_DX,
                                  const ShapeFunctionDerivativesType& rDN_DX_Deriv,
                                  const double Weight,
                                  const double WeightDeriv);

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

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

        KRATOS_CATCH("");
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY;

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Unaccessible methods
    ///@{

    VMSAdjointElement& operator=(VMSAdjointElement const& rOther);

    VMSAdjointElement(VMSAdjointElement const& rOther);

    ///@}

};  // class VMSAdjointElement

///@} // Kratos classes

///@name Input and output
///@{

/// Defines an input stream operator that does nothing.
template<unsigned int TDim>
inline std::istream& operator >>(std::istream& rIStream,
        VMSAdjointElement<TDim>& rThis) {
    return rIStream;
}

/// Defines an output stream operator that prints element info.
template<unsigned int TDim>
inline std::ostream& operator <<(std::ostream& rOStream,
        const VMSAdjointElement<TDim>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // FluidDynamicsApplication group

}// namespace Kratos

#endif // KRATOS_VMS_ADJOINT_ELEMENT_H_INCLUDED defined

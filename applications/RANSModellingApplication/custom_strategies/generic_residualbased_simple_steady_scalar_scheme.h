//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Suneth Warnakulasuriya
//                 Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME)
#define KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME

// Project includes
#include "containers/array_1d.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class GenericResidualBasedSimpleSteadyScalarScheme
    : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericResidualBasedSimpleSteadyScalarScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef Element::GeometryType GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    GenericResidualBasedSimpleSteadyScalarScheme(const double RelaxationFactor)
        : Scheme<TSparseSpace, TDenseSpace>(), mRelaxationFactor(RelaxationFactor)
    {
        KRATOS_INFO("GenericResidualBasedSimpleSteadyScalarScheme")
            << " Using residual based simple steady scheme with relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";
    }

    ~GenericResidualBasedSimpleSteadyScalarScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& rA,
                TSystemVectorType& rDx,
                TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentElement->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        AddRelaxation(rCurrentElement->GetGeometry(), LHS_Contribution,
                      RHS_Contribution, CurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                LocalSystemMatrixType& LHS_Contribution,
                                                LocalSystemVectorType& RHS_Contribution,
                                                Condition::EquationIdVectorType& EquationId,
                                                ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentCondition->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        AddRelaxation(rCurrentCondition->GetGeometry(), LHS_Contribution,
                      RHS_Contribution, CurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        CalculateSystemContributions(rCurrentElement, LHS_Contribution, rRHS_Contribution,
                                     rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution,
                                               rRHS_Contribution, rEquationId,
                                               rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType& rA,
                                   TSystemVectorType& rDx,
                                   TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        for (typename ModelPart::NodesContainerType::iterator itNode =
                 rModelPart.NodesBegin();
             itNode != rModelPart.NodesEnd(); itNode++)
            itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;

        double output;
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        const int number_of_elements = rModelPart.NumberOfElements();
#pragma omp parallel for private(output)
        for (int i = 0; i < number_of_elements; i++)
        {
            ModelPart::ElementsContainerType::iterator it_elem =
                rModelPart.ElementsBegin() + i;
            it_elem->Calculate(NODAL_AREA, output, CurrentProcessInfo);
        }

        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Operators
    ///@{
    void AddRelaxation(const GeometryType& rGeometry,
                       LocalSystemMatrixType& LHS_Contribution,
                       LocalSystemVectorType& RHS_Contribution,
                       ProcessInfo& CurrentProcessInfo)
    {
        if (LHS_Contribution.size1() == 0)
            return;

        const unsigned int NumNodes = rGeometry.PointsNumber();
        const unsigned int Dimension = rGeometry.WorkingSpaceDimension();

        Matrix Mass;
        this->CalculateLumpedMassMatrix(rGeometry, Mass);

        unsigned int DofIndex = 0;
        for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
        {
            const array_1d<double, 3>& rVel =
                rGeometry[iNode].FastGetSolutionStepValue(VELOCITY, 0);
            const double Area = rGeometry[iNode].FastGetSolutionStepValue(NODAL_AREA, 0);
            double VelNorm = 0.0;
            for (unsigned int d = 0; d < Dimension; ++d)
                VelNorm += rVel[d] * rVel[d];
            VelNorm = sqrt(VelNorm);
            double LocalDt;
            if (VelNorm != 0.0)
                LocalDt = pow(Area, 0.5) / VelNorm;
            else
                LocalDt = 1.0;

            Mass(DofIndex, DofIndex) *= 1.0 / (mRelaxationFactor * LocalDt);
            DofIndex++;
        }
        noalias(LHS_Contribution) += Mass;
    }

    void CalculateLumpedMassMatrix(const GeometryType& rGeometry,
                                   LocalSystemMatrixType& rLumpedMass) const
    {
        const unsigned int number_of_nodes = rGeometry.PointsNumber();

        if (rLumpedMass.size1() != number_of_nodes)
        {
            rLumpedMass.resize(number_of_nodes, number_of_nodes, false);
        }

        noalias(rLumpedMass) = ZeroMatrix(number_of_nodes, number_of_nodes);

        const double size_fraction = rGeometry.DomainSize() / number_of_nodes;

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            rLumpedMass(i, i) = size_fraction;
        }
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    double mRelaxationFactor;
    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater =
        TSparseSpace::CreateDofUpdater();

    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME defined */

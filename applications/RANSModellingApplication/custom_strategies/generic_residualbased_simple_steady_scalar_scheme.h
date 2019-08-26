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

//debugging
#include "input_output/vtk_output.h"

#include "custom_strategies/relaxed_dof_updater.h"

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

    void InitializeSolutionStep(ModelPart &r_model_part,
                                        TSystemMatrixType &A,
                                        TSystemVectorType &Dx,
                                        TSystemVectorType &b) override
    {
        BaseType::InitializeSolutionStep(r_model_part,A,Dx,b);
        if (TSparseSpace::Size(mPreviousB) != TSparseSpace::Size(b)) {
            TSparseSpace::Resize(mPreviousB, TSparseSpace::Size(b));
        }
        TSparseSpace::SetToZero(mPreviousB);
        mIterationCounter = 0;
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& rA,
                TSystemVectorType& rDx,
                TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx, mRelaxationFactor);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const GeometryType& r_geometry = rCurrentElement->GetGeometry();
        const double local_delta_time = this->GetGeometryLocalTimeStep(r_geometry);

        rCurrentElement->SetValue(DELTA_TIME, local_delta_time);

        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentElement->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        // AddRelaxation(rCurrentElement->GetGeometry(), local_delta_time,
        //               LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                LocalSystemMatrixType& LHS_Contribution,
                                                LocalSystemVectorType& RHS_Contribution,
                                                Condition::EquationIdVectorType& EquationId,
                                                ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const GeometryType& r_geometry = rCurrentCondition->GetGeometry();
        const double local_delta_time = this->GetGeometryLocalTimeStep(r_geometry);

        rCurrentCondition->SetValue(DELTA_TIME, local_delta_time);

        rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentCondition->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        // AddRelaxation(rCurrentCondition->GetGeometry(), local_delta_time,
        //               LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

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
    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType& rA,
                                   TSystemVectorType& rDx,
                                   TSystemVectorType& rb) override
    {
        // rModelPart.GetProcessInfo()[STEP] += 1;
        // mVtkOutput->PrintOutput();
    }


    void InitializeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType& rA,
                                   TSystemVectorType& rDx,
                                   TSystemVectorType& rb) override
    {
        // Debugging output
        Parameters default_parameters = Parameters(R"(
            {
                "model_part_name"                    : "FluidModelPart",
                "file_format"                        : "ascii",
                "output_precision"                   : 7,
                "output_control_type"                : "step",
                "output_frequency"                   : 1.0,
                "output_sub_model_parts"             : false,
                "folder_name"                        : "VTK_Output",
                "custom_name_prefix"                 : "",
                "save_output_files_in_folder"        : true,
                "nodal_solution_step_data_variables" : ["VELOCITY", "PRESSURE", "KINEMATIC_VISCOSITY", "TURBULENT_KINETIC_ENERGY", "TURBULENT_ENERGY_DISSIPATION_RATE", "TURBULENT_VISCOSITY", "VISCOSITY", "RANS_Y_PLUS"],
                "nodal_data_value_variables"         : [],
                "element_data_value_variables"       : [],
                "condition_data_value_variables"     : [],
                "gauss_point_variables"              : []
            })" );

            default_parameters["model_part_name"].SetString(rModelPart.Name());

            if (mVtkOutput == nullptr)
                mVtkOutput = new VtkOutput(rModelPart, default_parameters);
//         KRATOS_TRY;

//         for (typename ModelPart::NodesContainerType::iterator itNode =
//                  rModelPart.NodesBegin();
//              itNode != rModelPart.NodesEnd(); itNode++)
//             itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;

//         double output;
//         ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
//         const int number_of_elements = rModelPart.NumberOfElements();
// #pragma omp parallel for private(output)
//         for (int i = 0; i < number_of_elements; i++)
//         {
//             ModelPart::ElementsContainerType::iterator it_elem =
//                 rModelPart.ElementsBegin() + i;
//             it_elem->Calculate(NODAL_AREA, output, CurrentProcessInfo);
//         }

//         rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

//         KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Operators
    ///@{
    void AddRelaxation(const GeometryType& rGeometry,
                       const double LocalDeltaTime,
                       LocalSystemMatrixType& LHS_Contribution,
                       LocalSystemVectorType& RHS_Contribution,
                       ProcessInfo& CurrentProcessInfo)
    {
        if (LHS_Contribution.size1() == 0)
            return;

        Matrix Mass;
        this->CalculateLumpedMassMatrix(rGeometry, Mass);

        const unsigned int NumNodes = rGeometry.PointsNumber();
        const unsigned int Dimension = rGeometry.WorkingSpaceDimension();

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
                LocalDt = pow(Area, 1.0 / double(Dimension)) / VelNorm;
            else
                LocalDt = 1.0;

            Mass(iNode, iNode) *= 1.0 / (mRelaxationFactor * LocalDt);
        }

        noalias(LHS_Contribution) += Mass;

        // noalias(LHS_Contribution) += Mass * (1.0 / (mRelaxationFactor * LocalDeltaTime));
    }

    void CalculateLumpedMassMatrix(const GeometryType& rGeometry,
                                   LocalSystemMatrixType& rLumpedMass) const
    {
        const unsigned int number_of_nodes = rGeometry.PointsNumber();

        if (rLumpedMass.size1() != number_of_nodes)
            rLumpedMass.resize(number_of_nodes, number_of_nodes, false);

        rLumpedMass.clear();

        const double size_fraction = rGeometry.DomainSize() / number_of_nodes;

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rLumpedMass(i, i) = size_fraction;
    }

    double GetGeometryLocalTimeStep(const GeometryType& rGeometry) const
    {
        const unsigned int number_of_nodes = rGeometry.PointsNumber();

        double geometry_local_dt = 0.0;

        for (unsigned int i_node = 0; i_node < number_of_nodes; i_node++)
        {
            const array_1d<double, 3>& r_velocity =
                rGeometry[i_node].FastGetSolutionStepValue(VELOCITY);
            const double area = rGeometry[i_node].FastGetSolutionStepValue(NODAL_AREA);
            double velocity_magnitude = norm_2(r_velocity);

            const double local_dt =
                (velocity_magnitude > std::numeric_limits<double>::epsilon())
                    ? std::pow(area, 0.5) / velocity_magnitude
                    : 1.0;
            geometry_local_dt += local_dt;
        }

        return geometry_local_dt / number_of_nodes;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    TSystemVectorType mPreviousB;

    double mPreviousRelaxationFactor;

    unsigned int mIterationCounter = 0;

    VtkOutput* mVtkOutput;

    typedef RelaxedDofUpdater<TSparseSpace> DofUpdaterType;
    typedef typename DofUpdaterType::UniquePointer DofUpdaterPointerType;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;
    // typename TSparseSpace::DofUpdaterPointerType mpDofUpdater =
    //     TSparseSpace::CreateDofUpdater();

    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME defined */

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#include <algorithm>
#include <exception>
#include <sstream>

#include "testing/testing.h"
#include "containers/model.h"

#include "includes/define.h"
#include "includes/shared_pointers.h"
#include "includes/model_part.h"
#include "containers/pointer_vector.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "utilities/indirect_scalar.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/sensitivity_builder.h"
#include "utilities/adjoint_extensions.h"


namespace Kratos
{
namespace Testing
{
namespace
{
const double AlphaBossak = -0.3;
namespace Base
{
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

struct PrimalResults
{
    KRATOS_CLASS_POINTER_DEFINITION(PrimalResults);
    virtual ~PrimalResults() {};
    virtual void StoreCurrentSolutionStep(const ModelPart& rModelPart) = 0;
    virtual void LoadCurrentSolutionStep(ModelPart& rModelPart) const = 0;
};

class PrimalStrategy : public SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrimalStrategy);

    PrimalStrategy(ModelPart& rModelPart, PrimalResults::Pointer pPrimalResults)
        : SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>(rModelPart),
          mpPrimalResults(pPrimalResults)
    {
        auto p_scheme =
            Kratos::make_shared<ResidualBasedBossakDisplacementScheme<SparseSpaceType, LocalSpaceType>>(
                AlphaBossak);
        auto p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        auto p_conv_criteria =
            Kratos::make_shared<ResidualCriteria<SparseSpaceType, LocalSpaceType>>(
                1e-10, 1e-13);
        p_conv_criteria->SetEchoLevel(0);
        mpSolver = Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
            rModelPart, p_scheme, p_linear_solver);
    }

    void Initialize() override
    {
        mpSolver->Initialize();
        mpSolver->SetEchoLevel(0);
    }

    double Solve() override
    {
        auto result = mpSolver->Solve();
        mpPrimalResults->StoreCurrentSolutionStep(this->GetModelPart());
        return result;
    }

private:
    PrimalResults::Pointer mpPrimalResults;
    SolvingStrategyType::Pointer mpSolver;
};

class AdjointStrategy : public SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointStrategy);

    AdjointStrategy(ModelPart& rModelPart,
                    Kratos::shared_ptr<PrimalResults> pPrimalResults,
                    AdjointResponseFunction::Pointer pResponseFunction)
        : SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>(rModelPart),
          mpPrimalResults(pPrimalResults)
    {
        auto p_linear_solver =
            Kratos::make_shared<SkylineLUCustomScalarSolver<SparseSpaceType, LocalSpaceType>>();
        auto scheme_settings = Parameters{R"({ "alpha_bossak": )" + std::to_string(AlphaBossak) + " }"};
        auto p_adjoint_scheme =
            Kratos::make_shared<ResidualBasedAdjointBossakScheme<SparseSpaceType, LocalSpaceType>>(
                scheme_settings, pResponseFunction);
        mpSolver =
            Kratos::make_shared<ResidualBasedLinearStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>>(
                rModelPart, p_adjoint_scheme, p_linear_solver);
    }

    void Initialize() override
    {
        mpSolver->Initialize();
        mpSolver->SetEchoLevel(0);
    }

    double Solve() override
    {
        mpPrimalResults->LoadCurrentSolutionStep(this->GetModelPart());
        return mpSolver->Solve();
    }

private:
    PrimalResults::Pointer mpPrimalResults;
    SolvingStrategyType::Pointer mpSolver;
};
}

namespace NonLinearSpringMassDamper
{
/**
 * @class PrimalElement
 * @brief A system of two mass-spring-dampers for testing a second-order ode.
 * @details Taken from L.F. Fernandez, D.A. Tortorelli, Semi-analytical
 * sensitivity analysis for nonlinear transient problems.
 *
 *  |                _____                 _____
 *  |---[ Damper ]--|mass1|---[ Damper ]--|mass2|
 *  |-----/\/\/\----|_____|-----/\/\/\----|_____|
 *  |
 *
 * Spring force: fe = x + stiffness * x^3
 * Damper force: fc = damping * x'
 *
 * Momentum equations:
 * mass1 * acc1 + 2 * damping * vel1 - damping * vel2 + 2 * disp1 - disp2 + stiffness * disp1^3 - stiffness * (disp2 - disp1)^3 = 0,
 * mass2 * acc2 - damping * vel1 + damping * vel2 - disp1 + disp2 + stiffness * (disp2 - disp1)^3 = 0.
 */
class PrimalElement : public Element
{
public:
    typedef Kratos::intrusive_ptr<PrimalElement> Pointer;
    typedef Kratos::unique_ptr<PrimalElement> UniquePointer;

    static Pointer Create(Node<3>::Pointer pNode1)
    {
        auto nodes = PointerVector<Node<3>>{};
        nodes.push_back(pNode1);
        return Kratos::make_intrusive<PrimalElement>(nodes);
    }

    PrimalElement(const NodesArrayType& ThisNodes)
        : Element(0, ThisNodes)
    {
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult.resize(1);
        rResult[0] = this->GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    }

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rElementalDofList.resize(1);
        rElementalDofList[0] = this->GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    }

    void GetValuesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
    }

    void GetFirstDerivativesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X, Step);
    }

    void GetSecondDerivativesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_X, Step);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix.resize(1, 1, false);
        rLeftHandSideMatrix(0, 0) = stiffness;
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector.resize(1, false);
        const double& x1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        rRightHandSideVector(0) = -(stiffness * x1);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        rMassMatrix.resize(1, 1, false);
        rMassMatrix(0, 0) = mass1;
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        rDampingMatrix.resize(1, 1, false);
        rDampingMatrix(0, 0) = damping;
    }

private:
    const double mass1 = 1.;
    const double stiffness = 1.;
    const double damping = 0.04;
};

class AdjointElement : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        ThisExtensions(Element* pElement) : mpElement{pElement}
        {
        }

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            if (rVector.size() != 1)
            {
                rVector.resize(1);
            }
            rVector[0] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_X, Step);
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            if (rVector.size() != 1)
            {
                rVector.resize(1);
            }
            rVector[0] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_X, Step);
        }

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            if (rVector.size() != 1)
            {
                rVector.resize(1);
            }
            rVector[0] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_X, Step);
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            if (rVariables.size() != 1)
            {
                rVariables.resize(1);
            }
            rVariables[0] = &ADJOINT_VECTOR_2;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            if (rVariables.size() != 1)
            {
                rVariables.resize(1);
            }
            rVariables[0] = &ADJOINT_VECTOR_3;
        }

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override
        {
            if (rVariables.size() != 1)
            {
                rVariables.resize(1);
            }
            rVariables[0] = &AUX_ADJOINT_VECTOR_1;
        }
    };

public:
    typedef Kratos::intrusive_ptr<AdjointElement> Pointer;
    typedef Kratos::unique_ptr<AdjointElement> UniquePointer;


    static Pointer Create(Node<3>::Pointer pNode1)
    {
        auto nodes = PointerVector<Node<3>>{};
        nodes.push_back(pNode1);
        return Kratos::make_intrusive<AdjointElement>(nodes);
    }

    AdjointElement(const NodesArrayType& ThisNodes)
        : Element(0, ThisNodes), mPrimalElement(ThisNodes)
    {
        SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult.resize(1);
        rResult[0] = this->GetGeometry()[0].GetDof(ADJOINT_VECTOR_1_X).EquationId();
    }

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rElementalDofList.resize(1);
        rElementalDofList[0] = this->GetGeometry()[0].pGetDof(ADJOINT_VECTOR_1_X);
    }

    void GetValuesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_VECTOR_1_X, Step);
    }

    void GetFirstDerivativesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_VECTOR_2_X, Step);
    }

    void GetSecondDerivativesVector(Vector& values, int Step = 0) const override
    {
        values.resize(1);
        values[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ADJOINT_VECTOR_3_X, Step);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalElement.CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
        KRATOS_CATCH("");
    }

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        if (rDesignVariable == SCALAR_SENSITIVITY)
        {
            rOutput.resize(1, 1, false);
            const double& x1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
            rOutput(0, 0) = -(x1);
        }
        else
        {
            KRATOS_ERROR << "Invalid variable: " << rDesignVariable << std::endl;
        }
        KRATOS_CATCH("");
    }

    ///@}

private:
    PrimalElement mPrimalElement;
};

class ResponseFunction : public AdjointResponseFunction
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

        ResponseFunction(ModelPart& rModelPart) : mrModelPart(rModelPart)
        {
        }

        void CalculateGradient(const Element& rAdjointElement,
                               const Matrix& rResidualGradient,
                               Vector& rResponseGradient,
                               const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(1, false);
            rResponseGradient(0) = rAdjointElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
        }

        void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                               const Matrix& rResidualGradient,
                                               Vector& rResponseGradient,
                                               const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(1, false);
            rResponseGradient(0) = rAdjointElement.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
        }

        void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                                const Matrix& rResidualGradient,
                                                Vector& rResponseGradient,
                                                const ProcessInfo& rProcessInfo) override
        {
            rResponseGradient.resize(1, false);
            rResponseGradient(0) = 0.;
        }

        void CalculatePartialSensitivity(Element& rAdjointElement,
                                         const Variable<double>& rVariable,
                                         const Matrix& rSensitivityMatrix,
                                         Vector& rSensitivityGradient,
                                         const ProcessInfo& rProcessInfo) override
        {
            rSensitivityGradient.resize(1, false);
            rSensitivityGradient(0) = 0.;
        }

        void CalculatePartialSensitivity(Element& rAdjointElement,
                                         const Variable<array_1d<double, 3>>& rVariable,
                                         const Matrix& rSensitivityMatrix,
                                         Vector& rSensitivityGradient,
                                         const ProcessInfo& rProcessInfo) override
        {
            rSensitivityGradient.resize(1, false);
            rSensitivityGradient(0) = 0.;
        }

        double CalculateValue(ModelPart& rModelPart) override
        {
            const double& x1 =
                mrModelPart.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT_X);
            const double& v1 = mrModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY_X);
            return 0.5*(x1 * x1 + v1 * v1);
        }

    private:
        ModelPart& mrModelPart;
};

struct PrimalResults : Base::PrimalResults
{
    std::vector<double> time;
    std::vector<double> x1;
    std::vector<double> v1;
    std::vector<double> a1;
    
    void StoreCurrentSolutionStep(const ModelPart& rModelPart) override
    {
        this->time.push_back(rModelPart.GetProcessInfo()[TIME]);
        auto& node1 = rModelPart.GetNode(1);
        this->x1.push_back(node1.FastGetSolutionStepValue(DISPLACEMENT_X));
        this->v1.push_back(node1.FastGetSolutionStepValue(VELOCITY_X));
        const double acc1 =
            (1. - AlphaBossak) * node1.FastGetSolutionStepValue(ACCELERATION_X) +
            AlphaBossak * node1.FastGetSolutionStepValue(ACCELERATION_X, 1);
        this->a1.push_back(acc1);
    }

    void LoadCurrentSolutionStep(ModelPart& rModelPart) const override
    {
        const double current_time = rModelPart.GetProcessInfo()[TIME];
        if (current_time < 1e-8 || current_time > this->time.back() + 1e-8)
        {
            std::stringstream ss;
            ss << "Adjoint time = " << current_time
               << " outside of primal solution time range!\n";
            throw std::runtime_error{ss.str()};
        }
        auto it = std::find_if(
            this->time.begin(), this->time.end(),
            [current_time](const double& t) -> bool { return current_time <= t; });
        std::size_t pos = it - this->time.begin();
        if (pos == this->time.size())
        {
            --pos;
        }
        else if (pos > 0)
        {
            const auto t0 = this->time.at(pos - 1);
            const auto t1 = this->time.at(pos);
            if (std::abs(current_time - t0) < 1e-8)
            {
                pos = pos - 1;
            }
            else if (!(std::abs(current_time - t1) < 1e-8))
            {
                std::stringstream ss;
                ss << "Adjoint time = " << current_time
                   << " does not match primal solution time!\n";
                throw std::runtime_error{ss.str()};
            }
        }
        auto& node1 = rModelPart.GetNode(1);
        node1.FastGetSolutionStepValue(DISPLACEMENT_X) = this->x1.at(pos);
        node1.FastGetSolutionStepValue(VELOCITY_X) = this->v1.at(pos);
        node1.FastGetSolutionStepValue(ACCELERATION_X) = this->a1.at(pos);
    }
};

void InitializePrimalModelPart(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 1);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.SetBufferSize(1);
    for (auto& r_node : rModelPart.Nodes())
    {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
    }
    auto p_element =
        PrimalElement::Create(rModelPart.pGetNode(1));
    rModelPart.AddElement(p_element);
    auto& node1 = rModelPart.GetNode(1);
    node1.FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
    node1.FastGetSolutionStepValue(VELOCITY_X) = -0.75;
}

void InitializeAdjointModelPart(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 1);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_VECTOR_2);
    rModelPart.AddNodalSolutionStepVariable(ADJOINT_VECTOR_3);
    rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_VECTOR_1);
    rModelPart.AddNodalSolutionStepVariable(SCALAR_SENSITIVITY);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.SetBufferSize(1);
    for (auto& r_node : rModelPart.Nodes())
    {
        r_node.AddDof(ADJOINT_VECTOR_1_X);
    }
    auto p_adjoint_element =
        AdjointElement::Create(rModelPart.pGetNode(1));
    rModelPart.AddElement(p_adjoint_element);
}

}

} // unnamed namespace

KRATOS_TEST_CASE_IN_SUITE(ResidualBasedAdjointBossak_TwoMassSpringDamperSystem, KratosCoreFastSuite)
{
    namespace Nlsmd = NonLinearSpringMassDamper;
    // Solve the primal problem.
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("test");

    Nlsmd::InitializePrimalModelPart(model_part);
    auto p_results_data = Kratos::make_shared<Nlsmd::PrimalResults>();
    Base::PrimalStrategy solver(model_part, p_results_data);
    solver.Initialize();
    const double end_time = 20;
    const double start_time = 0.;
    const std::size_t N = 20;
    const double delta_time = (end_time - start_time) / N;
    model_part.CloneTimeStep(start_time - delta_time);
    model_part.CloneTimeStep(start_time);

    std::cout << std::endl << model_part.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT_X) << std::endl; 
    
    for (double current_time = start_time; current_time < end_time;)
    {
        current_time += delta_time;
        model_part.CloneTimeStep(current_time);
        solver.Solve();
        KRATOS_WATCH(model_part.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT_X)) 
    }

    // Solve the adjoint problem.
    ModelPart& adjoint_model_part = current_model.CreateModelPart("test_adjoint");

    Nlsmd::InitializeAdjointModelPart(adjoint_model_part);
    auto p_response_function =
        Kratos::make_shared<Nlsmd::ResponseFunction>(adjoint_model_part);
    Base::AdjointStrategy adjoint_solver(adjoint_model_part, p_results_data, p_response_function);
    adjoint_solver.Initialize();
    SensitivityBuilder sensitivity_builder(
        Parameters{R"(
            {
                "element_data_value_sensitivity_variables": ["SCALAR_SENSITIVITY"],
                "build_mode": "integrate",
                "nodal_solution_step_sensitivity_calculation_is_thread_safe" : true
            })"},
        adjoint_model_part, p_response_function);
    sensitivity_builder.Initialize();
    adjoint_model_part.CloneTimeStep(end_time + 2. * delta_time);
    adjoint_model_part.CloneTimeStep(end_time + delta_time);
    for (double current_time = end_time + delta_time; current_time >= start_time + 1.5 * delta_time;)
    {
        current_time -= delta_time;
        adjoint_model_part.CloneTimeStep(current_time);
        adjoint_solver.Solve();
        KRATOS_WATCH(model_part.GetNode(1).FastGetSolutionStepValue(ADJOINT_VECTOR_1_X))
        sensitivity_builder.UpdateSensitivities();
    }
    
    // Check.
    const double adjoint_sensitivity = adjoint_model_part.Elements().front().GetValue(SCALAR_SENSITIVITY);
    KRATOS_WATCH(adjoint_sensitivity)
   
}
}
}
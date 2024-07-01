//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#if !defined(KRATOS_BDF2_TURBULENT_SCHEME_DEM_COUPLED_H_INCLUDED )
#define  KRATOS_BDF2_TURBULENT_SCHEME_DEM_COUPLED_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "includes/define.h"
#include "includes/dof.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_strategies/schemes/bdf2_turbulent_scheme.h"
#include "swimming_dem_application_variables.h"
#include "linear_solvers/amgcl_solver.h"


namespace Kratos
{
///@addtogroup SwimmingDEMApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// A scheme for BDF2 time integration.
/**
 */
template<class TSparseSpace,class TDenseSpace>
class BDF2TurbulentSchemeDEMCoupled : public BDF2TurbulentScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BDF2TurbulentSchemeDEMCoupled
    KRATOS_CLASS_POINTER_DEFINITION(BDF2TurbulentSchemeDEMCoupled);
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;
    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixType MatrixType;
    typedef typename TSparseSpace::VectorType VectorType;


    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Dof<TDataType> TDofType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef CoordinateTransformationUtils<LocalSystemMatrixType, LocalSystemVectorType, double> RotationToolType;
    typedef typename RotationToolType::UniquePointer RotationToolPointerType;
    typedef Element::GeometryType GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BDF2TurbulentSchemeDEMCoupled()
    : BDF2TurbulentScheme<TSparseSpace, TDenseSpace>()
    , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor to use the formulation combined with a turbulence model.
    /**
     * The turbulence model is assumed to be implemented as a Kratos::Process.
     * The model's Execute() method wil be called at the start of each
     * non-linear iteration.
     * @param pTurbulenceModel pointer to the turbulence model
     */
    BDF2TurbulentSchemeDEMCoupled(Process::Pointer pTurbulenceModel)
        : BDF2TurbulentScheme<TSparseSpace, TDenseSpace>()
        , mpTurbulenceModel(pTurbulenceModel)
        , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor for periodic boundary conditions.
    /**
     * @param rPeriodicVar the variable used to store periodic pair indices.
     */
    BDF2TurbulentSchemeDEMCoupled(const Kratos::Variable<int>& rPeriodicVar)
        : BDF2TurbulentScheme<TSparseSpace, TDenseSpace>()
        , mrPeriodicIdVar(rPeriodicVar)
    {}


    /// Destructor.
    ~BDF2TurbulentSchemeDEMCoupled() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Set the time iteration coefficients
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        ProcessInfo CurrentProcessInfo = rModelPart.GetProcessInfo();

        this->SetTimeCoefficients(rModelPart.GetProcessInfo());

        // Base function initializes elements and conditions
        BaseType::InitializeSolutionStep(rModelPart,A,Dx,b);

        // Recalculate mesh velocity (to account for variable time step)
        const double tol = 1.0e-12;
        const double Dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double OldDt = rModelPart.GetProcessInfo().GetPreviousSolutionStepInfo(1)[DELTA_TIME];
        if(std::abs(Dt - OldDt) > tol) {
            const int n_nodes = rModelPart.NumberOfNodes();
            const Vector& BDFcoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];

#pragma omp parallel for
            for(int i_node = 0; i_node < n_nodes; ++i_node) {
                auto it_node = rModelPart.NodesBegin() + i_node;
                auto& rMeshVel = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
                const auto& rDisp0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                const auto& rDisp1 = it_node->FastGetSolutionStepValue(DISPLACEMENT,1);
                const auto& rDisp2 = it_node->FastGetSolutionStepValue(DISPLACEMENT,2);
                rMeshVel = BDFcoefs[0] * rDisp0 + BDFcoefs[1] * rDisp1 + BDFcoefs[2] * rDisp2;
            }
        }
        this->UpdateFluidFraction(rModelPart, CurrentProcessInfo);
    }

    void SetTimeCoefficients(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        //calculate the BDF coefficients
        double OldDt;
        double Dt = rCurrentProcessInfo[DELTA_TIME];
        double step = rCurrentProcessInfo[STEP];
        // Initialization of the previous delta time at the beginning of the simulation (when using adaptive delta time)
        if (rCurrentProcessInfo[MANUFACTURED] && step < 2){
            OldDt = rCurrentProcessInfo[DELTA_TIME];
        }
        else {
            OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];
        }

        double Rho = OldDt / Dt;
        double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

        Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
        BDFcoeffs.resize(3, false);

        BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
        BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
        BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)

        KRATOS_CATCH("");
    }

    void FullProjection(ModelPart& rModelPart)
    {
        const ProcessInfo ProcessInfo = rModelPart.GetProcessInfo();

        KRATOS_INFO_IF("ResidualBasedSimpleSteadyScheme", rModelPart.GetCommunicator().MyPID() == 0)<< "Computing OSS projections" << std::endl;

        const int number_of_nodes = rModelPart.NumberOfNodes();
        const int number_of_elements = rModelPart.NumberOfElements();

        unsigned int dimension = ProcessInfo[DOMAIN_SIZE];
        VectorType MassProjectionRHS = ZeroVector(number_of_nodes);
        VectorType MomentumProjectionRHS = ZeroVector(number_of_nodes*dimension);
        if (mMassMatrixAlreadyComputed == false){
            noalias(mGlobalDivProjMassMatrix) = ZeroMatrix(number_of_nodes,number_of_nodes);
            noalias(mGlobalAdvProjMassMatrix) = ZeroMatrix(number_of_nodes*dimension,number_of_nodes*dimension);
            #pragma omp for schedule(guided, 512)
            for (int e = 0; e < number_of_elements; e++){
            ModelPart::ElementsContainerType::iterator it_elem = rModelPart.ElementsBegin() + e;
            GeometryType r_geometry = it_elem->GetGeometry();
            unsigned int NumNodes = r_geometry.PointsNumber();
            GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
            GeometryType::IntegrationPointsArrayType r_integrations_points = r_geometry.IntegrationPoints( integration_method );
            auto r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
            Vector detJ_vector(r_number_integration_points);
            r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
            Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
            for (unsigned int g = 0; g < r_number_integration_points; g++){
                double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
                for (unsigned int i = 0; i < NumNodes; ++i){
                    for (unsigned int j= 0; j < NumNodes; ++j){
                        for (unsigned int d = 0; d < dimension; ++d){
                            mGlobalAdvProjMassMatrix(dimension*(r_geometry[i].Id()-1)+d,dimension*(r_geometry[j].Id()-1)+d) += Weight * NContainer(g,i) * NContainer(g,j);
                        }
                    mGlobalDivProjMassMatrix(r_geometry[i].Id()-1,r_geometry[j].Id()-1) += Weight * NContainer(g,i) * NContainer(g,j);
                    }
                }
            }
        }
        mMassMatrixAlreadyComputed = true;
        }

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++) {
        ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
        noalias(it_node->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
        it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        }
        array_1d<double, 3 > output;

        #pragma omp parallel for private(output)
        for (int i = 0; i < number_of_elements; i++) {
            ModelPart::ElementIterator it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(ADVPROJ,output,ProcessInfo);
        }

        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++) {
            ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
            MassProjectionRHS[i] += it_node->FastGetSolutionStepValue(DIVPROJ);
            array_1d<double,3>& AdvProj = it_node->FastGetSolutionStepValue(ADVPROJ);
            unsigned int row = i*dimension;
            for (unsigned int d = 0; d < dimension; ++d)
                MomentumProjectionRHS[row+d] += AdvProj[d];
        }

        VectorType MomProj = ZeroVector(number_of_nodes*dimension);
        VectorType MassProj = ZeroVector(number_of_nodes);

        AMGCLSolver<TSparseSpace, TDenseSpace > LinearSolver;
        LinearSolver.Solve(mGlobalAdvProjMassMatrix, MomProj, MomentumProjectionRHS);
        LinearSolver.Solve(mGlobalDivProjMassMatrix, MassProj, MassProjectionRHS);

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++) {
            ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;

            array_1d<double,3>& MomentumProjection = it_node->FastGetSolutionStepValue(ADVPROJ);
            unsigned int row = i*dimension;
            for (unsigned int d = 0; d < dimension; ++d)
                MomentumProjection[d] = MomProj[row+d];
            it_node->FastGetSolutionStepValue(DIVPROJ) = MassProj[i];
        }
    }

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        BaseType::InitializeNonLinIteration(rModelPart, A, Dx, b);

        if (mpTurbulenceModel != 0) mpTurbulenceModel->Execute();

        const ProcessInfo ProcessInfo = rModelPart.GetProcessInfo();

        //if orthogonal subscales are computed
        if (ProcessInfo[OSS_SWITCH] == 1.0)
        {
            //this->FullProjection(rModelPart);
            this->LumpedProjection(rModelPart);
        }

    }

    void FinalizeNonLinIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {

        BaseType::FinalizeNonLinIteration(rModelPart, A, Dx, b);

    }

    /// Store the iteration results as solution step variables and update acceleration after a Newton-Raphson iteration.
    /**
     * @param rModelPart fluid ModelPart
     * @param rDofSet DofSet containing the Newton-Raphson system degrees of freedom.
     * @param A Newton-Raphson system matrix (unused)
     * @param Dx Newton-Raphson iteration solution
     * @param b Newton-Raphson right hand side (unused)
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        double alpha = rModelPart.GetProcessInfo()[RELAXATION_ALPHA];

        TSparseSpace::InplaceMult(Dx, alpha);

        BDF2TurbulentScheme<TSparseSpace, TDenseSpace>::Update(rModelPart,rDofSet,A,Dx,b);

        KRATOS_CATCH("");
    }

    void UpdateFluidFraction(
        ModelPart& r_model_part,
        ProcessInfo& r_current_process_info)
    {
        BDF2TurbulentScheme<TSparseSpace, TDenseSpace>::SetTimeCoefficients(r_current_process_info);
        const Vector& BDFcoefs = r_current_process_info[BDF_COEFFICIENTS];
        double step = r_current_process_info[STEP];

        block_for_each(r_model_part.Nodes(), [&](Node& rNode)
        {
            double& fluid_fraction_0 = rNode.FastGetSolutionStepValue(FLUID_FRACTION);
            double& fluid_fraction_1 = rNode.FastGetSolutionStepValue(FLUID_FRACTION_OLD);
            double& fluid_fraction_2 = rNode.FastGetSolutionStepValue(FLUID_FRACTION_OLD_2);

            // This condition is needed to avoid large time variation of the porosity at the beginning of the simulation that can induce fluid instabilities
            // if (step <= 2){
            //     fluid_fraction_2 = fluid_fraction_0;
            //     fluid_fraction_1 = fluid_fraction_0;

            // }
            // else{
            rNode.FastGetSolutionStepValue(FLUID_FRACTION_RATE) = BDFcoefs[0] * fluid_fraction_0 + BDFcoefs[1] * fluid_fraction_1 + BDFcoefs[2] * fluid_fraction_2;
            //}

            rNode.GetSolutionStepValue(FLUID_FRACTION_OLD_2) = rNode.GetSolutionStepValue(FLUID_FRACTION_OLD);
            rNode.GetSolutionStepValue(FLUID_FRACTION_OLD) = rNode.GetSolutionStepValue(FLUID_FRACTION);
        });
    }

    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        BDF2TurbulentScheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(r_model_part, A, Dx, b);
        KRATOS_CATCH("")
    }

    /// Free memory allocated by this object.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "BDF2TurbulentSchemeDEMCoupled";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    /// Pointer to a turbulence model
    Process::Pointer mpTurbulenceModel = nullptr;

    RotationToolPointerType mpRotationTool = nullptr;

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    const Kratos::Variable<int>& mrPeriodicIdVar;

    MatrixType mGlobalDivProjMassMatrix;
    MatrixType mGlobalAdvProjMassMatrix;
    bool mMassMatrixAlreadyComputed = false;

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BDF2TurbulentSchemeDEMCoupled & operator=(BDF2TurbulentSchemeDEMCoupled const& rOther)
    {}

    /// Copy constructor.
    BDF2TurbulentSchemeDEMCoupled(BDF2TurbulentSchemeDEMCoupled const& rOther)
    {}

    ///@}

}; // Class BDF2TurbulentSchemeDEMCoupled

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TSparseSpace,class TDenseSpace>
inline std::istream& operator >>(std::istream& rIStream,BDF2TurbulentSchemeDEMCoupled<TSparseSpace,TDenseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpace,class TDenseSpace>
inline std::ostream& operator <<(std::ostream& rOStream,const BDF2TurbulentSchemeDEMCoupled<TSparseSpace,TDenseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_BDF2_TURBULENT_SCHEME_DEM_COUPLED_H_INCLUDED  defined

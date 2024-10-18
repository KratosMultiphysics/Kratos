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

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Dof<TDataType> TDofType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef CoordinateTransformationUtils<LocalSystemMatrixType, LocalSystemVectorType, double> RotationToolType;
    typedef typename RotationToolType::UniquePointer RotationToolPointerType;

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
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // Initialize containers
        const int n_nodes = rModelPart.NumberOfNodes();
        const int n_elems = rModelPart.NumberOfElements();
        const array_1d<double,3> zero_vect = ZeroVector(3);
#pragma omp parallel for firstprivate(zero_vect)
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto ind = rModelPart.NodesBegin() + i_node;
            noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = zero_vect; // "x"
            ind->FastGetSolutionStepValue(DIVPROJ) = 0.0; // "x"
            ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0; // "Ml"
        }

        // Newton-Raphson parameters
        const double RelTol = rModelPart.GetProcessInfo()[RELAXATION_ALPHA] * 1e-4 * rModelPart.NumberOfNodes();
        const double AbsTol = rModelPart.GetProcessInfo()[RELAXATION_ALPHA] * 1e-6 * rModelPart.NumberOfNodes();
        const unsigned int MaxIter = 100;

        // iteration variables
        unsigned int iter = 0;
        array_1d<double,3> dMomProj = zero_vect;
        double dMassProj = 0.0;

        double RelMomErr = 1000.0 * RelTol;
        double RelMassErr = 1000.0 * RelTol;
        double AbsMomErr = 1000.0 * AbsTol;
        double AbsMassErr = 1000.0 * AbsTol;

        while( ( (AbsMomErr > AbsTol && RelMomErr > RelTol) || (AbsMassErr > AbsTol && RelMassErr > RelTol) ) && iter < MaxIter)
        {
            // Reinitialize RHS
#pragma omp parallel for firstprivate(zero_vect)
            for (int i_node = 0; i_node < n_nodes; ++i_node)
            {
                auto ind = rModelPart.NodesBegin() + i_node;
                noalias(ind->GetValue(ADVPROJ)) = zero_vect; // "b"
                ind->GetValue(DIVPROJ) = 0.0; // "b"
                ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0; // Reset because Calculate will overwrite it
            }

            // Reinitialize errors
            RelMomErr = 0.0;
            RelMassErr = 0.0;
            AbsMomErr = 0.0;
            AbsMassErr = 0.0;

            // Compute new values
            array_1d<double, 3 > output;
#pragma omp parallel for private(output)
            for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
                auto it_elem = rModelPart.ElementsBegin() + i_elem;
                it_elem->Calculate(SUBSCALE_VELOCITY, output, rCurrentProcessInfo);
            }

            rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);

            // Update iteration variables
#pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; ++i_node) {
                auto ind = rModelPart.NodesBegin() + i_node;
                const double Area = ind->FastGetSolutionStepValue(NODAL_AREA); // Ml dx = b - Mc x
                dMomProj = ind->GetValue(ADVPROJ) / Area;
                dMassProj = ind->GetValue(DIVPROJ) / Area;

                RelMomErr += std::sqrt(std::pow(dMomProj[0],2) + std::pow(dMomProj[1],2) + std::pow(dMomProj[2],2));
                RelMassErr += std::fabs(dMassProj);

                auto& rMomRHS = ind->FastGetSolutionStepValue(ADVPROJ);
                double& rMassRHS = ind->FastGetSolutionStepValue(DIVPROJ);
                rMomRHS += dMomProj;
                rMassRHS += dMassProj;

                AbsMomErr += std::sqrt(std::pow(rMomRHS[0],2) + std::pow(rMomRHS[1],2) + std::pow(rMomRHS[2],2));
                AbsMassErr += std::fabs(rMassRHS);
            }

            if(AbsMomErr > 1e-10)
                RelMomErr /= AbsMomErr;
            else // If residual is close to zero, force absolute convergence to avoid division by zero errors
                RelMomErr = 1000.0;

            if(AbsMassErr > 1e-10)
                RelMassErr /= AbsMassErr;
            else
                RelMassErr = 1000.0;

            iter++;
        }

        KRATOS_INFO("BDF2TurbulentSchemeDEMCoupled") << "Performed OSS Projection in " << iter << " iterations" << std::endl;
    }

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (mpTurbulenceModel != 0) mpTurbulenceModel->Execute();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //if orthogonal subscales are computed
        if (CurrentProcessInfo[OSS_SWITCH] == 1.0)
        {
            //this->LumpedProjection(rModelPart);
            this->FullProjection(rModelPart);
        }

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {

        BaseType::FinalizeNonLinIteration(rModelPart, A, Dx, b);

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
            if (step <= 2){
                fluid_fraction_2 = fluid_fraction_0;
                fluid_fraction_1 = fluid_fraction_0;
            }

            rNode.FastGetSolutionStepValue(FLUID_FRACTION_RATE) = BDFcoefs[0] * fluid_fraction_0 + BDFcoefs[1] * fluid_fraction_1 + BDFcoefs[2] * fluid_fraction_2;

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

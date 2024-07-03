//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua, https://github.com/jgonzalezusua
//

#if !defined(KRATOS_BDF2_TURBULENT_SCHEME_DEM_COUPLED_H_INCLUDED )
#define  KRATOS_BDF2_TURBULENT_SCHEME_DEM_COUPLED_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <omp.h>
// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "utilities/parallel_utilities.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_strategies/schemes/bdf2_turbulent_scheme.h"
#include "swimming_dem_application_variables.h"
#include "linear_solvers/amgcl_solver.h"
#include "custom_utilities/binbased_DEM_fluid_coupled_mapping.h"


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
    //typedef BinBasedDEMFluidCoupledMapping<3, SphericParticle> BinBasedDEMFluidCoupledMapping3D;
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
    typedef GeometryData::SizeType  SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BDF2TurbulentSchemeDEMCoupled()
    : BDF2TurbulentScheme<TSparseSpace, TDenseSpace>()
    , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}


    BDF2TurbulentSchemeDEMCoupled(BinBasedDEMFluidCoupledMapping<3, SphericParticle>& rProjector)
    : BDF2TurbulentScheme<TSparseSpace, TDenseSpace>()
    , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    , mpProjector(&rProjector)
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
        //this->UpdateFluidFraction(rModelPart, CurrentProcessInfo);
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

    void AssembleMassMatrix(TSystemMatrixType& rA,
                          Matrix& rElementContribution,
                          std::vector<std::size_t>& rEquationId)
  {
    const SizeType local_size = rElementContribution.size1();

    for (unsigned int i_local = 0; i_local < local_size; i_local++)
    {
      unsigned int i_global = rEquationId[i_local];
      double* values_vector = rA.value_data().begin();
      std::size_t* index1_vector = rA.index1_data().begin();
      std::size_t* index2_vector = rA.index2_data().begin();
      size_t left_limit = index1_vector[i_global];
      //find the first entry
      size_t last_pos;
      unsigned int position = left_limit;
      while(rEquationId[0] != index2_vector[position])position++;

      last_pos = position;
      size_t last_found = rEquationId[0];
      double& r_a = values_vector[last_pos];
      const double& v_a = rElementContribution(i_local,0);
      AtomicAdd(r_a,  v_a);
      size_t posit = 0;
      unsigned int pos;
      for (unsigned int j=1; j<rEquationId.size(); j++) {
          unsigned int id_to_find = rEquationId[j];
          if(id_to_find > last_found) {
            pos = last_pos+1;
            while(id_to_find != index2_vector[pos]) pos++;
            posit = pos;
          } else if(id_to_find < last_found) {
            pos = last_pos-1;
            while(id_to_find != index2_vector[pos]) pos--;
            posit = pos;
          } else {
              posit = last_pos;
          }
          double& r = values_vector[posit];
          const double& v = rElementContribution(i_local,j);
          AtomicAdd(r,  v);
          last_found = id_to_find;
          last_pos = posit;
      }
    }
  }

  void ConstructMassMatrixStructure(TSystemMatrixType& rA,
                                    ModelPart& rModelPart,
                                    const SizeType rEquationSystemSize)
  {

      std::vector<std::unordered_set<IndexType> > indices(rEquationSystemSize);

      block_for_each(indices, [](std::unordered_set<IndexType>& rIndices){
          rIndices.reserve(40);
      });
      const int number_of_elements = rModelPart.NumberOfElements();
      std::vector<std::size_t> ids(3, 0);

      #pragma omp parallel firstprivate(ids,number_of_elements)
      {
          ProcessInfo& ProcessInfo = rModelPart.GetProcessInfo();
          // We repeat the same declaration for each thead
          std::vector<std::unordered_set<IndexType> > temp_indexes(rEquationSystemSize);

          #pragma omp for
          for (int index = 0; index < static_cast<int>(rEquationSystemSize); ++index)
              temp_indexes[index].reserve(30);

          // Element initial iterator
          const auto it_elem_begin = rModelPart.ElementsBegin();

          // We iterate over the elements
          #pragma omp for schedule(guided, 512) nowait
          for (int i_elem = 0; i_elem<number_of_elements; ++i_elem) {
              auto it_elem = it_elem_begin + i_elem;
              this->EquationId( *it_elem, ids, ProcessInfo);

              for (auto& id_i : ids) {
                  if (id_i < rEquationSystemSize) {
                      auto& row_indices = temp_indexes[id_i];
                      for (auto& id_j : ids)
                          if (id_j < rEquationSystemSize)
                              row_indices.insert(id_j);
                  }
              }
          }
          // Merging all the temporal indexes
          #pragma omp critical
          {
              for (int i = 0; i < static_cast<int>(temp_indexes.size()); ++i) {
                  indices[i].insert(temp_indexes[i].begin(), temp_indexes[i].end());
              }
          }
      }

      // Count the row sizes
      SizeType nnz = 0;
      for (IndexType i = 0; i < indices.size(); ++i)
          nnz += indices[i].size();

      rA = TSystemMatrixType(indices.size(), indices.size(), nnz);

      double* Avalues = rA.value_data().begin();
      std::size_t* Arow_indices = rA.index1_data().begin();
      std::size_t* Acol_indices = rA.index2_data().begin();

      // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
      Arow_indices[0] = 0;
      for (IndexType i = 0; i < rA.size1(); ++i)
          Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

      IndexPartition<std::size_t>(rA.size1()).for_each([&](std::size_t Index){
          const IndexType row_begin = Arow_indices[Index];
          const IndexType row_end = Arow_indices[Index + 1];
          IndexType k = row_begin;
          for (auto it = indices[Index].begin(); it != indices[Index].end(); ++it) {
              Acol_indices[k] = *it;
              Avalues[k] = 0.0;
              ++k;
          }

          std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
      });

      mGlobalProjMassMatrix.set_filled(indices.size() + 1, nnz);
  }

    void FullProjection(ModelPart& rModelPart)
    {
        const auto timer_projector = BuiltinTimer();
        KRATOS_INFO_IF("ResidualBasedBDF2SchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0)<< "Computing OSS projections" << std::endl;
        ProcessInfo ProcessInfo = rModelPart.GetProcessInfo();
        const int number_of_nodes = rModelPart.NumberOfNodes();
        const int number_of_elements = rModelPart.NumberOfElements();

        Timer::Start("FullProjection");
        const unsigned int dimension = ProcessInfo[DOMAIN_SIZE];
        const SizeType system_size = number_of_nodes*(dimension+1);
        VectorType ProjectionRHS = ZeroVector(system_size);
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        if (mMassMatrixAlreadyComputed == false){

            this->ConstructMassMatrixStructure(mGlobalProjMassMatrix,rModelPart,system_size);
            Matrix LocalLHS;
            std::vector<std::size_t> equation_id;
            #pragma omp parallel firstprivate(number_of_elements, LocalLHS, equation_id )
            {
            # pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < number_of_elements; k++){
            ModelPart::ElementsContainerType::iterator it_elem = el_begin + k;
            it_elem->EquationIdVector(equation_id, ProcessInfo);

            it_elem->Calculate(CONSISTENT_MASS_MATRIX,LocalLHS,ProcessInfo);

            this->AssembleMassMatrix(mGlobalProjMassMatrix,LocalLHS,equation_id);
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
            array_1d<double,3>& AdvProj = it_node->FastGetSolutionStepValue(ADVPROJ);
            unsigned int row = i*(dimension+1);
            for (unsigned int d = 0; d < dimension; ++d)
                ProjectionRHS[row+d] += AdvProj[d];
            ProjectionRHS[row+dimension] += it_node->FastGetSolutionStepValue(DIVPROJ);
        }

        VectorType Proj = ZeroVector(system_size);

        AMGCLSolver<TSparseSpace, TDenseSpace > LinearSolver;
        LinearSolver.Solve(mGlobalProjMassMatrix, Proj, ProjectionRHS);

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++) {
            ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;

            array_1d<double,3>& MomentumProjection = it_node->FastGetSolutionStepValue(ADVPROJ);
            unsigned int row = i*(dimension+1);
            for (unsigned int d = 0; d < dimension; ++d)
                MomentumProjection[d] = Proj[row+d];
            it_node->FastGetSolutionStepValue(DIVPROJ) = Proj[row+dimension];
        }
        Timer::Stop("FullProjection");
        KRATOS_INFO_IF("ResidualBasedBDF2SchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0) << "Residual projection time: " << timer_projector << std::endl;
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
            this->FullProjection(rModelPart);
            //this->LumpedProjection(rModelPart);
        }

    }

    void FinalizeNonLinIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        BaseType::FinalizeNonLinIteration(rModelPart, A, Dx, b);
        //mpProjector->UpdateHydrodynamicForces(rModelPart);

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

            rNode.FastGetSolutionStepValue(FLUID_FRACTION_RATE) = BDFcoefs[0] * fluid_fraction_0 + BDFcoefs[1] * fluid_fraction_1 + BDFcoefs[2] * fluid_fraction_2;

            fluid_fraction_2 = fluid_fraction_1;
            fluid_fraction_1 = fluid_fraction_0;
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
    BinBasedDEMFluidCoupledMapping<3,SphericParticle>::Pointer mpProjector = nullptr;
    MatrixType mGlobalProjMassMatrix;
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

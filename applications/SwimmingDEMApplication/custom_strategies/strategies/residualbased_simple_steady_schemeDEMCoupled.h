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


#if !defined(KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME_DEM_COUPLED )
#define  KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME_DEM_COUPLED


// External includes
#include "boost/smart_ptr.hpp"
#include <omp.h>
// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "linear_solvers/amgcl_solver.h"
#include "utilities/builtin_timer.h"


// Applications includes
#include "fluid_dynamics_application_variables.h"
#include "custom_strategies/schemes/residualbased_simple_steady_scheme.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TSparseSpace, class TDenseSpace >
class ResidualBasedSimpleSteadySchemeDEMCoupled : public ResidualBasedSimpleSteadyScheme<TSparseSpace, TDenseSpace> {
public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedSimpleSteadySchemeDEMCoupled);

  typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

  typedef LinearSolver<TSparseSpace,TDenseSpace> SolverType;

  //typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,SolverType> BuilderType;

  //typedef  SolverType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename TSparseSpace::MatrixType MatrixType;
  typedef typename TSparseSpace::VectorType VectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef Element::GeometryType GeometryType;

  typedef GeometryData::SizeType  SizeType;

  ///@}
  ///@name Life Cycle
  ///@{

  ResidualBasedSimpleSteadySchemeDEMCoupled(double VelocityRelaxationFactor,
                                  double PressureRelaxationFactor,
                                  unsigned int DomainSize)
      : ResidualBasedSimpleSteadyScheme<TSparseSpace, TDenseSpace>(VelocityRelaxationFactor, PressureRelaxationFactor, DomainSize),
        mVelocityRelaxationFactor(VelocityRelaxationFactor),
        mPressureRelaxationFactor(PressureRelaxationFactor),
        mRotationTool(DomainSize,DomainSize+1,SLIP)
  {}

  ResidualBasedSimpleSteadySchemeDEMCoupled(
      double VelocityRelaxationFactor,
      double PressureRelaxationFactor,
      unsigned int DomainSize,
      Process::Pointer pTurbulenceModel)
      : ResidualBasedSimpleSteadyScheme<TSparseSpace, TDenseSpace>(VelocityRelaxationFactor, PressureRelaxationFactor, DomainSize, pTurbulenceModel),
        mVelocityRelaxationFactor(VelocityRelaxationFactor),
        mPressureRelaxationFactor(PressureRelaxationFactor),
        mRotationTool(DomainSize,DomainSize+1,SLIP),
        mpTurbulenceModel(pTurbulenceModel)
  {}

  ~ResidualBasedSimpleSteadySchemeDEMCoupled() override {}


  ///@}
  ///@name Operators
  ///@{

  double GetVelocityRelaxationFactor() const
  {
    return mVelocityRelaxationFactor;
  }

  void SetVelocityRelaxationFactor(double factor)
  {
    mVelocityRelaxationFactor = factor;
  }

  double GetPressureRelaxationFactor() const
  {
    return mPressureRelaxationFactor;
  }

  void SetPressureRelaxationFactor(double factor)
  {
    mPressureRelaxationFactor = factor;
  }

  void InitializeNonLinIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {

      BaseType::InitializeNonLinIteration(rModelPart, A, Dx, b);

      ProcessInfo& ProcessInfo = rModelPart.GetProcessInfo();

      //if orthogonal subscales are computed
      if (ProcessInfo[OSS_SWITCH] == 1.0)
      {
        this->FullProjection(rModelPart);
        //this->LumpedProjection(rModelPart);
      }

    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                       TSystemMatrixType& rA,
                                       TSystemVectorType& rDx,
                                       TSystemVectorType& rb) override
    {
        BaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);
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
      KRATOS_INFO_IF("ResidualBasedSimpleSteadySchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0)<< "Computing OSS projections" << std::endl;

      ProcessInfo& ProcessInfo = rModelPart.GetProcessInfo();
      const int number_of_nodes = rModelPart.NumberOfNodes();
      const int number_of_elements = rModelPart.NumberOfElements();

      // NEW IMPLEMENTATION

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

      // OLD IMPLEMENTATION
      // unsigned int dimension = ProcessInfo[DOMAIN_SIZE];
      // VectorType MassProjectionRHS = ZeroVector(number_of_nodes);
      // VectorType MomentumProjectionRHS = ZeroVector(number_of_nodes*dimension);
      // if (mMassMatrixAlreadyComputed == false){
      //   mGlobalDivProjMassMatrix = ZeroMatrix(number_of_nodes,number_of_nodes);
      //   mGlobalAdvProjMassMatrix = ZeroMatrix(number_of_nodes*dimension,number_of_nodes*dimension);
      //   #pragma omp for schedule(guided, 512)
      //   for (int e = 0; e < number_of_elements; e++){
      //     ModelPart::ElementsContainerType::iterator it_elem = rModelPart.ElementsBegin() + e;
      //     GeometryType r_geometry = it_elem->GetGeometry();
      //     unsigned int NumNodes = r_geometry.PointsNumber();
      //     GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
      //     GeometryType::IntegrationPointsArrayType r_integrations_points = r_geometry.IntegrationPoints( integration_method );
      //     auto r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
      //     Vector detJ_vector(r_number_integration_points);
      //     r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
      //     Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
      //     for (unsigned int g = 0; g < r_number_integration_points; g++){
      //       double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
      //       for (unsigned int i = 0; i < NumNodes; ++i){
      //         for (unsigned int j= 0; j < NumNodes; ++j){
      //           for (unsigned int d = 0; d < dimension; ++d){
      //             mGlobalAdvProjMassMatrix(dimension*(r_geometry[i].Id()-1)+d,dimension*(r_geometry[j].Id()-1)+d) += Weight * NContainer(g,i) * NContainer(g,j);
      //           }
      //           mGlobalDivProjMassMatrix(r_geometry[i].Id()-1,r_geometry[j].Id()-1) += Weight * NContainer(g,i) * NContainer(g,j);
      //         }
      //       }
      //     }
      //   }
      //   mMassMatrixAlreadyComputed = true;
      // }

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
      KRATOS_INFO_IF("ResidualBasedSimpleSteadySchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0) << "Residual projection time: " << timer_projector << std::endl;
    }

void LumpedProjection(ModelPart& rModelPart)
    {

      KRATOS_INFO_IF("ResidualBasedSimpleSteadySchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0)
          << "Computing OSS projections" << std::endl;

      ProcessInfo& ProcessInfo = rModelPart.GetProcessInfo();

      const int number_of_nodes = rModelPart.NumberOfNodes();

      #pragma omp parallel for
      for (int i = 0; i < number_of_nodes; i++) {
        ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
        noalias(it_node->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
        it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
      }

      const int number_of_elements = rModelPart.NumberOfElements();
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
        if (it_node->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
          it_node->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
        const double area_inverse = 1.0 / it_node->FastGetSolutionStepValue(NODAL_AREA);
        it_node->FastGetSolutionStepValue(ADVPROJ) *= area_inverse;
        it_node->FastGetSolutionStepValue(DIVPROJ) *= area_inverse;
      }
    }
  ///@}

protected:

  ///@name Protected Operators
  ///@{

  ///@}

private:
  ///@name Member Variables
  ///@{

  double mVelocityRelaxationFactor;
  double mPressureRelaxationFactor;
  CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool = nullptr;
  Process::Pointer mpTurbulenceModel = nullptr;
  MatrixType mGlobalProjMassMatrix;
  bool mMassMatrixAlreadyComputed = false;
  typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

  ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME_DEM_COUPLED defined */
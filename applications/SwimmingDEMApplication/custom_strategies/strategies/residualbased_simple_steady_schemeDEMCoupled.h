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

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "linear_solvers/amgcl_solver.h"

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

  //typedef  SolverType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename TSparseSpace::MatrixType MatrixType;
  typedef typename TSparseSpace::VectorType VectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef Element::GeometryType GeometryType;

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

  void FullProjection(ModelPart& rModelPart)
    {

      KRATOS_INFO_IF("ResidualBasedSimpleSteadySchemeDEMCoupled", rModelPart.GetCommunicator().MyPID() == 0)<< "Computing OSS projections" << std::endl;

      ProcessInfo& ProcessInfo = rModelPart.GetProcessInfo();
      const int number_of_nodes = rModelPart.NumberOfNodes();
      const int number_of_elements = rModelPart.NumberOfElements();

      unsigned int dimension = ProcessInfo[DOMAIN_SIZE];
      VectorType MassProjectionRHS = ZeroVector(number_of_nodes);
      VectorType MomentumProjectionRHS = ZeroVector(number_of_nodes*dimension);
      if (mMassMatrixAlreadyComputed == false){
        mGlobalDivProjMassMatrix = ZeroMatrix(number_of_nodes,number_of_nodes);
        mGlobalAdvProjMassMatrix = ZeroMatrix(number_of_nodes*dimension,number_of_nodes*dimension);
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
  MatrixType mGlobalDivProjMassMatrix;
  MatrixType mGlobalAdvProjMassMatrix;
  bool mMassMatrixAlreadyComputed = false;
  typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

  ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME_DEM_COUPLED defined */
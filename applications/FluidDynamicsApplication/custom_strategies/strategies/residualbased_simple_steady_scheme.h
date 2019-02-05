//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//


#if !defined(KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME )
#define  KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME

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

namespace Kratos {

///@name Kratos Classes
///@{

template<class TSparseSpace, class TDenseSpace >
class ResidualBasedSimpleSteadyScheme : public Scheme<TSparseSpace, TDenseSpace> {
public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedSimpleSteadyScheme);

  typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

  typedef typename BaseType::DofsArrayType DofsArrayType;

  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

  typedef Element::GeometryType  GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  ResidualBasedSimpleSteadyScheme(double VelocityRelaxationFactor,
                                  double PressureRelaxationFactor,
                                  unsigned int DomainSize)
      : Scheme<TSparseSpace, TDenseSpace>(),
        mVelocityRelaxationFactor(VelocityRelaxationFactor),
        mPressureRelaxationFactor(PressureRelaxationFactor),
        mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0)
  {}

  ResidualBasedSimpleSteadyScheme(
      double VelocityRelaxationFactor,
      double PressureRelaxationFactor,
      unsigned int DomainSize,
      Process::Pointer pTurbulenceModel)
      : Scheme<TSparseSpace, TDenseSpace>(),
        mVelocityRelaxationFactor(VelocityRelaxationFactor),
        mPressureRelaxationFactor(PressureRelaxationFactor),
        mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0),
        mpTurbulenceModel(pTurbulenceModel)

  {}

  ~ResidualBasedSimpleSteadyScheme() override {}


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

  void Update(ModelPart& rModelPart,
                      DofsArrayType& rDofSet,
                      TSystemMatrixType& rA,
                      TSystemVectorType& rDx,
                      TSystemVectorType& rb) override
  {
    KRATOS_TRY;

    mRotationTool.RotateVelocities(rModelPart);

    mpDofUpdater->UpdateDofs(rDofSet,rDx);

    mRotationTool.RecoverVelocities(rModelPart);

    KRATOS_CATCH("");
  }

  void CalculateSystemContributions(
      Element::Pointer rCurrentElement,
      LocalSystemMatrixType& LHS_Contribution,
      LocalSystemVectorType& RHS_Contribution,
      Element::EquationIdVectorType& EquationId,
      ProcessInfo& CurrentProcessInfo) override
  {
    KRATOS_TRY;

    rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
    rCurrentElement->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

    Matrix SteadyLHS;
    rCurrentElement->CalculateLocalVelocityContribution(SteadyLHS, RHS_Contribution, CurrentProcessInfo);
    rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

    if (SteadyLHS.size1() != 0)
      noalias(LHS_Contribution) += SteadyLHS;

    AddRelaxation(rCurrentElement->GetGeometry(), LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

    // apply slip condition
    mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
    mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

    KRATOS_CATCH("");
  }

  void Condition_CalculateSystemContributions(
      Condition::Pointer rCurrentCondition,
      LocalSystemMatrixType& LHS_Contribution,
      LocalSystemVectorType& RHS_Contribution,
      Condition::EquationIdVectorType& EquationId,
      ProcessInfo& CurrentProcessInfo) override
  {
    KRATOS_TRY;

    rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
    rCurrentCondition->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

    Matrix SteadyLHS;
    rCurrentCondition->CalculateLocalVelocityContribution(SteadyLHS, RHS_Contribution, CurrentProcessInfo);
    rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

    if (SteadyLHS.size1() != 0)
      noalias(LHS_Contribution) += SteadyLHS;

    AddRelaxation(rCurrentCondition->GetGeometry(), LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

    // apply slip condition
    mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
    mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

    KRATOS_CATCH("");
  }

  void Calculate_RHS_Contribution(
      Element::Pointer rCurrentElement,
      LocalSystemVectorType& rRHS_Contribution,
      Element::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY;

    Matrix LHS_Contribution;
    CalculateSystemContributions(rCurrentElement,LHS_Contribution,
                                 rRHS_Contribution,rEquationId,rCurrentProcessInfo);

    KRATOS_CATCH("");
  }

  void Condition_Calculate_RHS_Contribution(
      Condition::Pointer rCurrentCondition,
      LocalSystemVectorType& rRHS_Contribution,
      Element::EquationIdVectorType& rEquationId,
      ProcessInfo& rCurrentProcessInfo) override
  {
    KRATOS_TRY;

    Matrix LHS_Contribution;
    Condition_CalculateSystemContributions(rCurrentCondition,LHS_Contribution,
                                           rRHS_Contribution,rEquationId,
                                           rCurrentProcessInfo);

    KRATOS_CATCH("");
  }

  void InitializeNonLinIteration(ModelPart& rModelPart,
                                         TSystemMatrixType& rA,
                                         TSystemVectorType& rDx,
                                         TSystemVectorType& rb) override
  {
    KRATOS_TRY;

    for (typename ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin();
         itNode != rModelPart.NodesEnd(); itNode++)
      itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;

    double output;
    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
    const int number_of_elements = rModelPart.NumberOfElements();
    #pragma omp parallel for private(output)
    for (int i = 0; i < number_of_elements; i++) {
      ModelPart::ElementsContainerType::iterator it_elem = rModelPart.ElementsBegin() + i;
      it_elem->Calculate(NODAL_AREA, output, CurrentProcessInfo);
    }

    rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

    if (mpTurbulenceModel != 0) // If not null
      mpTurbulenceModel->Execute();

    KRATOS_CATCH("");
  }

  void FinalizeNonLinIteration(ModelPart &rModelPart,
                                       TSystemMatrixType &rA,
                                       TSystemVectorType &rDx,
                                       TSystemVectorType &rb) override
  {
    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

    //if orthogonal subscales are computed
    if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {

      KRATOS_INFO_IF("ResidualBasedSimpleSteadyScheme", rModelPart.GetCommunicator().MyPID() == 0)
          << "Computing OSS projections" << std::endl;

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
        it_elem->Calculate(ADVPROJ,output,CurrentProcessInfo);
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
  }

  void FinalizeSolutionStep(ModelPart &rModelPart,
                            TSystemMatrixType &rA,
                            TSystemVectorType &rDx,
                            TSystemVectorType &rb) override
  {
    LocalSystemVectorType RHS_Contribution;
    LocalSystemMatrixType LHS_Contribution;
    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

    for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin();
         itNode != rModelPart.NodesEnd(); ++itNode)
    {
      itNode->FastGetSolutionStepValue(REACTION_X,0) = 0.0;
      itNode->FastGetSolutionStepValue(REACTION_Y,0) = 0.0;
      itNode->FastGetSolutionStepValue(REACTION_Z,0) = 0.0;
    }

    for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin();
         itElem != rModelPart.Elements().ptr_end(); ++itElem)
    {
      (*itElem)->InitializeNonLinearIteration(rCurrentProcessInfo);
      (*itElem)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
      Matrix SteadyLHS;
      (*itElem)->CalculateLocalVelocityContribution(SteadyLHS,RHS_Contribution,rCurrentProcessInfo);

      GeometryType& rGeom = (*itElem)->GetGeometry();
      unsigned int NumNodes = rGeom.PointsNumber();
      unsigned int Dimension = rGeom.WorkingSpaceDimension();
      unsigned int index = 0;

      for (unsigned int i = 0; i < NumNodes; i++)
      {
        rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
        rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
        if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
        index++; // skip pressure dof
      }
    }

    rModelPart.GetCommunicator().AssembleCurrentData(REACTION);
    Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
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
    this->CalculateLumpedMassMatrix(rGeometry,Mass);

    unsigned int DofIndex = 0;
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      const array_1d<double, 3>& rVel = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY,0);
      const double Area = rGeometry[iNode].FastGetSolutionStepValue(NODAL_AREA,0);
      double VelNorm = 0.0;
      for (unsigned int d = 0; d < Dimension; ++d)
        VelNorm += rVel[d] * rVel[d];
      VelNorm = sqrt(VelNorm);
      double LocalDt;
      if (VelNorm != 0.0)
        LocalDt = pow(Area, 1.0 / double(Dimension)) / VelNorm;
      else
        LocalDt = 1.0;

      for (unsigned int i = 0; i < Dimension; i++)
      {
        Mass(DofIndex,DofIndex) *= 1.0 / (mVelocityRelaxationFactor * LocalDt);
        DofIndex++;
      }
      DofIndex++; // pressure dof
    }
    noalias(LHS_Contribution) += Mass;

    // pressure relaxation
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      unsigned int BlockIndex = iNode * (Dimension + 1);
      LHS_Contribution(BlockIndex+Dimension,BlockIndex+Dimension) *= 1.0 / mPressureRelaxationFactor;
    }
  }

  void CalculateLumpedMassMatrix(
    const GeometryType& rGeometry,
    LocalSystemMatrixType& rLumpedMass) const
  {
    const unsigned int dimension = rGeometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = rGeometry.PointsNumber();

    const unsigned int nodal_block_size = dimension + 1;
    const unsigned int local_size = nodal_block_size * number_of_nodes;

    if (rLumpedMass.size1() != local_size) {
      rLumpedMass.resize(local_size,local_size,false);
    }

    noalias(rLumpedMass) = ZeroMatrix(local_size,local_size);

    const double size_fraction = rGeometry.DomainSize() / number_of_nodes;

    for (unsigned int i = 0; i < number_of_nodes; i++){
      const unsigned int node_block = i*nodal_block_size;
      const double lumped_mass = size_fraction * rGeometry[i].FastGetSolutionStepValue(DENSITY);
      for (unsigned int d = 0; d < dimension; d++) {
        rLumpedMass(node_block+d,node_block+d) = lumped_mass;
      }
    }
  }

  ///@}

private:
  ///@name Member Variables
  ///@{

  double mVelocityRelaxationFactor;
  double mPressureRelaxationFactor;
  CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;
  Process::Pointer mpTurbulenceModel;
  typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

  ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_SIMPLE_STEADY_SCHEME defined */

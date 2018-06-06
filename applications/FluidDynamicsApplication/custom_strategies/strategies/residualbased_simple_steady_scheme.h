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

    std::cout << "VelocityRelaxationFactor = " << mVelocityRelaxationFactor << std::endl;
    std::cout << "PressureRelaxationFactor = " << mPressureRelaxationFactor << std::endl;
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
    Matrix Mass;
    rCurrentElement->CalculateMassMatrix(Mass, CurrentProcessInfo);
    Matrix SteadyLHS;
    rCurrentElement->CalculateLocalVelocityContribution(SteadyLHS, RHS_Contribution, CurrentProcessInfo);
    rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

    if (SteadyLHS.size1() != 0)
      noalias(LHS_Contribution) += SteadyLHS;

    AddRelaxation(rCurrentElement, LHS_Contribution, RHS_Contribution, Mass, CurrentProcessInfo);

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
    Matrix Mass;
    rCurrentCondition->CalculateMassMatrix(Mass, CurrentProcessInfo);
    Matrix SteadyLHS;
    rCurrentCondition->CalculateLocalVelocityContribution(SteadyLHS, RHS_Contribution, CurrentProcessInfo);
    rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

    if (SteadyLHS.size1() != 0)
      noalias(LHS_Contribution) += SteadyLHS;

    AddRelaxation(rCurrentCondition, LHS_Contribution, RHS_Contribution, Mass, CurrentProcessInfo);

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
    for (typename ModelPart::ElementsContainerType::iterator itElem = rModelPart.ElementsBegin();
         itElem != rModelPart.ElementsEnd(); itElem++)
      itElem->Calculate(NODAL_AREA, output, CurrentProcessInfo);

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

      if (rModelPart.GetCommunicator().MyPID() == 0)
        std::cout << "Computing OSS projections" << std::endl;

      for (typename ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin();
           itNode != rModelPart.NodesEnd(); itNode++)
      {
        noalias(itNode->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);
        itNode->FastGetSolutionStepValue(DIVPROJ) = 0.0;
        itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
      }

      array_1d<double, 3 > output;
      for (typename ModelPart::ElementsContainerType::iterator itElem = rModelPart.ElementsBegin();
           itElem != rModelPart.ElementsEnd(); itElem++)
      {
        itElem->Calculate(ADVPROJ, output, CurrentProcessInfo);
      }
      rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
      rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
      rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

      for (typename ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin();
           itNode != rModelPart.NodesEnd(); itNode++)
      {
        if (itNode->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
          itNode->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
        const double Area = itNode->FastGetSolutionStepValue(NODAL_AREA);
        itNode->FastGetSolutionStepValue(ADVPROJ) /= Area;
        itNode->FastGetSolutionStepValue(DIVPROJ) /= Area;
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
  void AddRelaxation(Element::Pointer rCurrentElement,
                     LocalSystemMatrixType& LHS_Contribution,
                     LocalSystemVectorType& RHS_Contribution,
                     LocalSystemMatrixType& Mass,
                     ProcessInfo& CurrentProcessInfo)
  {
    if (Mass.size1() == 0)
      return;

    GeometryType& rGeom = rCurrentElement->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int Dimension = rGeom.WorkingSpaceDimension();
    const unsigned int VelocityBlockSize = NumNodes * Dimension;
    unsigned int DofIndex = 0;
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      array_1d<double, 3>& rVel = rGeom[iNode].FastGetSolutionStepValue(VELOCITY,0);
      double Area = rGeom[iNode].FastGetSolutionStepValue(NODAL_AREA,0);
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
        double LumpedMass = 0.0;
        for (unsigned int j = 0; j < VelocityBlockSize; j++)
        {
          LumpedMass += Mass(DofIndex,j);
          Mass(DofIndex,j) = 0.0;
        }
        // the relaxation factor defines the local cfl number
        Mass(DofIndex,DofIndex) = LumpedMass / (mVelocityRelaxationFactor * LocalDt);
        DofIndex++;
      }
      DofIndex++; // pressure dof
    }
    noalias(LHS_Contribution) += Mass;

    // pressure relaxation
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      unsigned int BlockIndex = iNode * (Dimension + 1);
      LHS_Contribution(BlockIndex+Dimension,BlockIndex+Dimension) /= mPressureRelaxationFactor;
    }
  }

  void AddRelaxation(Condition::Pointer rCurrentCondition,
                     LocalSystemMatrixType& LHS_Contribution,
                     LocalSystemVectorType& RHS_Contribution,
                     LocalSystemMatrixType& Mass,
                     ProcessInfo& CurrentProcessInfo)
  {
    if (Mass.size1() == 0)
      return;

    GeometryType& rGeom = rCurrentCondition->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int Dimension = rGeom.WorkingSpaceDimension();
    const unsigned int VelocityBlockSize = NumNodes * Dimension;
    unsigned int DofIndex = 0;
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      array_1d<double, 3>& rVel = rGeom[iNode].FastGetSolutionStepValue(VELOCITY,0);
      double Area = rGeom[iNode].FastGetSolutionStepValue(NODAL_AREA,0);
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
        double LumpedMass = 0.0;
        for (unsigned int j = 0; j < VelocityBlockSize; j++)
        {
          LumpedMass += Mass(DofIndex,j);
          Mass(DofIndex,j) = 0.0;
        }
        // the relaxation factor defines the local cfl number
        Mass(DofIndex,DofIndex) = LumpedMass / (mVelocityRelaxationFactor * LocalDt);
        DofIndex++;
      }
      DofIndex++; // pressure dof
    }
    noalias(LHS_Contribution) += Mass;

    // pressure relaxation
    for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
    {
      unsigned int BlockIndex = iNode * (Dimension + 1);
      LHS_Contribution(BlockIndex+Dimension,BlockIndex+Dimension) /= mPressureRelaxationFactor;
    }
//Vector OldDofValues;
    //rCurrentCondition->GetFirstDerivativesVector(OldDofValues, 0);
    //noalias(RHS_Contribution) -= prod(Mass, OldDofValues);
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

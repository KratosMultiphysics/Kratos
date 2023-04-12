//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_V_P_STRATEGY_H
#define KRATOS_V_P_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"

#include "custom_utilities/solver_settings.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include <stdio.h>
#include <cmath>

namespace Kratos
{

  ///@addtogroup PFEMFluidDynamicsApplication
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

  template <class TSparseSpace,
            class TDenseSpace,
            class TLinearSolver>
  class VPStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace>
  {
  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(VPStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> BaseType;

    typedef TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    VPStrategy(ModelPart &rModelPart,
               SolverSettingsType &rSolverConfig) : BaseType(rModelPart)
    {
      std::cout << "VPStrategy INITIALIZE STRATEGY" << std::endl;
      InitializeStrategy(rSolverConfig);
    }

    VPStrategy(ModelPart &rModelPart,
               typename TLinearSolver::Pointer pVelocityLinearSolver,
               typename TLinearSolver::Pointer pPressureLinearSolver,
               bool ReformDofSet = true,
               unsigned int DomainSize = 2) : BaseType(rModelPart)
    {
      KRATOS_TRY;
      KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~VPStrategy() {}

    virtual int Check() override
    {
      return false;
    }

    virtual bool SolveSolutionStep() override
    {
      return false;
    }

    virtual void FinalizeSolutionStep() override {}

    virtual void InitializeSolutionStep() override {}

    void UpdateTopology(ModelPart &rModelPart, unsigned int echoLevel)
    {
      KRATOS_TRY;

      this->CalculateDisplacementsAndPorosity();
      BaseType::MoveMesh();
      KRATOS_CATCH("");
    }

    void SetBlockedAndIsolatedFlags()
    {
      KRATOS_TRY;

      ModelPart &rModelPart = BaseType::GetModelPart();
      const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

#pragma omp parallel
      {
        ModelPart::ElementIterator ElemBegin;
        ModelPart::ElementIterator ElemEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
        for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
        {
          unsigned int numNodes = itElem->GetGeometry().size();
          std::vector<array_1d<double, 3>> nodesCoordinates;
          nodesCoordinates.resize(numNodes);

          (itElem)->Set(BLOCKED, false);
          (itElem)->Set(ISOLATED, false);

          unsigned int freeSurfaceNodes = 0;
          unsigned int freeSurfaceRigidNodes = 0;
          unsigned int rigidNodes = 0;
          unsigned int isolatedNodes = 0;
          for (unsigned int i = 0; i < numNodes; i++)
          {
            if (itElem->GetGeometry()[i].Is(FREE_SURFACE))
            {
              freeSurfaceNodes++;
              if (itElem->GetGeometry()[i].Is(RIGID))
              {
                freeSurfaceRigidNodes++;
              }
            }
            else if (itElem->GetGeometry()[i].Is(RIGID))
            {
              rigidNodes++;
            }
            nodesCoordinates[i] = itElem->GetGeometry()[i].Coordinates();
            ElementWeakPtrVectorType &neighb_elems = itElem->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
            if (neighb_elems.size() == 1)
            {
              isolatedNodes++;
            }
          }

          if (dimension == 3)
          {
            double a1 = 0; // slope x for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            double b1 = 0; // slope y for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            double c1 = 0; // slope z for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            a1 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[2][2] - nodesCoordinates[0][2]) - (nodesCoordinates[2][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
            b1 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[2][0] - nodesCoordinates[0][0]) - (nodesCoordinates[2][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
            c1 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[2][1] - nodesCoordinates[0][1]) - (nodesCoordinates[2][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
            double a2 = 0; // slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            double b2 = 0; // slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            double c2 = 0; // slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            a2 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[3][2] - nodesCoordinates[0][2]) - (nodesCoordinates[3][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
            b2 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[3][0] - nodesCoordinates[0][0]) - (nodesCoordinates[3][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
            c2 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[3][1] - nodesCoordinates[0][1]) - (nodesCoordinates[3][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
            double a3 = 0; // slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            double b3 = 0; // slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            double c3 = 0; // slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            a3 = (nodesCoordinates[1][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[1][2] - nodesCoordinates[2][2]);
            b3 = (nodesCoordinates[1][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[1][0] - nodesCoordinates[2][0]);
            c3 = (nodesCoordinates[1][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[1][1] - nodesCoordinates[2][1]);
            double a4 = 0; // slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            double b4 = 0; // slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            double c4 = 0; // slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            a4 = (nodesCoordinates[0][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[0][2] - nodesCoordinates[2][2]);
            b4 = (nodesCoordinates[0][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[0][0] - nodesCoordinates[2][0]);
            c4 = (nodesCoordinates[0][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[0][1] - nodesCoordinates[2][1]);

            double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
            double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
            double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
            double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
            double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
            double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

            if ((fabs(cosAngle12) > 0.99 || fabs(cosAngle13) > 0.99 || fabs(cosAngle14) > 0.99 || fabs(cosAngle23) > 0.99 || fabs(cosAngle24) > 0.99 || fabs(cosAngle34) > 0.99) && (freeSurfaceNodes == numNodes) && isolatedNodes > 1)
            {
              (itElem)->Set(BLOCKED, true);
              // std::cout << "in the strategy BLOCKED ELEMENT: " << (itElem)->Id() << std::endl;
            }
            else if ((fabs(cosAngle12) > 0.995 || fabs(cosAngle13) > 0.995 || fabs(cosAngle14) > 0.995 || fabs(cosAngle23) > 0.995 || fabs(cosAngle24) > 0.995 || fabs(cosAngle34) > 0.995) && (freeSurfaceNodes == numNodes) && isolatedNodes == 1)
            {
              (itElem)->Set(BLOCKED, true);
              // std::cout << "in the strategy BLOCKED ELEMENT: " << (itElem)->Id() << std::endl;
            }
            else if ((fabs(cosAngle12) > 0.999 || fabs(cosAngle13) > 0.999 || fabs(cosAngle14) > 0.999 || fabs(cosAngle23) > 0.999 || fabs(cosAngle24) > 0.999 || fabs(cosAngle34) > 0.999) && (freeSurfaceNodes == numNodes))
            {
              (itElem)->Set(BLOCKED, true);
              // std::cout << "in the strategy BLOCKED ELEMENT: " << (itElem)->Id() << std::endl;
            }
          }

          if (freeSurfaceNodes == numNodes && rigidNodes == 0 && isolatedNodes >= (numNodes - 1))
          {
            (itElem)->Set(ISOLATED, true);
            (itElem)->Set(BLOCKED, false);
          }
        }
      }
      KRATOS_CATCH("");
    }

    void CalculatePressureVelocity()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      unsigned int timeStep = rCurrentProcessInfo[STEP];

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
           i != rModelPart.NodesEnd(); ++i)
      {
        if (timeStep == 1)
        {
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0;
        }
        else
        {
          double &CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE, 0);
          double &PreviousPressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
          double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
          CurrentPressureVelocity = (CurrentPressure - PreviousPressure) / timeInterval;
        }
      }
    }

    void CalculatePressureAcceleration()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      unsigned int timeStep = rCurrentProcessInfo[STEP];

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
           i != rModelPart.NodesEnd(); ++i)
      {
        if (timeStep == 1)
        {
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0;
        }
        else
        {
          double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
          double &PreviousPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1);
          double &CurrentPressureAcceleration = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);
          CurrentPressureAcceleration = (CurrentPressureVelocity - PreviousPressureVelocity) / timeInterval;
        }
      }
    }

    virtual void CalculateTemporalVariables()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
           i != rModelPart.NodesEnd(); ++i)
      {

        array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
        array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

        array_1d<double, 3> &CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
        array_1d<double, 3> &PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

        /* if((i)->IsNot(ISOLATED) || (i)->Is(SOLID)){ */
        if ((i)->IsNot(ISOLATED) && ((i)->IsNot(RIGID) || (i)->Is(SOLID)))
        {
          UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity);
        }
        else if ((i)->Is(RIGID))
        {
          array_1d<double, 3> Zeros(3, 0.0);
          (i)->FastGetSolutionStepValue(ACCELERATION, 0) = Zeros;
          (i)->FastGetSolutionStepValue(ACCELERATION, 1) = Zeros;
        }
        else
        {
          (i)->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0.0;
          if ((i)->SolutionStepsDataHas(VOLUME_ACCELERATION))
          {
            array_1d<double, 3> &VolumeAcceleration = (i)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
            (i)->FastGetSolutionStepValue(ACCELERATION, 0) = VolumeAcceleration;
            (i)->FastGetSolutionStepValue(VELOCITY, 0) += VolumeAcceleration * rCurrentProcessInfo[DELTA_TIME];
          }
        }

        const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
        unsigned int timeStep = rCurrentProcessInfo[STEP];
        if (timeStep == 1)
        {
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0;
        }
        else
        {
          double &CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE, 0);
          double &PreviousPressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
          double &CurrentPressureVelocity = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
          double &CurrentPressureAcceleration = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);

          CurrentPressureAcceleration = CurrentPressureVelocity / timeInterval;

          CurrentPressureVelocity = (CurrentPressure - PreviousPressure) / timeInterval;

          CurrentPressureAcceleration += -CurrentPressureVelocity / timeInterval;
        }
      }
    }

    void CalculateAccelerations()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
           i != rModelPart.NodesEnd(); ++i)
      {

        array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
        array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

        array_1d<double, 3> &CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
        array_1d<double, 3> &PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

        /* if((i)->IsNot(ISOLATED) || (i)->Is(SOLID)){ */
        if ((i)->IsNot(ISOLATED) && ((i)->IsNot(RIGID) || (i)->Is(SOLID)))
        {
          UpdateAccelerations(CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity);
        }
        else if ((i)->Is(RIGID))
        {
          array_1d<double, 3> Zeros(3, 0.0);
          (i)->FastGetSolutionStepValue(ACCELERATION, 0) = Zeros;
          (i)->FastGetSolutionStepValue(ACCELERATION, 1) = Zeros;
        }
        else
        {
          (i)->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0) = 0.0;
          (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1) = 0.0;
          if ((i)->SolutionStepsDataHas(VOLUME_ACCELERATION))
          {
            array_1d<double, 3> &VolumeAcceleration = (i)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
            (i)->FastGetSolutionStepValue(ACCELERATION, 0) = VolumeAcceleration;
            (i)->FastGetSolutionStepValue(VELOCITY, 0) += VolumeAcceleration * rCurrentProcessInfo[DELTA_TIME];
          }
        }
      }
    }

    inline void UpdateAccelerations(array_1d<double, 3> &CurrentAcceleration,
                                    const array_1d<double, 3> &CurrentVelocity,
                                    array_1d<double, 3> &PreviousAcceleration,
                                    const array_1d<double, 3> &PreviousVelocity)
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double Dt = rCurrentProcessInfo[DELTA_TIME];
      noalias(CurrentAcceleration) = 2.0 * (CurrentVelocity - PreviousVelocity) / Dt - PreviousAcceleration;
    }

    virtual void CalculateDisplacementsAndPorosity()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
           i != rModelPart.NodesEnd(); ++i)
      {
        if ((i)->IsNot(PFEMFlags::EULERIAN_INLET))

        {
          array_1d<double, 3> &CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 0);
          array_1d<double, 3> &PreviousVelocity = (i)->FastGetSolutionStepValue(VELOCITY, 1);

          array_1d<double, 3> &CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
          array_1d<double, 3> &PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);

          if(!i->IsFixed(DISPLACEMENT_X))
            CurrentDisplacement[0] = 0.5 * TimeStep * (CurrentVelocity[0] + PreviousVelocity[0]) + PreviousDisplacement[0];

          if(!i->IsFixed(DISPLACEMENT_Y))
            CurrentDisplacement[1] = 0.5 * TimeStep * (CurrentVelocity[1] + PreviousVelocity[1]) + PreviousDisplacement[1];

          if(!i->IsFixed(DISPLACEMENT_Z))
            CurrentDisplacement[2] = 0.5 * TimeStep * (CurrentVelocity[2] + PreviousVelocity[2]) + PreviousDisplacement[2];

          // currentFluidFractionRate = (currentFluidFraction - previousFluidFraction)/TimeStep;
        }
      }
    }

    virtual void UpdateStressStrain() {}

    virtual void Clear() override {}

    ///@}
    ///@name Access
    ///@{

    virtual void SetEchoLevel(int Level) override
    {
      BaseType::SetEchoLevel(Level);
    }

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
      buffer << "VPStrategy";
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "VPStrategy";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected Life Cycle
    ///@{

    ///@}
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

    /// Calculate the coefficients for time iteration.
    /**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME variables.
     */

    virtual bool SolveMomentumIteration(unsigned int it, unsigned int maxIt, bool &fixedTimeStep, double &velocityNorm)
    {
      return false;
    }

    virtual bool SolveContinuityIteration(unsigned int it, unsigned int maxIt, double &NormP)
    {
      return false;
    }

    void ComputeErrorL2Norm(double tensilStressSign) // tensilStressSign = 1.0 for FIC, tensilStressSign = -1.0 for FS
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double currentTime = rCurrentProcessInfo[TIME];
      const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

      long double sumErrorL2Velocity = 0;
      long double sumErrorL2VelocityX = 0;
      long double sumErrorL2VelocityY = 0;
      long double sumErrorL2Pressure = 0;
      long double sumErrorL2TauXX = 0;
      long double sumErrorL2TauYY = 0;
      long double sumErrorL2TauXY = 0;

#pragma omp parallel
      {
        ModelPart::ElementIterator ElemBegin;
        ModelPart::ElementIterator ElemEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
        for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
        {

          Element::GeometryType &geometry = itElem->GetGeometry();
          long double nodalArea = 0;

          if (dimension == 2)
          {
            nodalArea = geometry.Area() / 3.0;
          }
          else if (dimension == 3)
          {
            nodalArea = geometry.Volume() * 0.25;
          }

          long double bariPosX = 0;
          long double bariPosY = 0;

          long double eleErrorL2Velocity = 0;
          long double eleErrorL2VelocityX = 0;
          long double eleErrorL2VelocityY = 0;
          long double eleErrorL2Pressure = 0;

          // ShapeFunctionDerivativesArrayType DN_DX;
          Matrix NContainer;
          NContainer = geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

          const Vector &N = row(NContainer, 0);
          const unsigned int NumNodes = geometry.size();

          double elementalPressure = N[0] * geometry(0)->FastGetSolutionStepValue(PRESSURE);
          double elementalVelocityX = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_X);
          double elementalVelocityY = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_Y);
          ;

          for (unsigned int i = 1; i < NumNodes; i++)
          {
            elementalPressure += N[i] * geometry(i)->FastGetSolutionStepValue(PRESSURE);
            elementalVelocityX += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_X);
            elementalVelocityY += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_Y);
          }

          for (unsigned int i = 0; i < geometry.size(); i++)
          {

            const long double nodalPosX = geometry(i)->X();
            const long double nodalPosY = geometry(i)->Y();

            bariPosX += nodalPosX / 3.0;
            bariPosY += nodalPosY / 3.0;
          }

          const long double posX = bariPosX;
          const long double posY = bariPosY;
          long double expectedVelocityX = pow(posX, 2) * (1.0 - posX) * (1.0 - posX) * (2.0 * posY - 6.0 * pow(posY, 2) + 4.0 * pow(posY, 3));
          long double expectedVelocityY = -pow(posY, 2) * (1.0 - posY) * (1.0 - posY) * (2.0 * posX - 6.0 * pow(posX, 2) + 4.0 * pow(posX, 3));
          long double expectedPressure = -tensilStressSign * posX * (1.0 - posX);

          eleErrorL2VelocityX = elementalVelocityX - expectedVelocityX;
          eleErrorL2VelocityY = elementalVelocityY - expectedVelocityY;
          eleErrorL2Pressure = elementalPressure - expectedPressure;

          sumErrorL2VelocityX += pow(eleErrorL2VelocityX, 2) * geometry.Area();
          sumErrorL2VelocityY += pow(eleErrorL2VelocityY, 2) * geometry.Area();
          sumErrorL2Pressure += pow(eleErrorL2Pressure, 2) * geometry.Area();

          const long double tauXX = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XX);
          const long double tauYY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_YY);
          const long double tauXY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XY);

          long double expectedTauXX = 2.0 * (-4.0 * (1.0 - bariPosX) * bariPosX * (-1.0 + 2.0 * bariPosX) * bariPosY * (1.0 - 3.0 * bariPosY + 2.0 * pow(bariPosY, 2)));
          long double expectedTauYY = 2.0 * (4.0 * bariPosX * (1.0 - 3.0 * bariPosX + 2.0 * pow(bariPosX, 2)) * (1.0 - bariPosY) * bariPosY * (-1.0 + 2.0 * bariPosY));
          long double expectedTauXY = (2.0 * (1.0 - 6.0 * bariPosY + 6.0 * pow(bariPosY, 2)) * (1.0 - bariPosX) * (1.0 - bariPosX) * pow(bariPosX, 2) - 2.0 * (1.0 - 6.0 * bariPosX + 6.0 * pow(bariPosX, 2)) * (1.0 - bariPosY) * (1 - bariPosY) * pow(bariPosY, 2));

          long double nodalErrorTauXX = tauXX - expectedTauXX;
          long double nodalErrorTauYY = tauYY - expectedTauYY;
          long double nodalErrorTauXY = tauXY - expectedTauXY;

          sumErrorL2TauXX += pow(nodalErrorTauXX, 2) * geometry.Area();
          sumErrorL2TauYY += pow(nodalErrorTauYY, 2) * geometry.Area();
          sumErrorL2TauXY += pow(nodalErrorTauXY, 2) * geometry.Area();
        }
      }

      long double errorL2Velocity = sqrt(sumErrorL2Velocity);
      long double errorL2VelocityX = sqrt(sumErrorL2VelocityX);
      long double errorL2VelocityY = sqrt(sumErrorL2VelocityY);
      long double errorL2Pressure = sqrt(sumErrorL2Pressure);
      long double errorL2TauXX = sqrt(sumErrorL2TauXX);
      long double errorL2TauYY = sqrt(sumErrorL2TauYY);
      long double errorL2TauXY = sqrt(sumErrorL2TauXY);

      std::ofstream myfileVelocity;
      myfileVelocity.open("errorL2VelocityFile.txt", std::ios::app);
      myfileVelocity << currentTime << "\t" << errorL2Velocity << "\n";
      myfileVelocity.close();

      std::ofstream myfileVelocityX;
      myfileVelocityX.open("errorL2VelocityXFile.txt", std::ios::app);
      myfileVelocityX << currentTime << "\t" << errorL2VelocityX << "\n";
      myfileVelocityX.close();

      std::ofstream myfileVelocityY;
      myfileVelocityY.open("errorL2VelocityYFile.txt", std::ios::app);
      myfileVelocityY << currentTime << "\t" << errorL2VelocityY << "\n";
      myfileVelocityY.close();

      std::ofstream myfilePressure;
      myfilePressure.open("errorL2PressureFile.txt", std::ios::app);
      myfilePressure << currentTime << "\t" << errorL2Pressure << "\n";
      myfilePressure.close();

      std::ofstream myfileTauXX;
      myfileTauXX.open("errorL2TauXXFile.txt", std::ios::app);
      myfileTauXX << currentTime << "\t" << errorL2TauXX << "\n";
      myfileTauXX.close();

      std::ofstream myfileTauYY;
      myfileTauYY.open("errorL2TauYYFile.txt", std::ios::app);
      myfileTauYY << currentTime << "\t" << errorL2TauYY << "\n";
      myfileTauYY.close();

      std::ofstream myfileTauXY;
      myfileTauXY.open("errorL2TauXYFile.txt", std::ios::app);
      myfileTauXY << currentTime << "\t" << errorL2TauXY << "\n";
      myfileTauXY.close();
    }

    void ComputeErrorL2NormCasePoiseuille()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double currentTime = rCurrentProcessInfo[TIME];
      const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

      double sumErrorL2VelocityTheta = 0;
      double sumErrorL2TauTheta = 0;

      double r_in = 0.2;
      double R_out = 0.5;
      double kappa = r_in / R_out;
      double omega = 0.5;
      double viscosity = 100.0;

#pragma omp parallel
      {
        ModelPart::ElementIterator ElemBegin;
        ModelPart::ElementIterator ElemEnd;
        OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElemBegin, ElemEnd);
        for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
        {

          Element::GeometryType &geometry = itElem->GetGeometry();
          long double nodalArea = 0;

          if (dimension == 2)
          {
            nodalArea = geometry.Area() / 3.0;
          }
          else if (dimension == 3)
          {
            nodalArea = geometry.Volume() * 0.25;
          }

          long double bariPosX = 0;
          long double bariPosY = 0;

          long double eleErrorL2Velocity = 0;
          long double eleErrorL2VelocityX = 0;
          long double eleErrorL2VelocityY = 0;
          long double eleErrorL2Pressure = 0;

          Matrix NContainer;
          NContainer = geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

          const Vector &N = row(NContainer, 0);
          //  itElem->EvaluateInPoint(elementalPressure,PRESSURE,N);
          const unsigned int NumNodes = geometry.size();

          double elementalPressure = N[0] * geometry(0)->FastGetSolutionStepValue(PRESSURE);
          double elementalVelocityX = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_X);
          double elementalVelocityY = N[0] * geometry(0)->FastGetSolutionStepValue(VELOCITY_Y);
          ;

          for (unsigned int i = 1; i < NumNodes; i++)
          {
            elementalPressure += N[i] * geometry(i)->FastGetSolutionStepValue(PRESSURE);
            elementalVelocityX += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_X);
            elementalVelocityY += N[i] * geometry(i)->FastGetSolutionStepValue(VELOCITY_Y);
          }

          for (unsigned int i = 0; i < geometry.size(); i++)
          {

            // index = i*dimension;
            const long double nodalPosX = geometry(i)->X();
            const long double nodalPosY = geometry(i)->Y();

            bariPosX += nodalPosX / 3.0;
            bariPosY += nodalPosY / 3.0;
          }

          const long double posX = bariPosX;
          const long double posY = bariPosY;
          const double rPos = sqrt(pow(posX, 2) + pow(posY, 2));
          const double cosalfa = posX / rPos;
          const double sinalfa = posY / rPos;
          const double sin2alfa = 2.0 * cosalfa * sinalfa;
          const double cos2alfa = 1.0 - 2.0 * pow(sinalfa, 2);

          double expectedVelocityTheta = pow(kappa, 2) * omega * R_out / (1.0 - pow(kappa, 2)) * (R_out / rPos - rPos / R_out);
          double computedVelocityTheta = sqrt(pow(elementalVelocityX, 2) + pow(elementalVelocityY, 2));
          double nodalErrorVelocityTheta = computedVelocityTheta - expectedVelocityTheta;

          const long double tauXX = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XX);
          const long double tauYY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_YY);
          const long double tauXY = 0; // itElem->GetValue(ELEMENTAL_DEVIATORIC_STRESS_XY);

          double expectedTauTheta = (2.0 * viscosity * pow(kappa, 2) * omega * pow(R_out, 2)) / (1.0 - pow(kappa, 2)) / pow(rPos, 2);
          double computedTauTheta = (tauXX - tauYY) * sin2alfa / 2.0 - tauXY * cos2alfa;
          double nodalErrorTauTheta = computedTauTheta - expectedTauTheta;

          sumErrorL2VelocityTheta += pow(nodalErrorVelocityTheta, 2) * geometry.Area();
          sumErrorL2TauTheta += pow(nodalErrorTauTheta, 2) * geometry.Area();
        }
      }

      double errorL2VelocityTheta = sqrt(sumErrorL2VelocityTheta);
      double errorL2TauTheta = sqrt(sumErrorL2TauTheta);

      std::ofstream myfileVelocity;
      myfileVelocity.open("errorL2Poiseuille.txt", std::ios::app);
      myfileVelocity << currentTime << "\t" << errorL2VelocityTheta << "\t" << errorL2TauTheta << "\n";
      myfileVelocity.close();
    }

    double ComputeVelocityNorm()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const int n_nodes = rModelPart.NumberOfNodes();

      double NormV = 0.00;

#pragma omp parallel for reduction(+ \
                                   : NormV)
      for (int i_node = 0; i_node < n_nodes; ++i_node)
      {
        const auto it_node = rModelPart.NodesBegin() + i_node;
        const auto &r_vel = it_node->FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < 3; ++d)
        {
          NormV += r_vel[d] * r_vel[d];
        }
      }
      NormV = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormV);
      NormV = sqrt(NormV);

      const double zero_tol = 1.0e-12;
      if (NormV < zero_tol)
        NormV = 1.00;

      return NormV;
    }

    double ComputePressureNorm()
    {
      ModelPart &rModelPart = BaseType::GetModelPart();
      const int n_nodes = rModelPart.NumberOfNodes();

      double NormP = 0.00;

#pragma omp parallel for reduction(+ \
                                   : NormP)
      for (int i_node = 0; i_node < n_nodes; ++i_node)
      {
        const auto it_node = rModelPart.NodesBegin() + i_node;
        const double Pr = it_node->FastGetSolutionStepValue(PRESSURE);
        NormP += Pr * Pr;
      }
      NormP = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(NormP);
      NormP = sqrt(NormP);

      const double zero_tol = 1.0e-12;
      if (NormP < zero_tol)
        NormP = 1.00;

      return NormP;
    }

    virtual bool CheckVelocityConvergence(const double NormDv, double &errorNormDv)
    {
      return false;
    }

    virtual bool CheckPressureConvergence(const double NormDp, double &errorNormDp, double &NormP)
    {
      return false;
    }

    virtual bool FixTimeStepMomentum(const double DvErrorNorm, bool &fixedTimeStep)
    {
      return false;
    }

    virtual bool CheckMomentumConvergence(const double DvErrorNorm, bool &fixedTimeStep)
    {
      return false;
    }

    virtual bool FixTimeStepContinuity(const double DvErrorNorm, bool &fixedTimeStep)
    {
      return false;
    }

    virtual bool CheckContinuityConvergence(const double DvErrorNorm, bool &fixedTimeStep)
    {
      return false;
    }

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

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    // Fractional step index.
    /*  1 : Momentum step (calculate fractional step velocity)
     * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
     * 4 : Pressure step
     * 5 : Computation of projections
     * 6 : End of step velocity
     */
    //    unsigned int mStepId;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    virtual void InitializeStrategy(SolverSettingsType &rSolverConfig)
    {
      KRATOS_TRY;
      KRATOS_CATCH("");
    }

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
    VPStrategy &operator=(VPStrategy const &rOther) {}

    /// Copy constructor.
    VPStrategy(VPStrategy const &rOther) {}

    ///@}

  }; /// Class VPStrategy

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}

  ///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_V_P_STRATEGY_H

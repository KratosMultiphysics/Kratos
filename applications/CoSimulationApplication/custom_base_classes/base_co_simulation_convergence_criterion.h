//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#if !defined(KRATOS_CO_SIMULATION_BASE_CONVERGENCE_CRITERION_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_CONVERGENCE_CRITERION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "includes/model_part.h"

// Application includes
#include "includes/deprecated_variables.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

class CoSimulationBaseConvergenceCriterion
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseConvergenceCriterion);
    typedef double TDataType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseConvergenceCriterion(TDataType RatioTolerance,
                                         TDataType AbsTolerance)
    {
        mRelativeTolerance = RatioTolerance;
        mAbsoluteTolerance = AbsTolerance;
        mInitialNorm = -1;
    }

    virtual ~CoSimulationBaseConvergenceCriterion()
    {
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    virtual bool IsConverged(ModelPart &rModelPart, std::string variableName)
    {
        // Initialize
        TDataType DifferenceNorm = 0.0;
        unsigned int NumDofs = 0;

        // Set a partition for OpenMP
        unsigned int NumNodes = rModelPart.Nodes().size();
        typename ModelPart::NodeIterator itNodeBegin = rModelPart.NodesBegin();

        if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(variableName + "_X"))
        {
        // Loop over Dofs
#pragma omp parallel for reduction(+ \
                                   : DifferenceNorm, NumDofs)
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                ModelPart::NodeIterator node_i = itNodeBegin + i;
                Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(variableName);
                array_1d<double, 3> &nodalResult0 = node_i->FastGetSolutionStepValue(nodalResultVariable, 0);
                array_1d<double, 3> &nodalResult1 = node_i->FastGetSolutionStepValue(nodalResultVariable, 1);

                array_1d<double, 3> nodalResult = nodalResult0 - nodalResult1;

                DifferenceNorm += nodalResult[0] * nodalResult[0] + nodalResult[1] * nodalResult[1] + nodalResult[2] * nodalResult[2];
            }
        }
        else
        {

#pragma omp parallel for reduction(+ \
                                   : DifferenceNorm, NumDofs)
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                ModelPart::NodeIterator node_i = itNodeBegin + i;
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(variableName);
                double &nodalResult0 = node_i->FastGetSolutionStepValue(nodalResultVariable, 0);
                double &nodalResult1 = node_i->FastGetSolutionStepValue(nodalResultVariable, 1);
                double nodalResult = nodalResult0 - nodalResult1;
                DifferenceNorm += nodalResult * nodalResult;
            }
        }

        if (DifferenceNorm == 0.0)
            DifferenceNorm = 1.0;

        if (mInitialNorm < 0)
            mInitialNorm = DifferenceNorm;

        TDataType NormRatio = sqrt(DifferenceNorm / mInitialNorm);

        TDataType DifferenceNormAbs = sqrt(DifferenceNorm) / static_cast<TDataType>(NumDofs);

        std::cout << "CONVERGENCE CHECK:" << std::endl;
        std::cout << variableName << " .: ratio = " << NormRatio << "; exp.ratio = " << mRelativeTolerance << " abs = " << DifferenceNormAbs << " exp.abs = " << mAbsoluteTolerance << std::endl;

        if ((NormRatio <= mRelativeTolerance || DifferenceNormAbs <= mAbsoluteTolerance))
        {
            if (rModelPart.GetCommunicator().MyPID() == 0)
            {
                std::cout << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
            }
            return true;
        }
        else
        {
            return false;
        }
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
    TDataType mAbsoluteTolerance;
    TDataType mRelativeTolerance;
    TDataType mInitialNorm;

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

    ///@}
};
}
#endif //KRATOS_CO_SIMULATION_BASE_CONVERGENCE_CRITERION_H_INCLUDED
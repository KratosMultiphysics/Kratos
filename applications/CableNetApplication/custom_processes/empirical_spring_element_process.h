//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//
//



#ifndef EMPIRICAL_SPRING_ELEMENT_PROCESS_H
#define EMPIRICAL_SPRING_ELEMENT_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "geometries/line_3d_2.h"
#include "includes/model_part.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) EmpiricalSpringElementProcess
    : public Process
{
  public:

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmpiricalSpringElementProcess);

    /// Constructor.
    EmpiricalSpringElementProcess(ModelPart &rModelPart,
     Parameters InputParameters, const Vector &FittedPolynomial):mrModelPart(rModelPart),mParameters(InputParameters),mrFittedPoly(FittedPolynomial)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"           : "example_part",
            "computing_model_part_name" : "computing_domain",
            "node_ids"                  : [1,2],
            "element_id"                : 1,
            "property_id"               : 1,
            "displacement_data"         : [0.0,1.0,2.0,3.0],
            "force_data"                : [0.0,1.0,2.0,3.0],
            "polynomial_order"          : 3
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);

        KRATOS_ERROR_IF(mParameters["node_ids"].size()!=2) << "only two nodes for each spring allowed !" << std::endl;
        KRATOS_ERROR_IF(mParameters["displacement_data"].size()!=mParameters["force_data"].size()) << "only two nodes for each spring allowed !" << std::endl;


        KRATOS_CATCH("")
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        this->CreateEmpiricalSpringElement();
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }



    void CreateEmpiricalSpringElement() const
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

  protected:


  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    const Vector& mrFittedPoly;

}; // Class

}; // namespace

#endif

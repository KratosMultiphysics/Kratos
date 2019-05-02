//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ModelPartSetProperties, KratosSolidMechanicsFastSuite)
{
  Model current_model;

  ModelPart& model_part = current_model.CreateModelPart("Main");

  ModelPart& sub_model_part = model_part.CreateSubModelPart("Part1");

  Properties::Pointer NewProperty = Kratos::make_shared<Properties>(0);

  model_part.AddProperties(NewProperty);
    
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
  KRATOS_CHECK_EQUAL(sub_model_part.NumberOfProperties(), 0);

  sub_model_part.SetProperties(model_part.pProperties());
  
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
  KRATOS_CHECK_EQUAL(sub_model_part.NumberOfProperties(), 1);

  model_part.RemoveSubModelPart("Part1)");
  
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartGetProperties, KratosSolidMechanicsFastSuite)
{
  Model current_model;

  ModelPart& model_part = current_model.CreateModelPart("Main");

  ModelPart& sub_model_part = model_part.CreateSubModelPart("Part1");

  Properties::Pointer MainProperty = model_part.pGetProperties(-1);
    
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
  KRATOS_CHECK_EQUAL(sub_model_part.NumberOfProperties(), 0);

  Properties::Pointer SubProperty = sub_model_part.pGetProperties(-1);
  
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
  KRATOS_CHECK_EQUAL(sub_model_part.NumberOfProperties(), 1);  
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveProperties, KratosSolidMechanicsFastSuite)
{
  Model current_model;

  ModelPart& model_part = current_model.CreateModelPart("Main");

  ModelPart& sub_model_part = model_part.CreateSubModelPart("Part1");

  Properties::Pointer MainProperty = model_part.pGetProperties(0);
  
  sub_model_part.SetProperties(model_part.pProperties());
  
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);
  KRATOS_CHECK_EQUAL(sub_model_part.NumberOfProperties(), 1);

  model_part.RemoveSubModelPart("Part1)");
  
  KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 1);    
}

}
}  // namespace Kratos.

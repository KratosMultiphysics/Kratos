//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes
#include <algorithm>

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"

#include "delaunay_meshing_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, NODE_TO_ELEMENT,        0 );
KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, ELEMENT_TO_NODE,        1 );
KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, ELEMENT_TO_ELEMENT,     2 );
KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, INITIALIZE_MASTER_CONDITION,        3 );
KRATOS_CREATE_LOCAL_FLAG( MeshDataTransferUtilities, MASTER_ELEMENT_TO_MASTER_CONDITION, 4 );

void MeshDataTransferUtilities::TransferData(ModelPart& rModelPart,
                                             const Element & rReferenceElement,
                                             PointPointerVector &list_of_new_centers,
                                             std::vector<Geometry<Node<3> > >& list_of_new_vertices,
                                             Flags Options)

{
  KRATOS_TRY

  ModelPart::MeshesContainerType Meshes = rModelPart.GetMeshes();


  if(Options.Is(MeshDataTransferUtilities::NODE_TO_ELEMENT))
  {
    TransferNodalValuesToElements(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices);
  }
  else
  {
    if(Options.Is(MeshDataTransferUtilities::ELEMENT_TO_ELEMENT)){
      TransferElementalValuesToElements(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices);
    }
    else
    {
      if(Options.Is(MeshDataTransferUtilities::ELEMENT_TO_NODE)){
        TransferElementalValuesToNodes(rModelPart,rReferenceElement,list_of_new_centers,list_of_new_vertices);
      }
    }
  }

  KRATOS_CATCH( "" )
}


void MeshDataTransferUtilities::InitializeBoundaryData(Condition* rCurrentCondition,
                                                       const TransferParameters& rTransferVariables,
                                                       const ProcessInfo& rCurrentProcessInfo)
{

  KRATOS_TRY

  const unsigned int dimension  = rCurrentProcessInfo[SPACE_DIMENSION];
  const unsigned int voigt_size = dimension * (dimension +1) * 0.5; //axisymmetric, processinfo is needed


  BoundaryVariables rVariables;
  rVariables.Initialize(dimension,voigt_size);

  this->TransferInitialBoundaryData(rCurrentCondition, rTransferVariables, rVariables);

  KRATOS_CATCH( "" )

}

void MeshDataTransferUtilities::TransferInitialBoundaryData(Condition* rCurrentCondition,
                                                            const TransferParameters& rTransferVariables,
                                                            BoundaryVariables& rVariables)
{

  KRATOS_TRY

  //double
  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.DoubleVariables[i]),rVariables.DoubleVariable);
  }

  //array_1d
  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.Array1DVariables[i]),rVariables.Array1DVariable);
  }

  //Vector
  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.VectorVariables[i]),rVariables.VectorVariable);
  }

  //Matrix
  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.MatrixVariables[i]),rVariables.MatrixVariable);
  }


  KRATOS_CATCH( "" )

}

//*******************************************************************************************
//*******************************************************************************************

void MeshDataTransferUtilities::TransferCurrentBoundaryData(Element* rCurrentElement,
                                                            Condition* rCurrentCondition,
                                                            const TransferParameters& rTransferVariables,
                                                            BoundaryVariables& rVariables,
                                                            BoundaryVariableArrays& rVariableArrays,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //double
  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
  {
    rCurrentElement->GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),rVariableArrays.DoubleVariableArray,rCurrentProcessInfo);

    //if there is more than one integration point, an average or an interpolation is need
    rVariables.DoubleVariable = rVariableArrays.DoubleVariableArray[0];
    for(unsigned int j=1; j<rVariableArrays.array_size; ++j)
    {
      rVariables.DoubleVariable  += rVariableArrays.DoubleVariableArray[j];
    }
    rVariables.DoubleVariable *= (1.0/double(rVariableArrays.array_size));
    rCurrentCondition->SetValue(*(rTransferVariables.DoubleVariables[i]),rVariables.DoubleVariable);
  }

  //array_1d
  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
  {
    rCurrentElement->GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),rVariableArrays.Array1DVariableArray,rCurrentProcessInfo);

    //if there is more than one integration point, an average or an interpolation is need
    rVariables.Array1DVariable = rVariableArrays.Array1DVariableArray[0];
    for(unsigned int j=1; j<rVariableArrays.array_size; ++j)
    {
      rVariables.Array1DVariable += rVariableArrays.Array1DVariableArray[j];
    }
    rVariables.Array1DVariable *= (1.0/double(rVariableArrays.array_size));
    rCurrentCondition->SetValue(*(rTransferVariables.Array1DVariables[i]),rVariables.Array1DVariable);
  }

  //Vector
  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
  {

    rCurrentElement->GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),rVariableArrays.VectorVariableArray,rCurrentProcessInfo);

    //if there is more than one integration point, an average or an interpolation is need
    rVariables.VectorVariable = rVariableArrays.VectorVariableArray[0];
    for(unsigned int j=1; j<rVariableArrays.array_size; ++j)
    {
      rVariables.VectorVariable  += rVariableArrays.VectorVariableArray[j];
    }
    rVariables.VectorVariable *= (1.0/double(rVariableArrays.array_size));
    rCurrentCondition->SetValue(*(rTransferVariables.VectorVariables[i]),rVariables.VectorVariable);
  }

  //Matrix
  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
  {

    rCurrentElement->GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),rVariableArrays.MatrixVariableArray,rCurrentProcessInfo);

    //if there is more than one integration point, an average or an interpolation is need
    rVariables.MatrixVariable = rVariableArrays.MatrixVariableArray[0];
    for(unsigned int j=1; j<rVariableArrays.array_size; ++j)
    {
      rVariables.MatrixVariable  += rVariableArrays.MatrixVariableArray[i];
    }
    rVariables.MatrixVariable *= (1.0/double(rVariableArrays.array_size));
    rCurrentCondition->SetValue(*(rTransferVariables.MatrixVariables[i]),rVariables.MatrixVariable);
  }


  KRATOS_CATCH( "" )
}


void MeshDataTransferUtilities::TransferBoundaryData(Condition::Pointer rCurrentCondition,
                                                     Condition::Pointer rReferenceCondition,
                                                     const TransferParameters& rTransferVariables)
{
  KRATOS_TRY

  //double
  for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.DoubleVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.DoubleVariables[i])) );
  }

  //array_1d
  for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.Array1DVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.Array1DVariables[i])) );
  }

  //Vector
  for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.VectorVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.VectorVariables[i])) );
  }

  //Matrix
  for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
  {
    rCurrentCondition->SetValue(*(rTransferVariables.MatrixVariables[i]), rReferenceCondition->GetValue(*(rTransferVariables.MatrixVariables[i])) );
  }


  KRATOS_CATCH( "" )
}

void MeshDataTransferUtilities::TransferBoundaryData(Element* rCurrentElement,
                                                     Condition* rCurrentCondition,
                                                     const TransferParameters& rTransferVariables,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  const unsigned int dimension  = rCurrentProcessInfo[SPACE_DIMENSION];
  const unsigned int voigt_size = dimension * (dimension +1) * 0.5;

  BoundaryVariables Variables;
  Variables.Initialize(dimension,voigt_size);


  //initialize to zero
  this->TransferInitialBoundaryData(rCurrentCondition,rTransferVariables,Variables);


  unsigned int integration_points_number = (rCurrentElement->pGetGeometry())->IntegrationPointsNumber(rCurrentElement->GetIntegrationMethod());

  BoundaryVariableArrays VariableArrays;
  VariableArrays.Initialize(integration_points_number);

  //transfer element values to condition
  this->TransferCurrentBoundaryData(rCurrentElement,rCurrentCondition,rTransferVariables,Variables,VariableArrays,rCurrentProcessInfo);


  KRATOS_CATCH( "" )

}

void MeshDataTransferUtilities::TransferBoundaryData(const TransferParameters& rTransferVariables,
                                                     ModelPart& rModelPart)
{
  KRATOS_TRY


  if(rTransferVariables.Options.Is(MeshDataTransferUtilities::INITIALIZE_MASTER_CONDITION))
  {

    //std::cout<<"  TRANSFER INITIALIZE MASTER CONDITION "<<std::endl;

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const unsigned int dimension  = rCurrentProcessInfo[SPACE_DIMENSION];
    const unsigned int voigt_size = dimension * (dimension +1) * 0.5; //axisymmetric, processinfo is needed

    BoundaryVariables Variables;
    Variables.Initialize(dimension,voigt_size);

    //initialize to zero all skin master-conditions = (used instead of master-nodes)
    for(auto& i_cond : rModelPart.Conditions())
    {
      this->TransferInitialBoundaryData(&i_cond, rTransferVariables, Variables);
    }

    //std::cout<<"  TRANSFER DONE "<<std::endl;
  }


  if(rTransferVariables.Options.Is(MeshDataTransferUtilities::MASTER_ELEMENT_TO_MASTER_CONDITION))
  {

    std::cout<<"  TRANSFER MASTER_ELEMENT_TO_MASTER_CONDITION "<<std::endl;

    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const unsigned int dimension  = rCurrentProcessInfo[SPACE_DIMENSION];
    const unsigned int voigt_size = dimension * (dimension +1) * 0.5; //axisymmetric, processinfo is needed

    BoundaryVariables Variables;
    Variables.Initialize(dimension,voigt_size);


    //initialize to zero all skin master-conditions = (used instead of master-nodes)
    for(auto& i_cond : rModelPart.Conditions())
    {
      this->TransferInitialBoundaryData(&i_cond, rTransferVariables, Variables);
    }


    BoundaryVariableArrays VariableArrays;

    //store that information in all body skin if there is a Contact Condition;
    for(auto& i_cond : rModelPart.Conditions())
    {

      if(i_cond.Is(CONTACT) && i_cond.Is(ACTIVE)){

        //std::cout<<" Transfer: Cond: "<<i_cond.Id()<<" is Active "<<std::endl;

        Element*           MasterElement   = i_cond.GetValue(MASTER_ELEMENT).get();
        Condition*         MasterCondition = i_cond.GetValue(MASTER_CONDITION).get();

        unsigned int integration_points_number = (MasterElement->pGetGeometry())->IntegrationPointsNumber(MasterElement->GetIntegrationMethod());

        VariableArrays.Initialize(integration_points_number);

        this->TransferCurrentBoundaryData(MasterElement,MasterCondition,rTransferVariables,Variables,VariableArrays,rCurrentProcessInfo);


        //std::cout<<" MasterCond: "<<MasterCondition->Id()<<" is Active "<<std::endl;

      }

    }


    std::cout<<"  TRANSFER DONE "<<std::endl;

  }


  KRATOS_CATCH( "" )

}


void MeshDataTransferUtilities::TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
                                                              ModelPart& rModelPart)
{

  KRATOS_TRY

  //std::cout<<" [ Data Transfer NODE to ELEMENT ] :"<<std::endl;

  double alpha = 0.25; //[0,1] //smoothing level of the Jacobian

  Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin()->GetGeometry();
  GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
  unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );

  std::vector<double> NodesDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
  std::vector<Vector> NodesVectorVariableArray (integration_points_number);
  std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

  std::vector<double> ElementDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
  std::vector<Vector> ElementVectorVariableArray (integration_points_number);
  std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);


  ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

  for(auto& i_elem : rModelPart.Elements())
  {

    Geometry<Node<3> > & rGeometry = i_elem.GetGeometry();
    IntegrationMethod =  rGeometry.GetDefaultIntegrationMethod();
    integration_points_number = rGeometry.IntegrationPointsNumber( IntegrationMethod );

    const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( IntegrationMethod );

    //shape functions
    Vector N;

    //double
    for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
    {

      //elemental value
      i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);

      std::fill(NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0 );

      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        N = row( Ncontainer, j );

        //nodal value
        for( unsigned int k=0 ; k<rGeometry.size(); ++k)
        {
          NodesDoubleVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i]));
        }

        NodesDoubleVariableArray[j] *= (alpha);
        NodesDoubleVariableArray[j] += (1-alpha) * ElementDoubleVariableArray[j];

      }


      i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),NodesDoubleVariableArray,CurrentProcessInfo);

    }

    //array_1d
    for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
    {

      //elemental value
      i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);

      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        NodesArray1DVariableArray[j].clear();
      }


      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        N = row( Ncontainer, j );

        //nodal value
        for( unsigned int k=0 ; k<rGeometry.size(); ++k)
        {
          NodesArray1DVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i]));
        }

        NodesArray1DVariableArray[j] *= (alpha);
        NodesArray1DVariableArray[j] += (1-alpha) * ElementArray1DVariableArray[j];

      }


      i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),NodesArray1DVariableArray,CurrentProcessInfo);

    }

    //Vector
    for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
    {

      //elemental value
      i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);


      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        if( ElementVectorVariableArray[j].size() != 0 )
          NodesVectorVariableArray[j] = ZeroVector(ElementVectorVariableArray[j].size());
        else
          NodesVectorVariableArray[j] = ZeroVector(0); //¿value?
      }

      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        N = row( Ncontainer, j );

        //nodal value
        for( unsigned int k=0 ; k<rGeometry.size(); ++k)
        {
          NodesVectorVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]));
        }

        NodesVectorVariableArray[j] *= (alpha);
        NodesVectorVariableArray[j] += (1-alpha) * ElementVectorVariableArray[j];

      }


      i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),NodesVectorVariableArray,CurrentProcessInfo);

    }


    //Matrix

    for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
    {

      //elemental value
      i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);


      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        if( ElementMatrixVariableArray[j].size1() != 0 && ElementMatrixVariableArray[j].size2() != 0 )
          NodesMatrixVariableArray[j] = ZeroMatrix(ElementMatrixVariableArray[j].size1(),ElementMatrixVariableArray[j].size2());
        else
          NodesMatrixVariableArray[j] = ZeroMatrix(0,0); //¿value?
      }


      for(unsigned int j=0; j<integration_points_number; ++j)
      {
        N = row( Ncontainer, j );

        //nodal value
        for( unsigned int k=0 ; k<rGeometry.size(); ++k)
        {
          NodesMatrixVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]));
        }

        NodesMatrixVariableArray[j] *= (alpha);
        NodesMatrixVariableArray[j] += (1-alpha) * ElementMatrixVariableArray[j];

      }


      i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),NodesMatrixVariableArray,CurrentProcessInfo);

    }


  }

  //std::cout<<" [ Finished NODE to ELEMENT Transfer ]"<<std::endl;

  KRATOS_CATCH( "" )
}


//*******************************************************************************************
//*******************************************************************************************

//KRATOS MESH INPUT
void MeshDataTransferUtilities::TransferNodalValuesToElements(const TransferParameters& rTransferVariables,
                                                              const Variable<double>& rCriticalVariable,
                                                              const double& rCriticalValue,
                                                              ModelPart& rModelPart)
{

  KRATOS_TRY


  //std::cout<<" [ Data Transfer NODE to ELEMENT ] : based on critical values of "<<rCriticalVariable<<std::endl;
  double alpha = 0.25; //[0,1] //smoothing level of the Jacobian

  Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin()->GetGeometry();
  GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
  unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );


  std::vector<double> ComputedValues(integration_points_number);
  double computed_value=0;
  double critical_value=0;

  ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

  for(auto& i_elem : rModelPart.Elements())
  {
    i_elem.GetValueOnIntegrationPoints(rCriticalVariable,ComputedValues,CurrentProcessInfo);

    computed_value = ComputedValues[0];

    for(unsigned int j=1; j<integration_points_number; ++j)
    {
      computed_value += ComputedValues[j];
    }

    computed_value *= i_elem.GetGeometry().Area()/double(integration_points_number);

    critical_value = rCriticalValue;

    if( computed_value > critical_value )
    {
      for(unsigned int i = 0; i<i_elem.GetGeometry().size(); ++i)
      {
        i_elem.GetGeometry()[i].Set(TO_REFINE);
      }
    }
  }

  std::vector<double> NodesDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
  std::vector<Vector> NodesVectorVariableArray (integration_points_number);
  std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

  std::vector<double> ElementDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
  std::vector<Vector> ElementVectorVariableArray (integration_points_number);
  std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);

  int counter = 0;

  for(auto& i_elem : rModelPart.Elements())
  {

    Geometry<Node<3> >& rGeometry = i_elem.GetGeometry();
    IntegrationMethod = rGeometry.GetDefaultIntegrationMethod();
    integration_points_number = rGeometry.IntegrationPointsNumber( IntegrationMethod );
    const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( IntegrationMethod );

    //shape functions
    Vector N;

    bool apply_smoothing = false;
    for( unsigned int k=0 ; k<rGeometry.size(); ++k)
    {
      if(rGeometry[k].Is(TO_REFINE))
        apply_smoothing = true;
    }

    if( apply_smoothing ){

      counter++;
      //double
      for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
      {

        //elemental value
        i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);

        std::fill(NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0 );

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          N = row( Ncontainer, j );

          //nodal value
          for( unsigned int k=0 ; k<rGeometry.size(); ++k)
          {
            NodesDoubleVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i]));
            //std::cout<<" Node var ["<<k<<"] "<<rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i]))<<std::endl;
          }

          NodesDoubleVariableArray[j] *= (alpha);
          NodesDoubleVariableArray[j] += (1-alpha) * ElementDoubleVariableArray[j];

        }

        //std::cout<<" transfer ["<<i_elem.Id()<<"] "<<NodesDoubleVariableArray[0]<<" element "<<ElementDoubleVariableArray[0]<<std::endl;

        i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),NodesDoubleVariableArray,CurrentProcessInfo);

      }

      //array_1d
      for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
      {

        //elemental value
        i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);
        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          NodesArray1DVariableArray[j].clear();
        }

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          N = row( Ncontainer, j );

          //nodal value
          for( unsigned int k=0 ; k<rGeometry.size(); ++k)
          {
            NodesArray1DVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i]));
          }

          NodesArray1DVariableArray[j] *= (alpha);
          NodesArray1DVariableArray[j] += (1-alpha) * ElementArray1DVariableArray[j];

        }


        i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),NodesArray1DVariableArray,CurrentProcessInfo);

      }

      //Vector
      for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
      {

        //elemental value
        i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          if(ElementVectorVariableArray[j].size() != 0)
            NodesVectorVariableArray[j] = ZeroVector(ElementVectorVariableArray[j].size());
          else
            NodesVectorVariableArray[j] = ZeroVector(0); //¿value?
        }

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          N = row( Ncontainer, j );

          //nodal value
          for( unsigned int k=0 ; k<rGeometry.size(); ++k)
          {
            NodesVectorVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]));
          }

          NodesVectorVariableArray[j] *= (alpha);
          NodesVectorVariableArray[j] += (1-alpha) * ElementVectorVariableArray[j];

        }


        i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),NodesVectorVariableArray,CurrentProcessInfo);

      }


      //Matrix

      for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
      {

        //elemental value
        i_elem.GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          if(ElementMatrixVariableArray[j].size1() !=0 && ElementMatrixVariableArray[j].size2() != 0)
            NodesMatrixVariableArray[j]= ZeroMatrix(ElementMatrixVariableArray[j].size1(),ElementMatrixVariableArray[j].size2());
          else
            NodesMatrixVariableArray[j] = ZeroMatrix(0,0); //¿value?
        }

        for(unsigned int j=0; j<integration_points_number; ++j)
        {
          N = row( Ncontainer, j );

          //nodal value
          for( unsigned int k=0 ; k<rGeometry.size(); ++k)
          {
            NodesMatrixVariableArray[j] += N [k] * rGeometry[k].FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]));
          }

          NodesMatrixVariableArray[j] *= (alpha);
          NodesMatrixVariableArray[j] += (1-alpha) * ElementMatrixVariableArray[j];

        }


        i_elem.SetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),NodesMatrixVariableArray,CurrentProcessInfo);

      }

    }
  }

  for(auto& i_node : rModelPart.Nodes())
  {
    i_node.Reset(TO_REFINE);
  }

  //std::cout<<" [ Finished NODE to ELEMENT Transfer ] : ( Performed "<<counter<<" transfers of "<<rModelPart.NumberOfElements()<<" possible )"<<std::endl;

  KRATOS_CATCH( "" )
}

//KRATOS MESH INPUT
void MeshDataTransferUtilities::TransferElementalValuesToNodes( const TransferParameters& rTransferVariables,
                                                                ModelPart& rModelPart)
{

  KRATOS_TRY

  //std::cout<<" [ Data Transfer ELEMENT to NODE ] :"<<std::endl;

  ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
  NodesContainerType& rNodes      = rModelPart.Nodes();

  Geometry<Node<3> >& rGeom = rModelPart.ElementsBegin()->GetGeometry();
  GeometryData::IntegrationMethod IntegrationMethod =  rGeom.GetDefaultIntegrationMethod();
  unsigned int integration_points_number = rGeom.IntegrationPointsNumber( IntegrationMethod );

  std::vector<double> NodesDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > NodesArray1DVariableArray (integration_points_number);
  std::vector<Vector> NodesVectorVariableArray (integration_points_number);
  std::vector<Matrix> NodesMatrixVariableArray (integration_points_number);

  std::vector<double> ElementDoubleVariableArray (integration_points_number);
  std::vector<array_1d<double,3> > ElementArray1DVariableArray (integration_points_number);
  std::vector<Vector> ElementVectorVariableArray (integration_points_number);
  std::vector<Matrix> ElementMatrixVariableArray (integration_points_number);

  unsigned int buffer_size = rModelPart.GetBufferSize();
  VariablesList& variables_list = rModelPart.GetNodalSolutionStepVariablesList();


  for(auto& i_node : rNodes)
  {

    if(i_node.IsNot(RIGID)){

      //fill variables that are non assigned vectors
      this->FillVectorData(variables_list, i_node);

      ElementWeakPtrVectorType& nElements = i_node.GetValue(NEIGHBOUR_ELEMENTS);

      double Area         = 0;
      double ElementArea  = 0;
      std::fill( NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0);

      //Initialize variables

      //double
      std::fill( NodesDoubleVariableArray.begin(), NodesDoubleVariableArray.end(), 0.0);


      //Array1D
      for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
      {
        NodesArray1DVariableArray[i].clear();
      }


      //Vector
      for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
      {
        (rModelPart.Elements().begin())->GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);
        if( ElementVectorVariableArray[i].size() != 0)
          NodesVectorVariableArray[i] = ZeroVector(ElementVectorVariableArray[i].size());
        else
          NodesVectorVariableArray[i] = ZeroVector(0); //¿value?
      }

      //Matrix
      for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
      {
        (rModelPart.Elements().begin())->GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);
        if( ElementMatrixVariableArray[i].size1() != 0  && ElementMatrixVariableArray[i].size2() != 0 )
          NodesMatrixVariableArray[i] = ZeroMatrix(ElementMatrixVariableArray[i].size1(),ElementMatrixVariableArray[i].size2());
        else
          NodesMatrixVariableArray[i] = ZeroMatrix(0,0); //¿value?
      }


      for(auto& i_nelem : nElements)
      {

        Geometry<Node<3> >& rGeometry = i_nelem.GetGeometry();
        ElementArea = rGeometry.Area();
        Area += ElementArea;

        //double
        for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
        {
          //elemental value
          i_nelem.GetValueOnIntegrationPoints(*(rTransferVariables.DoubleVariables[i]),ElementDoubleVariableArray,CurrentProcessInfo);
          for(unsigned int j=0; j<integration_points_number; ++j)
          {
            NodesDoubleVariableArray[i] += ElementDoubleVariableArray[j] * ElementArea/double(integration_points_number);
          }
        }

        //Array1D
        for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
        {
          //elemental value
          i_nelem.GetValueOnIntegrationPoints(*(rTransferVariables.Array1DVariables[i]),ElementArray1DVariableArray,CurrentProcessInfo);
          for(unsigned int j=0; j<integration_points_number; ++j)
          {
            NodesArray1DVariableArray[i] += ElementArray1DVariableArray[j] * ElementArea/double(integration_points_number);
          }
        }

        //Vector
        for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
        {
          //elemental value
          i_nelem.GetValueOnIntegrationPoints(*(rTransferVariables.VectorVariables[i]),ElementVectorVariableArray,CurrentProcessInfo);
          for(unsigned int j=0; j<integration_points_number; ++j)
          {
            NodesVectorVariableArray[i] += ElementVectorVariableArray[j] * ElementArea/double(integration_points_number);
          }
        }

        //Matrix
        for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
        {
          //elemental value
          i_nelem.GetValueOnIntegrationPoints(*(rTransferVariables.MatrixVariables[i]),ElementMatrixVariableArray,CurrentProcessInfo);
          for(unsigned int j=0; j<integration_points_number; ++j)
          {
            NodesMatrixVariableArray[i] += ElementMatrixVariableArray[j] * ElementArea/double(integration_points_number);
          }
        }

      }

      if(Area!=0){
        //double
        for(unsigned int i=0; i<rTransferVariables.DoubleVariables.size(); ++i)
        {
          NodesDoubleVariableArray[i] /= Area;

          //std::cout<<" node ["<<i_node.Id()<<"] "<<NodesDoubleVariableArray[i]<<" Area "<<Area<<std::endl;
          if( i_node.SolutionStepsDataHas(*(rTransferVariables.DoubleVariables[i]))){
            i_node.FastGetSolutionStepValue(*(rTransferVariables.DoubleVariables[i])) = NodesDoubleVariableArray[i];
          }
          else{
            std::cout<<" ERROR TR: Something Wrong in node ["<<i_node.Id()<<"] : variable "<<*(rTransferVariables.DoubleVariables[i])<<" was not defined "<<std::endl;
          }
        }
        //Array1D
        for(unsigned int i=0; i<rTransferVariables.Array1DVariables.size(); ++i)
        {
          NodesArray1DVariableArray[i] /= Area;

          if( i_node.SolutionStepsDataHas(*(rTransferVariables.Array1DVariables[i]))){
            i_node.FastGetSolutionStepValue(*(rTransferVariables.Array1DVariables[i])) = NodesArray1DVariableArray[i];
          }
          else{
            std::cout<<" ERROR TR: Something Wrong in node ["<<i_node.Id()<<"] : variable "<<*(rTransferVariables.Array1DVariables[i])<<" was not defined "<<std::endl;
          }
        }
        //Vector
        for(unsigned int i=0; i<rTransferVariables.VectorVariables.size(); ++i)
        {
          NodesVectorVariableArray[i] /= Area;

          if( i_node.SolutionStepsDataHas(*(rTransferVariables.VectorVariables[i]))){
            i_node.FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i])) = NodesVectorVariableArray[i];
            //fill buffer if empty
            for(unsigned int step = 1; step<buffer_size; ++step)
            {
              if(i_node.FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step).size() == 0){
                i_node.FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step) = NodesVectorVariableArray[i];
                i_node.FastGetSolutionStepValue(*(rTransferVariables.VectorVariables[i]), step).clear();
              }

            }
          }
          else{
            std::cout<<" ERROR TR: Something Wrong in node ["<<i_node.Id()<<"] : variable "<<*(rTransferVariables.VectorVariables[i])<<" was not defined "<<std::endl;
          }

        }
        //Matrix
        for(unsigned int i=0; i<rTransferVariables.MatrixVariables.size(); ++i)
        {
          NodesMatrixVariableArray[i] /= Area;

          if( i_node.SolutionStepsDataHas(*(rTransferVariables.MatrixVariables[i]))){
            i_node.FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i])) = NodesMatrixVariableArray[i];
            //fill buffer if empty
            for(unsigned int step = 1; step<buffer_size; ++step)
            {
              if(i_node.FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).size1() == 0 &&
                 i_node.FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).size2() == 0 ){
                i_node.FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step) = NodesMatrixVariableArray[i];
                i_node.FastGetSolutionStepValue(*(rTransferVariables.MatrixVariables[i]), step).clear();
              }
            }
          }
          else{
            std::cout<<" ERROR TR: Something Wrong in node ["<<i_node.Id()<<"] : variable "<<*(rTransferVariables.MatrixVariables[i])<<" was not defined "<<std::endl;
          }

        }
      }
      else{
        std::cout<<" ERROR TR: Something Wrong in node ["<<i_node.Id()<<"] : Area = 0 (neighbours: "<<nElements.size()<<") "<<std::endl;
      }


    }

  }


  //std::cout<<" [ Finished ELEMENT to NODE Transfer ]"<<std::endl;

  KRATOS_CATCH( "" )
}


//CENTERS AND NODES INPUT
void MeshDataTransferUtilities::TransferNodalValuesToElements(ModelPart& rModelPart,
                                                              const Element& rReferenceElement,
                                                              PointPointerVector &list_of_new_centers,
                                                              std::vector<Geometry<Node<3> > >& list_of_new_vertices)
{

  KRATOS_TRY

  std::cout<<" [ Data Transfer NODE to ELEMENT: NOT IMPLEMENTED YET ]"<<std::endl;
  std::cout<<" [ Finished NODE to ELEMENT Transfer: NOT IMPLEMENTED YET ]"<<std::endl;

  KRATOS_CATCH( "" )

}


//CENTERS AND NODES INPUT
void MeshDataTransferUtilities::TransferElementalValuesToNodes(ModelPart& rModelPart,
                                                               const Element & rReferenceElement,
                                                               PointPointerVector &list_of_new_centers,
                                                               std::vector<Geometry<Node<3> > >& list_of_new_vertices)
{

  KRATOS_TRY

  std::cout<<" [ Data Transfer ELEMENT to NODE: NOT IMPLEMENTED YET ]"<<std::endl;
  std::cout<<" [ Finished ELEMENT to NODE Transfer: NOT IMPLEMENTED YET ]"<<std::endl;

  KRATOS_CATCH( "" )

}


//CENTERS AND NODES INPUT
void MeshDataTransferUtilities::TransferElementalValuesToElements(ModelPart& rModelPart,
                                                                  const Element & rReferenceElement,
                                                                  PointPointerVector &list_of_new_centers,
                                                                  std::vector<Geometry<Node<3> > >& list_of_new_vertices)
{

  KRATOS_TRY

  //std::cout<<" [ Data Transfer ELEMENT to ELEMENT ]"<<std::endl;

  //definitions for spatial search
  typedef Node<3>                                  PointType;
  typedef Node<3>::Pointer                  PointPointerType;
  typedef std::vector<PointPointerType>   PointPointerVector;
  //typedef std::vector<PointType>             PointTypeVector;
  typedef PointPointerVector::iterator         PointIterator;
  typedef std::vector<double>                 DistanceVector;
  typedef std::vector<double>::iterator     DistanceIterator;
  typedef Bucket<3, PointType, PointPointerVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
  typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
  //definitions for spatial search


  ElementsContainerType& rPreElements = rModelPart.Elements();

  //creating an auxiliary list for the pre integration points
  PointPointerVector list_of_pre_centers;


  //find the center and "radius" of the element
  double xc=0;
  double yc=0;
  double zc=0;
  double radius=0;

  // CREATE LIST OF PREVIOUS ELEMENT CENTERS
  for(auto& i_elem : rPreElements)
  {
    PointsArrayType& vertices=i_elem.GetGeometry().Points();

    if( vertices.size() == 3 ){
      CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(),
                                      vertices[1].X(), vertices[1].Y(),
                                      vertices[2].X(), vertices[2].Y(),
                                      xc,yc,radius);
    }
    else if( vertices.size() == 4 ){
      CalculateCenterAndSearchRadius( vertices[0].X(), vertices[0].Y(), vertices[0].Z(),
                                      vertices[1].X(), vertices[1].Y(), vertices[1].Z(),
                                      vertices[2].X(), vertices[2].Y(), vertices[2].Z(),
                                      vertices[3].X(), vertices[3].Y(), vertices[3].Z(),
                                      xc,yc,zc,radius);
    }
    else{
      KRATOS_THROW_ERROR( std::logic_error, "Wrong Number of Nodes for the Element in Transfer",*this );
    }

    unsigned int id= i_elem.Id();
    PointPointerType p_center =  Kratos::make_intrusive<PointType>(id,xc,yc,zc);

    // if ((*i_elem.base())->GetOptions().Is(Element::THERMAL))
    //   std::cout<<" is thermal "<<std::endl;

    // if ((*i_elem.base())->GetOptions().Is(Element::MECHANICAL))
    //   std::cout<<" is mechanical "<<std::endl;

    list_of_pre_centers.push_back( p_center );

    //std::cout<<" id pre elems "<<list_of_pre_centers.back()->Id()<<std::endl;
  }

  //std::cout<<" list of pre centers "<<list_of_pre_centers.size()<<" list of new centers "<<list_of_new_centers.size()<<std::endl;


  // SEARCH PREVIOUS CENTER NEARER TO THE NEW CENTER

  //creating an auxiliary list for the pre integration points
  unsigned int   bucket_size = 40;
  KdtreeType     nodes_tree(list_of_pre_centers.begin(),list_of_pre_centers.end(),bucket_size);
  double         ResultDistance;

  //make a loop on temporal elements
  ElementsContainerType temporal_elements;
  temporal_elements.reserve(rModelPart.Elements().size());
  temporal_elements.swap(rModelPart.Elements());


  int count=0;
  int moved_transfers = 0;

  unsigned int on_distance = 0;
  unsigned int on_connectivities = 0;
  unsigned int same_selection = 0;

  for(auto& i_center : list_of_new_centers)
  {
    count++;

    //std::cout<<" INTEGRATION POINT  "<<count<<std::endl;

    //find the nearest integration point to the new integration point (work_point)
    PointType& work_point = *i_center;


    // find integration point (Element) which values will be transferred
    PointPointerType result_point = nodes_tree.SearchNearestPoint(work_point,ResultDistance);

    ElementsContainerType::iterator pe; //element corresponding to the result point

    bool extra_search = false;
    if(list_of_new_vertices[i_center->Id()-1].size() == 4 && extra_search) //try no improve transfer in 3D trying to find another result point
    {
      // find integration points in radius (Candidate Elements)
      xc = 0;
      yc = 0;
      zc = 0;
      radius = 0;
      if( list_of_new_vertices[i_center->Id()-1].size() == 3 ){
        CalculateCenterAndSearchRadius( list_of_new_vertices[i_center->Id()-1][0].X(), list_of_new_vertices[i_center->Id()-1][0].Y(),
                                        list_of_new_vertices[i_center->Id()-1][1].X(), list_of_new_vertices[i_center->Id()-1][1].Y(),
                                        list_of_new_vertices[i_center->Id()-1][2].X(), list_of_new_vertices[i_center->Id()-1][2].Y(),
                                        xc,yc,radius);
      }
      else if( list_of_new_vertices[i_center->Id()-1].size() == 4 ){
        CalculateCenterAndSearchRadius( list_of_new_vertices[i_center->Id()-1][0].X(), list_of_new_vertices[i_center->Id()-1][0].Y(), list_of_new_vertices[i_center->Id()-1][0].Z(),
                                        list_of_new_vertices[i_center->Id()-1][1].X(), list_of_new_vertices[i_center->Id()-1][1].Y(), list_of_new_vertices[i_center->Id()-1][1].Z(),
                                        list_of_new_vertices[i_center->Id()-1][2].X(), list_of_new_vertices[i_center->Id()-1][2].Y(), list_of_new_vertices[i_center->Id()-1][2].Z(),
                                        list_of_new_vertices[i_center->Id()-1][3].X(), list_of_new_vertices[i_center->Id()-1][3].Y(), list_of_new_vertices[i_center->Id()-1][3].Z(),
                                        xc,yc,zc,radius);
      }
      else{
        KRATOS_THROW_ERROR( std::logic_error, "Wrong Number of Nodes for the Element in Transfer",*this );
      }

      //std::cout<<" Radius "<<radius<<std::endl;

      unsigned int max_points = 50;
      unsigned int num_points_in_radius = 0;
      PointPointerVector list_of_points_in_radius(max_points);
      DistanceVector list_of_distances(max_points);

      num_points_in_radius = nodes_tree.SearchInRadius(work_point,radius,list_of_points_in_radius.begin(),list_of_distances.begin(),max_points);

      std::vector<int> coincident_nodes;
      int coincident_node = 0;

      for(PointIterator i_candidate = list_of_points_in_radius.begin(); i_candidate!=list_of_points_in_radius.begin()+ num_points_in_radius; ++i_candidate)
      {
        ElementsContainerType::iterator ce = temporal_elements.find( (*i_candidate)->Id() );
        PointsArrayType& vertices          = ( *ce.base() )->GetGeometry().Points();
        coincident_node = 0;
        for(unsigned int i=0; i<list_of_new_vertices[i_center->Id()-1].size(); ++i)
        {
          for(unsigned int j=0; j<vertices.size(); ++j)
          {
            if( list_of_new_vertices[i_center->Id()-1][i].Id() == vertices[j].Id() )
            {
              coincident_node+=1;
            }
          }
        }

        coincident_nodes.push_back(coincident_node);
      }

      //check if there is an element with more coincident connectivities
      int coincident = 0;
      int candidate = 0;
      for(unsigned int j=0; j<coincident_nodes.size(); ++j)
      {
        if(coincident_nodes[j]>coincident){
          candidate = j;
          coincident = coincident_nodes[j];
        }

      }

      //check if is the only one
      //std::sort(coincident_nodes.begin(), coincident_nodes.end(), std::greater<int>());
      int num_candidates = 0;
      std::vector<int> coincident_candidates;
      for(unsigned int j=0; j<coincident_nodes.size(); ++j)
      {
        if( coincident == coincident_nodes[j] ){
          num_candidates += 1;
          coincident_candidates.push_back(j);
        }
      }

      if(num_candidates > 1){
        int distance_candidates = list_of_distances[coincident_candidates[0]];
        candidate = coincident_candidates[0];
        for(unsigned int j=1; j<coincident_candidates.size(); ++j)
        {
          if(list_of_distances[coincident_candidates[j]] < distance_candidates)
            candidate = coincident_candidates[j];
        }
      }

      unsigned int result_id = result_point->Id();

      //select point
      if( result_point->Id() == (*(list_of_points_in_radius.begin()+candidate))->Id() ){
        same_selection += 1;
      }
      else{

        // standart case no radius is needed
        // if( ResultDistance > list_of_distances[candidate] ){
        // 	result_id = (*(list_of_points_in_radius.begin()+candidate))->Id();
        // 	on_connectivities += 1;
        // }
        // else{
        // 	result_id = result_point->Id();
        // 	on_distance += 1;
        // }

        // alternative case if connectivities coincides have a bigger weight
        if( num_candidates > 1 ){
          result_id = result_point->Id();
          on_distance += 1;
        }
        else{
          result_id = (*(list_of_points_in_radius.begin()+candidate))->Id();
          on_connectivities += 1;
        }

      }


      //element selected
      pe = temporal_elements.find( result_id );

      //std::cout<<" ResultPointDistance "<<result_point->Id()<<" ResultPointConnectivities "<<(*(list_of_points_in_radius.begin()+candidate))->Id()<<std::endl;

    }
    else{

      //result point is a previous center
      pe = temporal_elements.find( result_point->Id() );
    }

    Element::Pointer PreviousElement   = ( *pe.base() );
    PointsArrayType& vertices          = (*pe.base())->GetGeometry().Points();

    if( vertices.size() < 3 )
      std::cout<<" Vertices have not the size "<<std::endl;

    if( fabs(ResultDistance) > 1e-30 ){

      // // Connectivities Ids
      // std::cout<<" New Id: "<<i_center->Id()<<" Connectivities : ["<<list_of_new_vertices[i_center->Id()-1][0].Id();
      // for(unsigned int pn=1; pn<list_of_new_vertices[i_center->Id()-1].size(); ++pn){
      //   std::cout<<", "<<list_of_new_vertices[i_center->Id()-1][pn].Id();
      // }
      // std::cout<<"] "<<std::endl;

      // std::cout<<" Pre Id: "<<(*pe.base())->Id()<<" Connectivities : ["<<vertices[0].Id();
      // for(unsigned int pn=1; pn<vertices.size(); ++pn){
      //   std::cout<<", "<<vertices[pn].Id();
      // }
      // std::cout<<"] "<<std::endl;

      // // Connectivities coordinates
      // // std::cout<<" New Id: "<<i_center->Id()<<" Connectivities : ["<<list_of_new_vertices[i_center->Id()-1][0].Coordinates();
      // // for(unsigned int pn=1; pn<vertices.size(); ++pn){
      // //   std::cout<<", "<<list_of_new_vertices[i_center->Id()-1][pn].Coordinates();
      // // }
      // // std::cout<<"] "<<std::endl;

      // // std::cout<<" Pre Id: "<<(*pe.base())->Id()<<" Connectivities : ["<<vertices[0].Coordinates();
      // // for(unsigned int pn=1; pn<vertices.size(); ++pn){
      // //   std::cout<<", "<<vertices[pn].Coordinates();
      // // }
      // // std::cout<<"] "<<std::endl;

      // std::cout<<" New Center :"<<work_point<<std::endl;
      // std::cout<<" Old Center :"<<*result_point<<std::endl;
      // std::cout<<" [Distance: "<<ResultDistance<<"]"<<std::endl;

      moved_transfers ++;
    }


    //***************************************************
    // CREATE ELEMENT WITH THE TRANSFERED VALUES :: START

    //check size and id
    if( i_center->Id()-1 >= list_of_new_vertices.size() )
      KRATOS_ERROR << "list of new vertices out of bounds " << i_center->Id()-1 << " > " <<list_of_new_vertices.size()<<std::endl;

    // inside of the loop create the element, set the variables, and push_back the new element to the model part
    Element::Pointer new_element = (*pe.base())->Clone(i_center->Id(), list_of_new_vertices[i_center->Id()-1]);

    //set transfer variables (commented December 2019:: adding MODEL_PART_NAME to elements, assigned when cloned)
    //new_element->SetValue(MODEL_PART_NAME,vertices[0].GetValue(MODEL_PART_NAME)); //MODEL_PART_NAME set as a variable
    new_element->AssignFlags(*(*pe.base()));

    // In case of interaction of two or more fluids with different properties, the property
    // of the new element is retrivied from the nodes using the PROPERTY_ID variable.
    if (rModelPart.NumberOfProperties() > 1) {
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        GeometryType& r_geometry = new_element->GetGeometry();
        std::vector<int> array_of_properties;
        // unsigned int max_count = 1, curr_count = 1;

        for (unsigned int i = 0; i < list_of_new_vertices[i_center->Id() - 1].size(); i++) {
            if (r_geometry[i].IsNot(RIGID)) {
                array_of_properties.push_back(r_geometry[i].FastGetSolutionStepValue(PROPERTY_ID, 0));
            }
        }
        std::sort(array_of_properties.begin(), array_of_properties.end());
        unsigned int property_id = array_of_properties[0];

        // for (unsigned int i = 0; i < array_of_properties.size(); i++) {
        //    if (array_of_properties[i+1] == array_of_properties[i]) {
        //        curr_count++;
        //    } else {
        //        if (curr_count > max_count) {
        //            max_count = curr_count;
        //            property_id = array_of_properties[i];
        //        }
        //        curr_count = 1;
        //    }
        //}
        // if (curr_count > max_count) {
        //    property_id = array_of_properties.back();
        //}
        Properties::Pointer p_new_property = rModelPart.pGetProperties(property_id);
        new_element->SetProperties(p_new_property);
    }

    //check
    //new_element->PrintInfo(std::cout);

    //setting new elements
    (rModelPart.Elements()).push_back(new_element);

    // CREATE ELEMENT WITH THE TRANSFERED VALUES :: END
    //***************************************************


    //ELEMENT TRANSFER CHECK//
    // std::cout<<" Write Stress [element: "<<new_element->Id()<<"]: ";
    // std::cout.flush();
    // std::vector<Matrix >  StressMatrix;
    // Matrix Stress;
    // Stress.clear();
    // StressMatrix.push_back(Stress);
    // ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
    // new_element->GetValueOnIntegrationPoints(CAUCHY_STRESS_TENSOR,StressMatrix,CurrentProcessInfo);
    // std::cout<<StressMatrix[0]<<std::endl;
    //ELEMENT TRANSFER CHECK//
  }

  temporal_elements.clear();
  //std::cout<<" on distances "<<on_distance<<" on connectitivies "<<on_connectivities<<" same selection "<<same_selection<<std::endl;
  //std::cout<<" [ MOVED TRANSFERS: "<<moved_transfers<<" ]"<<std::endl;
  //std::cout<<" [ Finished ELEMENT to ELEMENT Transfer ]"<<std::endl;


  KRATOS_CATCH( "" )
}


VariablesListDataValueContainer MeshDataTransferUtilities::InterpolateVariables( Geometry<Node<3> > &geom,
                                                                                 const std::vector<double>& N,
                                                                                 VariablesList& rVariablesList,
                                                                                 Node<3>::Pointer pnode,
                                                                                 double alpha )
{

  KRATOS_TRY

  //Copy Variables List
  // std::cout<<" node["<<pnode->Id()<<"] Data "<<(pnode)->SolutionStepData()<<std::endl;

  VariablesListDataValueContainer PreviousVariablesListData = (pnode)->SolutionStepData();

  // std::cout<<" CopiedData "<<PreviousVariablesListData<<std::endl;


  Interpolate( geom, N, rVariablesList, pnode, alpha);

  VariablesListDataValueContainer CurrentVariablesListData = (pnode)->SolutionStepData();

  (pnode)->SolutionStepData() = PreviousVariablesListData;

  // std::cout<<" PreviousVariables "<<PreviousVariablesListData<<std::endl;
  // std::cout<<" CurrentVariables "<<CurrentVariablesListData<<std::endl;

  return CurrentVariablesListData;

  KRATOS_CATCH( "" )

}


void MeshDataTransferUtilities::FillVectorData(VariablesList& rVariablesList,
                                               Node<3>& rNode)
{

  KRATOS_TRY

  unsigned int buffer_size = rNode.GetBufferSize();

  for(const auto& i_variable : rVariablesList)
  {
    std::string variable_name = i_variable.Name();
    if(KratosComponents<Variable<Vector > >::Has(variable_name))
    {
      //std::cout<<"Vector"<<std::endl;
      const Variable<Vector>& variable = KratosComponents<Variable<Vector > >::Get(variable_name);
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //getting the data of the solution step
        Vector& node_data = rNode.FastGetSolutionStepValue(variable, step);

        if( node_data.size() == 0 ){
          node_data = ZeroVector(1);
        }
      }
    }

  }

  KRATOS_CATCH( "" )
}

void MeshDataTransferUtilities::Interpolate( Geometry<Node<3> > &geom,
                                             const std::vector<double>& N,
                                             VariablesList& rVariablesList,
                                             Node<3>::Pointer pnode,
                                             double alpha )
{

  KRATOS_TRY

  unsigned int buffer_size = pnode->GetBufferSize();

  bool all_null = true;
  for(unsigned int i=0; i<N.size(); ++i)
    if(N[i]!=0)
      all_null = false;

  if( all_null )
    KRATOS_THROW_ERROR( std::logic_error,"SOMETHING is wrong with the Interpolation Functions", "" )

  for(const auto& i_variable : rVariablesList)
  {
    std::string variable_name = i_variable.Name();
    double data;
    if(KratosComponents<Variable<double> >::Has(variable_name))
    {
      const Variable<double>& variable = KratosComponents<Variable<double> >::Get(variable_name);
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //getting the data of the solution step
        double& node_data = pnode->FastGetSolutionStepValue(variable, step);

        std::vector<double* > nodes_data;
        for(unsigned int i=0; i<geom.size(); ++i)
          nodes_data.push_back(&geom[i].FastGetSolutionStepValue(variable, step));

        if(alpha != 1 ){

          data = node_data * (1-alpha);
          for(unsigned int i=0; i<geom.size(); ++i)
            data += (alpha) * (N[i]*(*nodes_data[i]));

        }
        else{

          data = (N[0]*(*nodes_data[0]));
          for(unsigned int i=1; i<geom.size(); ++i)
            data += (N[i]*(*nodes_data[i]));
        }

        node_data = data;
      }
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
      //std::cout<<"array1d"<<std::endl;
      const Variable<array_1d<double, 3> >& variable = KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name);
      array_1d<double, 3> data;
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //getting the data of the solution step
        array_1d<double, 3>& node_data = pnode->FastGetSolutionStepValue(variable, step);

        std::vector<array_1d<double, 3>* > nodes_data;
        for(unsigned int i=0; i<geom.size(); ++i)
          nodes_data.push_back(&geom[i].FastGetSolutionStepValue(variable, step));

        if(alpha != 1 ){

          data = node_data * (1-alpha);
          for(unsigned int i=0; i<geom.size(); ++i)
            data += (alpha) * (N[i]*(*nodes_data[i]));

        }
        else{

          data = (N[0]*(*nodes_data[0]));
          for(unsigned int i=1; i<geom.size(); ++i)
            data += (N[i]*(*nodes_data[i]));
        }

        node_data = data;
      }

    }
    else if(KratosComponents<Variable<Matrix > >::Has(variable_name))
    {
      //std::cout<<"Matrix"<<std::endl;
      const Variable<Matrix>& variable = KratosComponents<Variable<Matrix > >::Get(variable_name);
      Matrix data;
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //getting the data of the solution step
        Matrix& node_data = pnode->FastGetSolutionStepValue(variable, step);

        std::vector<Matrix* > nodes_data;
        for(unsigned int i=0; i<geom.size(); ++i)
          nodes_data.push_back(&geom[i].FastGetSolutionStepValue(variable, step));

        if( node_data.size1() > 0 && node_data.size2() ){
          bool same_size = true;

          for(unsigned int i=0; i<geom.size(); ++i)
            if(node_data.size1() != (*nodes_data[i]).size1() && node_data.size2() != (*nodes_data[i]).size2())
              same_size = false;

          if( same_size ) {

            if(alpha != 1 ){
              data = node_data * (1-alpha);
              for(unsigned int i=0; i<geom.size(); ++i)
                data += (alpha) * (N[i]*(*nodes_data[i]));

            }
            else{

              data = (N[0]*(*nodes_data[0]));
              for(unsigned int i=1; i<geom.size(); ++i)
                data += (N[i]*(*nodes_data[i]));
            }

            node_data = data;

          }
        }
      }

    }
    else if(KratosComponents<Variable<Vector > >::Has(variable_name))
    {
      //std::cout<<"Vector"<<std::endl;
      const Variable<Vector>& variable = KratosComponents<Variable<Vector > >::Get(variable_name);
      Vector data;
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //getting the data of the solution step
        Vector& node_data = pnode->FastGetSolutionStepValue(variable, step);

        std::vector<Vector* > nodes_data;
        for(unsigned int i=0; i<geom.size(); ++i)
          nodes_data.push_back(&geom[i].FastGetSolutionStepValue(variable, step));

        // std::cout<<" node ["<<pnode->Id()<<"]"<<std::endl;
        // std::cout<<" variable "<<variable<<std::endl;
        // std::cout<<" node data "<<node_data<<std::endl;
        if( node_data.size() > 0 ){
          bool same_size = true;

          for(unsigned int i=0; i<geom.size(); ++i)
            if(node_data.size() != (*nodes_data[i]).size())
              same_size = false;

          if( same_size ) {

            if(alpha != 1 ){

              data = node_data * (1-alpha);
              for(unsigned int i=0; i<geom.size(); ++i)
                data += (alpha) * (N[i]*(*nodes_data[i]));

            }
            else{

              data = (N[0]*(*nodes_data[0]));
              for(unsigned int i=1; i<geom.size(); ++i)
                data += (N[i]*(*nodes_data[i]));

            }
            node_data = data;

          }
        }

      }
    }
    else if(KratosComponents<Variable<std::string > >::Has(variable_name))
    {
      //std::cout<<"string"<<std::endl;
      //NO INTERPOLATION

      const Variable<std::string>& variable = KratosComponents<Variable<std::string> >::Get(variable_name);
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //assign data from the first node
        pnode->FastGetSolutionStepValue(variable, step) = geom[0].FastGetSolutionStepValue(variable, step);
      }
    }
    else if(KratosComponents<Variable<int > >::Has(variable_name))
    {
      //std::cout<<"int"<<std::endl;
      //NO INTERPOLATION
      const Variable<int>& variable = KratosComponents<Variable<int> >::Get(variable_name);
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //assign data from the first node
        pnode->FastGetSolutionStepValue(variable, step) = geom[0].FastGetSolutionStepValue(variable, step);
      }
    }
    else if(KratosComponents<Variable<bool > >::Has(variable_name))
    {
      //std::cout<<"bool"<<std::endl;
      //NO INTERPOLATION
      const Variable<bool>& variable = KratosComponents<Variable<bool> >::Get(variable_name);
      for(unsigned int step = 0; step<buffer_size; ++step)
      {
        //assign data from the first node
        pnode->FastGetSolutionStepValue(variable, step) = geom[0].FastGetSolutionStepValue(variable, step);
      }
    }
  }

  KRATOS_CATCH( "" )

}

//doubles only
void MeshDataTransferUtilities::InterpolateData( Geometry<Node<3> >& geom,
                                                 const std::vector<double>& N,
                                                 unsigned int step_data_size,
                                                 Node<3>::Pointer pnode,
                                                 double alpha )
{

  KRATOS_TRY

  unsigned int buffer_size = pnode->GetBufferSize();

  //alpha [0,1] //smoothing level of the interpolation
  double data;
  for(unsigned int step = 0; step<buffer_size; ++step)
  {
    //getting the data of the solution step
    double* step_data = (pnode)->SolutionStepData().Data(step);

    std::vector<double* > nodes_data;
    for(unsigned int i=0; i<geom.size(); ++i)
    {
      nodes_data.push_back(geom[i].SolutionStepData().Data(step));
    }

    //copying this data in the position of the vector we are interested in
    for(unsigned int j= 0; j<step_data_size; ++j)
    {
      data = step_data[j] * (1-alpha);
      for(unsigned int i=0; i<geom.size(); ++i)
        data += (alpha) * (N[i]*nodes_data[i][j]);
      step_data[j] = data;
    }
  }

  KRATOS_CATCH( "" )
}


VariablesListDataValueContainer MeshDataTransferUtilities::InterpolateVariablesData( Geometry<Node<3> >& geom,
                                                                                     const std::vector<double>& N,
                                                                                     unsigned int step_data_size,
                                                                                     Node<3>::Pointer pnode,
                                                                                     double alpha )

{
  KRATOS_TRY

  //Copy Variables List
  VariablesListDataValueContainer PreviousVariablesListData = (pnode)->SolutionStepData();

  InterpolateData( geom, N, step_data_size, pnode, alpha );

  VariablesListDataValueContainer CurrentVariablesListData = (pnode)->SolutionStepData();

  (pnode)->SolutionStepData() = PreviousVariablesListData;

  //std::cout<<" PreviousVariables "<<PreviousVariablesListData<<std::endl;
  //std::cout<<" CurrentVariables "<<CurrentVariablesListData<<std::endl;

  return CurrentVariablesListData;

  KRATOS_CATCH( "" )

}


}  // namespace Kratos.

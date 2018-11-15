//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#include "includes/model_part.h"

#if !defined(KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_H_INCLUDED )
#define  KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_H_INCLUDED


// System includes


// External includes


// Project includes
#include "delaunay_meshing_application_variables.h"
#include "custom_utilities/mesher_utilities.hpp"

namespace Kratos
{

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

/// Short class definition.
/// Tool to evaluate the normals on nodes based on the normals of a set of surface conditions
class BoundaryNormalsCalculationUtilities
{
 public:

  ///@name Type Definitions
  ///@{

  typedef ModelPart::NodesContainerType               NodesArrayType;
  typedef ModelPart::ElementsContainerType     ElementsContainerType;
  typedef ModelPart::ConditionsContainerType ConditionsContainerType;
  typedef ModelPart::MeshType                               MeshType;

  ///@}
  ///@name Life Cycle
  ///@{
  ///@}
  ///@name Operators
  ///@{
  ///@}
  ///@name Operations
  ///@{


  /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area).
  /** This is done on the base of the Conditions provided which should be
   * understood as the surface elements of the area of interest.
   * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
   * @param dimension Spatial dimension (2 or 3)
   * @note Use this fuction instead of its overload taking a Conditions array for MPI applications,
   * as it will take care of communication between partitions.
   */

  /// Calculates the area normal (unitary vector oriented as the normal).

  void CalculateUnitBoundaryNormals(ModelPart& rModelPart, int EchoLevel = 0)
  {

    mEchoLevel = EchoLevel;

    if( !rModelPart.IsSubModelPart() ){

      this->ResetBodyNormals(rModelPart); //clear boundary normals

      //FLUID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //standard assignation // fails in sharp edges angle<90
          AddNormalsToNodes(*i_mp);
        }
      }
      //SOLID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //standard assignation // fails in sharp edges angle<90
          AddNormalsToNodes(*i_mp);
        }
      }
      //RIGID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(RIGID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //standard assignation // fails in sharp edges angle<90
          AddNormalsToNodes(*i_mp);
        }
      }
    }
    else{

      this->ResetBodyNormals(rModelPart); //clear boundary normals

      CalculateBoundaryNormals(rModelPart);

      //standard assignation // fails in sharp edges angle<90
      AddNormalsToNodes(rModelPart);

    }


    // For MPI: correct values on partition boundaries
    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

  }


  /// Calculates the area normal (unitary vector oriented as the normal) and weight the normal to shrink
  void CalculateWeightedBoundaryNormals(ModelPart& rModelPart, int EchoLevel = 0)
  {

    mEchoLevel = EchoLevel;

    if( !rModelPart.IsSubModelPart() ){

      this->ResetBodyNormals(rModelPart); //clear boundary normals

      //FLUID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //assignation for solid boundaries : Unity Normals on nodes and Shrink_Factor on nodes
          AddWeightedNormalsToNodes(*i_mp);
        }
      }
      //SOLID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //assignation for solid boundaries : Unity Normals on nodes and Shrink_Factor on nodes
          AddWeightedNormalsToNodes(*i_mp);
        }
      }
      //RIGID Domains
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin(); i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(RIGID) && i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) && i_mp->IsNot(CONTACT) ){

          CalculateBoundaryNormals(*i_mp);

          //assignation for solid boundaries : Unity Normals on nodes and Shrink_Factor on nodes
          AddWeightedNormalsToNodes(*i_mp);
        }
      }
    }
    else{

      this->ResetBodyNormals(rModelPart); //clear boundary normals

      CalculateBoundaryNormals(rModelPart);

      //assignation for solid boundaries: Unity Normals on nodes and Shrink_Factor on nodes
      AddWeightedNormalsToNodes(rModelPart);
    }

    // For MPI: correct values on partition boundaries
    rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
    rModelPart.GetCommunicator().AssembleCurrentData(SHRINK_FACTOR);
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

  int mEchoLevel;

  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{
  ///@}

  /// Calculates the normals of the BOUNDARY for a given model part
  void CalculateBoundaryNormals(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    if( rModelPart.NumberOfConditions() && this->CheckConditionsLocalSpace(rModelPart, dimension-1) ){

      if(mEchoLevel > 0)
        std::cout<<"   ["<<rModelPart.Name()<<"] (BC)"<<std::endl;

      ConditionsContainerType& rConditions = rModelPart.Conditions();

      this->CalculateBoundaryNormals(rConditions);

    }
    else if( rModelPart.NumberOfElements() ){

      if( this->CheckElementsLocalSpace(rModelPart, dimension) ){

        if(mEchoLevel > 0)
          std::cout<<"   ["<<rModelPart.Name()<<"] (BVE)"<<std::endl;

        this->CalculateVolumeBoundaryNormals(rModelPart);

      }
      else if( this->CheckElementsLocalSpace(rModelPart, dimension-1) ){

        if(mEchoLevel > 0)
          std::cout<<"   ["<<rModelPart.Name()<<"] (BE)"<<std::endl;

        ElementsContainerType& rElements = rModelPart.Elements();

        this->CalculateBoundaryNormals(rElements);

      }

    }

    KRATOS_CATCH( "" )

  }


  //this function adds the Contribution of one of the geometries
  //to the corresponding nodes
  static void CalculateUnitNormal2D(ConditionsContainerType::iterator it, array_1d<double,3>& An)
  {
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    //Attention: normal criterion changed for line conditions (vs line elements)
    An[0] =    pGeometry[1].Y()-pGeometry[0].Y();
    An[1] =  -(pGeometry[1].X()-pGeometry[0].X());
    An[2] =    0.00;

    array_1d<double,3>& normal = (it)->GetValue(NORMAL);
    noalias(normal) = An/norm_2(An);
  }

  static void CalculateUnitNormal3D(ConditionsContainerType::iterator it, array_1d<double,3>& An,
                                     array_1d<double,3>& v1,array_1d<double,3>& v2 )
  {
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);

    array_1d<double,3>& normal = (it)->GetValue(NORMAL);

    noalias(normal) = An/norm_2(An);

  }


  //this function adds the Contribution of one of the geometries
  //to the corresponding nodes
  static void CalculateUnitNormal2D(ElementsContainerType::iterator it, array_1d<double,3>& An)
  {
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    if(pGeometry.size()<2){
      std::cout<<" Warning 2D geometry with only "<<pGeometry.size()<<" node :: multiple normal definitions "<<std::endl;
      (it)->GetValue(NORMAL).clear();
    }
    else{

      An[0] = -(pGeometry[1].Y()-pGeometry[0].Y());
      An[1] =   pGeometry[1].X()-pGeometry[0].X();
      An[2] =   0.00;

      array_1d<double,3>& normal = (it)->GetValue(NORMAL);
      noalias(normal) = An/norm_2(An);

    }

  }

  static void CalculateUnitNormal3D(ElementsContainerType::iterator it, array_1d<double,3>& An,
                                     array_1d<double,3>& v1,array_1d<double,3>& v2 )
  {
    Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

    if(pGeometry.size()<3){
      std::cout<<" Warning 3D geometry with only "<<pGeometry.size()<<" nodes :: multiple normal definitions "<<std::endl;
      (it)->GetValue(NORMAL).clear();
    }
    else{

      v1[0] = pGeometry[1].X() - pGeometry[0].X();
      v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
      v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

      v2[0] = pGeometry[2].X() - pGeometry[0].X();
      v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
      v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

      MathUtils<double>::CrossProduct(An,v1,v2);
      An *= 0.5;

      array_1d<double,3>& normal = (it)->GetValue(NORMAL);

      noalias(normal) = An/norm_2(An);

    }
  }


  void ResetBodyNormals(ModelPart& rModelPart)
  {
    KRATOS_TRY

    //resetting the normals
    for(NodesArrayType::iterator in = rModelPart.NodesBegin();
        in !=rModelPart.NodesEnd(); ++in)
    {
      (in->GetSolutionStepValue(NORMAL)).clear();
    }

    KRATOS_CATCH( "" )
  }

  void CheckBodyNormals(ModelPart& rModelPart)
  {
    KRATOS_TRY

    //resetting the normals
    for(NodesArrayType::iterator in = rModelPart.NodesBegin();
        in !=rModelPart.NodesEnd(); ++in)
    {
      std::cout<<" ID: "<<in->Id()<<" normal: "<<(in->GetSolutionStepValue(NORMAL))<<std::endl;
    }

    KRATOS_CATCH( "" )
  }

  /// Calculates the "area normal" (vector oriented as the normal with a magnitude proportional to the area).
  /** This is done on the base of the Conditions provided which should be
   * understood as the surface elements of the area of interest.
   * @param rConditions A set of conditions defining the "skin" of a model
   * @param dimension Spatial dimension (2 or 3)
   * @note This function is not recommended for distributed (MPI) runs, as
   * the user has to ensure that the calculated normals are assembled between
   * processes. The overload of this function that takes a ModelPart is
   * preferable in this case, as it performs the required communication.
   */
  void CalculateBoundaryNormals(ConditionsContainerType& rConditions)
  {
    KRATOS_TRY

    //resetting the normals
    // for(ConditionsContainerType::iterator it =  rConditions.begin();
    //     it != rConditions.end(); ++it)
    //   {
    //     Element::GeometryType& rNodes = it->GetGeometry();
    //     for(unsigned int in = 0; in<rNodes.size(); ++in)
    //       ((rNodes[in]).GetSolutionStepValue(NORMAL)).clear();
    //   }

    const unsigned int dimension = (rConditions.begin())->GetGeometry().WorkingSpaceDimension();

    //std::cout<<" condition geometry: "<<(rConditions.begin())->GetGeometry()<<std::endl;

    //calculating the normals and storing on the conditions
    array_1d<double,3> An;
    if(dimension == 2)
    {
      for(ConditionsContainerType::iterator it =  rConditions.begin(); it !=rConditions.end(); ++it)
      {
        if(it->IsNot(CONTACT) && it->Is(BOUNDARY) )
          CalculateUnitNormal2D(it,An);
      }

    }
    else if(dimension == 3)
    {
      array_1d<double,3> v1;
      array_1d<double,3> v2;
      for(ConditionsContainerType::iterator it =  rConditions.begin(); it !=rConditions.end(); ++it)
      {
        //calculate the normal on the given condition
        if(it->IsNot(CONTACT) && it->Is(BOUNDARY)){
          CalculateUnitNormal3D(it,An,v1,v2);
        }
      }
    }

    KRATOS_CATCH( "" )
  }


  /// Calculates the "area normal" (vector oriented as the normal with a magnitude proportional to the area).
  /** This is done on the base of the Elements provided which should be
   * understood as the surface elements of the area of interest.
   * @param rElements A set of elements defining the "skin" of a model
   * @param dimension Spatial dimension (2 or 3)
   * @note This function is not recommended for distributed (MPI) runs, as
   * the user has to ensure that the calculated normals are assembled between
   * processes. The overload of this function that takes a ModelPart is
   * preferable in this case, as it performs the required communication.
   */
  void CalculateBoundaryNormals(ElementsContainerType& rElements)

  {
    KRATOS_TRY

    //resetting the normals
    // for(ElementsContainerType::iterator it =  rElements.begin(); it != rElements.end(); ++it)
    //   {
    //     Element::GeometryType& rNodes = it->GetGeometry();
    //     for(unsigned int in = 0; in<rNodes.size(); ++in)
    //       ((rNodes[in]).GetSolutionStepValue(NORMAL)).clear();
    //   }

    const unsigned int dimension = (rElements.begin())->GetGeometry().WorkingSpaceDimension();

    //std::cout<<" element geometry: "<<(rElements.begin())->GetGeometry()<<std::endl;

    //calculating the normals and storing on elements
    array_1d<double,3> An;
    if(dimension == 2)
    {
      for(ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); ++it)
      {
        if(it->IsNot(CONTACT)){
          it->Set(BOUNDARY); //give an error in set flags (for the created rigid body)
          CalculateUnitNormal2D(it,An);
        }
      }
    }
    else if(dimension == 3)
    {
      array_1d<double,3> v1;
      array_1d<double,3> v2;
      for(ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); ++it)
      {
        //calculate the normal on the given surface element
        if(it->IsNot(CONTACT)){
          it->Set(BOUNDARY); //give an error in set flags (for the created rigid body)
          CalculateUnitNormal3D(it,An,v1,v2);
        }
      }
    }

    KRATOS_CATCH( "" )
  }


  //Check if the mesh has elements with the spatial dimension
  bool CheckElementsDimension(ModelPart& rModelPart, unsigned int dimension)
  {
    KRATOS_TRY

    ElementsContainerType& rElements = rModelPart.Elements();

    ElementsContainerType::iterator it =  rElements.begin();

    if( (it)->GetGeometry().Dimension() == dimension ){
      return true;
    }
    else{
      return false;
    }

    KRATOS_CATCH( "" )
  }


  //Check if the mesh has boundary conditions with the spatial dimension
  bool CheckConditionsDimension(ModelPart& rModelPart, unsigned int dimension)
  {
    KRATOS_TRY

    ConditionsContainerType& rConditions = rModelPart.Conditions();

    ConditionsContainerType::iterator it =  rConditions.begin();

    if( (it)->GetGeometry().Dimension() == dimension ){
      return true;
    }
    else{
      return false;
    }

    KRATOS_CATCH( "" )
  }


  //Check if the elements local space is equal to the spatial dimension
  bool CheckElementsLocalSpace(ModelPart& rModelPart, unsigned int dimension)
  {
    KRATOS_TRY

    ElementsContainerType& rElements = rModelPart.Elements();

    ElementsContainerType::iterator it =  rElements.begin();

    if( (it)->GetGeometry().LocalSpaceDimension() == dimension ){
      return true;
    }
    else{
      return false;
    }

    KRATOS_CATCH( "" )
  }


  //Check if the conditions local space is equal to the spatial dimension
  bool CheckConditionsLocalSpace(ModelPart& rModelPart, unsigned int dimension)
  {
    KRATOS_TRY

    ConditionsContainerType& rConditions = rModelPart.Conditions();

    ConditionsContainerType::iterator it =  rConditions.begin();

    if( (it)->GetGeometry().LocalSpaceDimension() == dimension ){
      return true;
    }
    else{
      return false;
    }

    KRATOS_CATCH( "" )
  }


  /// Calculates the normals of the BOUNDARY nodes given a mesh
  //  using a consistent way: SOTO & CODINA
  //  fails in sharp edges angle<90
  void CalculateVolumeBoundaryNormals(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    //Reset normals
    ModelPart::NodesContainerType&    rNodes = rModelPart.Nodes();
    ModelPart::ElementsContainerType& rElems = rModelPart.Elements();

    //Check if the neigbours search is already done and set
    bool neighsearch=false;
    unsigned int number_of_nodes = rElems.begin()->GetGeometry().PointsNumber();
    for(unsigned int i=0; i<number_of_nodes; ++i)
      if( (rElems.begin()->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS)).size() > 1 )
        neighsearch=true;

    if( !neighsearch )
      std::cout<<" WARNING :: Neighbour Search Not PERFORMED "<<std::endl;

    for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
    {
      (in->GetSolutionStepValue(NORMAL)).clear();

      if(!neighsearch){
        //*************  Neigbours of nodes search  ************//
        WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
        rE.erase(rE.begin(),rE.end() );
        in->Reset(BOUNDARY);
        //*************  Neigbours of nodes search ************//
      }
    }


    if(!neighsearch){
      //*************  Neigbours of nodes search ************//
      //add the neighbour elements to all the nodes in the mesh
      for(ModelPart::ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ++ie)
      {
        Element::GeometryType& pGeom = ie->GetGeometry();
        for(unsigned int i = 0; i < pGeom.size(); ++i)
        {
          (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
        }
      }
      //*************  Neigbours of nodes search ************//
    }

    //calculating the normals and storing it on nodes

    Vector An(3);
    Element::IntegrationMethod mIntegrationMethod = Element::GeometryDataType::GI_GAUSS_1; //one gauss point
    int PointNumber = 0; //one gauss point
    Matrix J;
    Matrix InvJ;
    double detJ;
    Matrix DN_DX;

    unsigned int assigned       = 0;
    unsigned int not_assigned   = 0;
    unsigned int boundary_nodes = 0;
    //int boundarycounter=0;
    for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); ++in)
    {
      noalias(An) = ZeroVector(3);

      //if(in->Is(BOUNDARY)){

      WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);

      for(WeakPointerVector<Element >::iterator ie= rE.begin(); ie!=rE.end(); ++ie)
      {

        Element::GeometryType& rGeometry = ie->GetGeometry();

        if( rGeometry.EdgesNumber() > 1 &&  rGeometry.LocalSpaceDimension() == dimension ){


          //********** Compute the element integral ******//
          const Element::GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( mIntegrationMethod );
          const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = rGeometry.ShapeFunctionsLocalGradients( mIntegrationMethod );


          J.resize(dimension, dimension, false);
          J = rGeometry.Jacobian( J, PointNumber , mIntegrationMethod );

          InvJ.clear();
          detJ=0;
          //Calculating the inverse of the jacobian and the parameters needed
          MathUtils<double>::InvertMatrix( J, InvJ, detJ);


          //Compute cartesian derivatives for one gauss point
          DN_DX = prod( DN_De[PointNumber] , InvJ );

          double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;


          for(unsigned int i = 0; i < rGeometry.size(); ++i)
          {
            if(in->Id() == rGeometry[i].Id()){

              for(unsigned int d=0; d<dimension; ++d)
              {
                An[d] += DN_DX(i,d) * IntegrationWeight;
              }
            }
          }

          //********** Compute the element integral ******//

        }

        if(norm_2(An)>1e-12){
          noalias(in->FastGetSolutionStepValue(NORMAL)) = An/norm_2(An);
          assigned +=1;
          if(!neighsearch){
            in->Set(BOUNDARY);
          }
        }
        else{
          (in->FastGetSolutionStepValue(NORMAL)).clear();
          //std::cout<<" ERROR: normal not set "<<std::endl;
          not_assigned +=1;
          if(!neighsearch){
            in->Set(BOUNDARY,false);
          }
        }

      }

      if(in->Is(BOUNDARY))
        boundary_nodes +=1;

      //boundarycounter++;
      //}
    }


    if(mEchoLevel > 0)
      std::cout<<"  [ Boundary_Normals  (Mesh Nodes:"<<rNodes.size()<<")[Boundary nodes: "<<boundary_nodes<<" (SET:"<<assigned<<" / NOT_SET:"<<not_assigned<<")] ]"<<std::endl;


    //std::cout<<" Boundary COUNTER "<<boundarycounter<<std::endl;
    KRATOS_CATCH( "" )

  }

  //**************************************************************************
  //**************************************************************************

  void AddNormalsToNodes(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    if( rModelPart.NumberOfConditions() && this->CheckConditionsDimension(rModelPart, dimension-1) ){

      ConditionsContainerType& rConditions = rModelPart.Conditions();

      //adding the normals to the nodes
      for(ConditionsContainerType::iterator it = rConditions.begin(); it !=rConditions.end(); ++it)
      {
        Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
        double coeff = 1.00/pGeometry.size();
        const array_1d<double,3>& An = it->GetValue(NORMAL);

        for(unsigned int i = 0; i<pGeometry.size(); ++i)
        {
          noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * An;
        }
      }

    }
    else if( rModelPart.NumberOfElements() && this->CheckElementsDimension(rModelPart, dimension-1) ){

      ElementsContainerType& rElements = rModelPart.Elements();

      //adding the normals to the nodes
      for(ElementsContainerType::iterator it = rElements.begin(); it !=rElements.end(); ++it)
      {
        Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
        double coeff = 1.00/pGeometry.size();
        const array_1d<double,3>& An = it->GetValue(NORMAL);

        for(unsigned int i = 0; i<pGeometry.size(); ++i)
        {
          noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * An;
        }
      }

    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  void AddWeightedNormalsToNodes(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    ModelPart::NodesContainerType& rNodes = rModelPart.Nodes();

    unsigned int MaxNodeId = MesherUtilities::GetMaxNodeId(rModelPart);
    std::vector<int> NodeNeighboursIds(MaxNodeId+1);
    std::fill( NodeNeighboursIds.begin(), NodeNeighboursIds.end(), 0 );

    if(mEchoLevel > 1)
      std::cout<<"   ["<<rModelPart.Name()<<"] [conditions:"<<rModelPart.NumberOfConditions()<<", elements:"<<rModelPart.NumberOfElements()<<"] dimension: "<<dimension<<std::endl;


    if( rModelPart.NumberOfConditions() && this->CheckConditionsLocalSpace(rModelPart, dimension-1) ){

      if(mEchoLevel > 0)
        std::cout<<"   ["<<rModelPart.Name()<<"] (C)"<<std::endl;

      //add the neighbour boundary conditions to all the nodes in the mesh
      ModelPart::ConditionsContainerType& rConditions = rModelPart.Conditions();

      std::vector<WeakPointerVector<Condition> > Neighbours(rNodes.size()+1);

      unsigned int id = 1;
      for(ModelPart::ConditionsContainerType::iterator i_cond = rConditions.begin(); i_cond!=rConditions.end(); ++i_cond)
      {
        if(i_cond->IsNot(CONTACT) && i_cond->Is(BOUNDARY)){

          Condition::GeometryType& pGeometry = i_cond->GetGeometry();

          if( mEchoLevel > 2 )
            std::cout<<" Condition ID "<<i_cond->Id()<<" id "<<id<<std::endl;

          for(unsigned int i = 0; i < pGeometry.size(); ++i)
          {
            if( mEchoLevel > 2 ){
              if(NodeNeighboursIds.size()<=pGeometry[i].Id())
                std::cout<<" Shrink node in geom "<<pGeometry[i].Id()<<" number of nodes "<<NodeNeighboursIds.size()<<std::endl;
            }

            if(NodeNeighboursIds[pGeometry[i].Id()]==0){
              NodeNeighboursIds[pGeometry[i].Id()]=id;
              Neighbours[id].push_back( Condition::WeakPointer( *(i_cond.base()) ) );
              id++;
            }
            else{
              Neighbours[NodeNeighboursIds[pGeometry[i].Id()]].push_back( Condition::WeakPointer( *(i_cond.base()) ) );
            }

          }

        }
      }


      //**********Set Boundary Nodes Only
      if( id > 1 ){
        ModelPart::NodesContainerType::iterator nodes_begin = rNodes.begin();
        ModelPart::NodesContainerType  BoundaryNodes;

        if( rModelPart.Is(FLUID) ){
          for(unsigned int i = 0; i<rNodes.size(); ++i)
          {
            if((nodes_begin + i)->Is(BOUNDARY) && (nodes_begin + i)->IsNot(RIGID) && NodeNeighboursIds[(nodes_begin+i)->Id()]!=0){
              BoundaryNodes.push_back( *((nodes_begin+i).base()) );
            }
          }
        }
        else{

          for(unsigned int i = 0; i<rNodes.size(); ++i)
          {
            if((nodes_begin + i)->Is(BOUNDARY) && NodeNeighboursIds[(nodes_begin+i)->Id()]!=0){
              BoundaryNodes.push_back( *((nodes_begin+i).base()) );
            }
          }
        }

        ComputeBoundaryShrinkage<Condition>( BoundaryNodes, Neighbours, NodeNeighboursIds, dimension);
      }

    }
    else if( rModelPart.NumberOfElements() && this->CheckElementsLocalSpace(rModelPart, dimension-1)){

      if(mEchoLevel > 0)
        std::cout<<"   ["<<rModelPart.Name()<<"] (E) "<<std::endl;

      //add the neighbour boundary elements to all the nodes in the mesh
      ModelPart::ElementsContainerType& rElements = rModelPart.Elements();

      std::vector<WeakPointerVector<Element> > Neighbours(rNodes.size()+1);

      unsigned int id = 1;
      for(ModelPart::ElementsContainerType::iterator i_elem = rElements.begin(); i_elem!=rElements.end(); ++i_elem)
      {
        if(i_elem->IsNot(CONTACT) && i_elem->Is(BOUNDARY)){

          Condition::GeometryType& pGeometry = i_elem->GetGeometry();

          if( mEchoLevel > 2 )
            std::cout<<" Element ID "<<i_elem->Id()<<" id "<<id<<std::endl;

          for(unsigned int i = 0; i < pGeometry.size(); ++i)
          {
            if( mEchoLevel > 2 ){
              if(NodeNeighboursIds.size()<=pGeometry[i].Id())
                std::cout<<" Shrink node in geom "<<pGeometry[i].Id()<<" number of nodes "<<NodeNeighboursIds.size()<<" Ids[id] "<<NodeNeighboursIds[pGeometry[i].Id()]<<std::endl;
            }

            if(NodeNeighboursIds[pGeometry[i].Id()]==0){
              NodeNeighboursIds[pGeometry[i].Id()]=id;
              Neighbours[id].push_back( Element::WeakPointer( *(i_elem.base()) ) );
              id++;
            }
            else{
              Neighbours[NodeNeighboursIds[pGeometry[i].Id()]].push_back( Element::WeakPointer( *(i_elem.base()) ) );
            }

          }

        }
      }

      //**********Set Boundary Nodes Only
      if( id > 1 ){
        ModelPart::NodesContainerType::iterator nodes_begin = rNodes.begin();
        ModelPart::NodesContainerType  BoundaryNodes;

        for(unsigned int i = 0; i<rNodes.size(); ++i)
        {
          if((nodes_begin + i)->Is(BOUNDARY) && NodeNeighboursIds[(nodes_begin+i)->Id()]!=0 ){
            BoundaryNodes.push_back( *((nodes_begin+i).base()) );
          }
        }

        ComputeBoundaryShrinkage<Element>( BoundaryNodes, Neighbours, NodeNeighboursIds, dimension );
      }

    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  template<class TClassType>
  void ComputeBoundaryShrinkage(ModelPart::NodesContainerType& rNodes, const std::vector<WeakPointerVector<TClassType> >& rNeighbours, const std::vector<int>& rNodeNeighboursIds, const unsigned int& dimension )
  {
    KRATOS_TRY

    //********** Construct the normals in order to have a good direction and length for applying a contraction
    ModelPart::NodesContainerType::iterator NodesBegin = rNodes.begin();
    int NumberOfNodes = rNodes.size();

    int not_assigned = 0;

    #pragma omp parallel for reduction(+:not_assigned)
    for(int pn=0; pn<NumberOfNodes; ++pn)
    {
      ModelPart::NodesContainerType::iterator iNode = NodesBegin + pn;
      unsigned int Id = rNodeNeighboursIds[(iNode)->Id()];

      double MeanCosinus=0;

      // are integer values that will be used as quotients
      double TotalFaces=0;
      double ProjectionValue=0;
      double NumberOfNormals=0;

      int SingleFaces=0;
      int TipAcuteAngles=0;

      array_1d<double,3>&  rNormal = (iNode)->FastGetSolutionStepValue(NORMAL);
      double&        rShrinkFactor = (iNode)->FastGetSolutionStepValue(SHRINK_FACTOR);

      //Initialize
      unsigned int NumberOfNeighbourNormals = rNeighbours[Id].size();
      if( NumberOfNeighbourNormals != 0 )
      {
        noalias(rNormal) = ZeroVector(3);
        rShrinkFactor = 0;
      }

      if( mEchoLevel > 1 ){

        std::cout<<" Id "<<Id<<" normals size "<<NumberOfNeighbourNormals<<" normal "<<rNormal<<" shrink "<<rShrinkFactor<<std::endl;
        for (unsigned int i_norm=0; i_norm<NumberOfNeighbourNormals; ++i_norm)//loop over node neighbour faces
        {
          std::cout<<" normal ["<<i_norm<<"]["<<(iNode)->Id()<<"]: "<<rNeighbours[Id][i_norm].GetValue(NORMAL)<<std::endl;
        }
      }

      Vector FaceNormals(NumberOfNeighbourNormals);
      noalias(FaceNormals) = ZeroVector(NumberOfNeighbourNormals);
      Vector TipNormals(NumberOfNeighbourNormals);
      noalias(TipNormals) = ZeroVector(NumberOfNeighbourNormals);

      //Check coincident faces-normals
      //Boundary normal is the reference normal:
      //Neighbour normals : FaceNormals[i] = 0 (not coincident and not similar normal) 1 (coincident) 2 (similar) -> other(quotient or similar normals)
      for(unsigned int i_norm=0; i_norm<NumberOfNeighbourNormals; ++i_norm)//loop over node neighbour faces
      {

        const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions-elements

        if( FaceNormals[i_norm] != 1 ){ //if is not marked as coincident

          double CloseProjection=0;
          for (unsigned int j_norm=i_norm+1; j_norm<NumberOfNeighbourNormals; ++j_norm)//loop over node neighbour faces
          {

            const array_1d<double,3>& rNormalVector = rNeighbours[Id][j_norm].GetValue(NORMAL); //conditions-elements

            ProjectionValue = inner_prod(rEntityNormal,rNormalVector);

            if( FaceNormals[j_norm]==0 ){ //in order to not have coincident normals
              //coincident normal
              if( ProjectionValue > 0.995 ){
                FaceNormals[j_norm] = 1;
              }
              //close normal
              else if( ProjectionValue > 0.955 ){
                CloseProjection    += 1;
                FaceNormals[j_norm] = 2;
                FaceNormals[i_norm] = 2;
              }
            }

            if( ProjectionValue < -0.005 ){
              TipAcuteAngles     +=1;
              TipNormals[j_norm] +=1;
              TipNormals[i_norm] +=1;
            }

          }//end for the j_norm for

          for (unsigned int j_norm=i_norm; j_norm<NumberOfNeighbourNormals; ++j_norm)//loop over node neighbour faces
          {
            if(FaceNormals[j_norm]==2 && CloseProjection>0)
              FaceNormals[j_norm]=(1.0/(CloseProjection+1));
          }

          if( FaceNormals[i_norm]!=0 ){ //it is assigned
            rNormal += rEntityNormal*FaceNormals[i_norm]; //outer normal
            NumberOfNormals += FaceNormals[i_norm];
          }
          else{
            rNormal += rEntityNormal; //outer normal
            NumberOfNormals += 1;
          }

          TotalFaces+=1;
        }

      }

      //Modify direction of tip normals in 3D  --------- start ----------

      if(dimension==3 && NumberOfNormals>=3 && TipAcuteAngles>0){

        std::vector< array_1d<double,3> > NormalsTriad(3);
        SingleFaces=0;

        if(NumberOfNormals==3){ //Deterministic solution when there are 3 no-coplanar planes

          for (unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          {
            if(TipNormals[i_norm]>=1 && FaceNormals[i_norm]==0){ //tip normal and no coincident, no similar face normal

              const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions
              NormalsTriad[SingleFaces]=rEntityNormal;
              SingleFaces+=1;
            }
          }

        }
        else if (NumberOfNormals>3) { //Reduction to the deterministic solution with 3 planes

          std::vector<int> SignificativeNormals(NumberOfNeighbourNormals);
          std::fill(SignificativeNormals.begin(), SignificativeNormals.end(), -1 );

          double MaxProjection=0;
          double Projection=0;

          //get the positive largest projection between planes
          for(unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          {
            MaxProjection=std::numeric_limits<double>::min();

            if(TipNormals[i_norm]>=1 && FaceNormals[i_norm]!=1){ //tip normal and no coincident face normal

              const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions

              for (unsigned int j_norm=0;j_norm<NumberOfNeighbourNormals;++j_norm) //loop over node neighbour faces to check the most coplanar ones
              {
                if(TipNormals[j_norm]>=1 && FaceNormals[j_norm]!=1 && i_norm!=j_norm){ //tip normal and no coincident face normal

                  const array_1d<double,3>& rNormalVector = rNeighbours[Id][j_norm].GetValue(NORMAL); //conditions
                  Projection=inner_prod(rEntityNormal,rNormalVector);

                  if(MaxProjection<Projection && Projection>0.3){ //most coplanar normals with a projection largest than 0.3
                    MaxProjection=Projection;
                    SignificativeNormals[i_norm]=j_norm; //set in SignificativeNormals
                  }

                }
              }
            }
          }

          // std::cout<<" SignificativeNormals [ ";
          // for(unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          // {
          //   std::cout <<SignificativeNormals[i_norm] <<" ";
          // }
          // std::cout<<"]"<<std::endl;
          // std::cout<<" Face: "<<FaceNormals<<" Tip: "<<TipNormals<<" Id: "<<(iNode)->Id()<<" NumberOfNormals "<<NumberOfNormals<<" SingleFaces "<<SingleFaces<<std::endl;

          //get the most obtuse planes normals
          for (unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          {

            if(SingleFaces<3){

              if(TipNormals[i_norm]>=1 && FaceNormals[i_norm]!=1){ //tip normal and no coincident face normal

                array_1d<double,3> EntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions

                if(FaceNormals[i_norm]!=0) // similar face normal
                  EntityNormal*=FaceNormals[i_norm];

                for (unsigned int j_norm=i_norm+1;j_norm<NumberOfNeighbourNormals;++j_norm)//loop over node neighbour faces to check the most obtuse planes
                {
                  if(TipNormals[j_norm]>=1 && FaceNormals[j_norm]!=1){ //tip normal and no coincident face normal

                    if(i_norm==(unsigned int)SignificativeNormals[j_norm]){ //i_norm is a significative normal with a close coincident projection

                      const array_1d<double,3>& rNormalVector = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions

                      if(FaceNormals[i_norm]!=0){ // similar face normal
                        EntityNormal+=rNormalVector*FaceNormals[i_norm]; //product with the quotient
                      }
                      else{
                        EntityNormal+=rNormalVector;
                      }

                      if(norm_2(EntityNormal)!=0)
                        EntityNormal=EntityNormal/norm_2(EntityNormal);

                      FaceNormals[j_norm]=1; //mark as coincident in order to not use this plane again

                    }
                  }
                }

                NormalsTriad[SingleFaces]=EntityNormal;
                ++SingleFaces;

              }
            }
          }

        }

        if(SingleFaces<3){

          for (unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          {
            if(SingleFaces==3)
              break;

            if(TipNormals[i_norm]>=1 && FaceNormals[i_norm]!=0 && FaceNormals[i_norm]!=1){ //tip normal and similar face normal
              const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions
              NormalsTriad[SingleFaces]=rEntityNormal;
              ++SingleFaces;
            }

            if(TipNormals[i_norm]==0 && FaceNormals[i_norm]==0 && SingleFaces>0){ //if only two of the tip normals are acute one is not set
              const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions
              NormalsTriad[SingleFaces]=rEntityNormal;
              ++SingleFaces;
            }

          }

        }

        if(SingleFaces==2){

          for (unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
          {
            if(SingleFaces==3)
              break;

            if(TipNormals[i_norm]==0 && FaceNormals[i_norm]!=0 && FaceNormals[i_norm]!=1 && SingleFaces>0){ //if only two of the tip normals are acute one is not set
              NormalsTriad[SingleFaces]=rNormal;
              ++SingleFaces;
            }
          }
        }


        if(SingleFaces==3){

          array_1d<double,3> CrossProductN;
          MathUtils<double>::CrossProduct(CrossProductN,NormalsTriad[1],NormalsTriad[2]);
          double Projection =inner_prod(NormalsTriad[0],CrossProductN);

          if(fabs(Projection) > 1e-15){
            rNormal = CrossProductN;
            MathUtils<double>::CrossProduct(CrossProductN,NormalsTriad[2],NormalsTriad[0]);
            rNormal += CrossProductN;
            MathUtils<double>::CrossProduct(CrossProductN,NormalsTriad[0],NormalsTriad[1]);
            rNormal += CrossProductN;

            rNormal /= Projection;  //intersection of three planes
          }
          else{
            rNormal = 0.5 * (NormalsTriad[1] + NormalsTriad[2]);
            // std::cout<<" Projection of the THREE planes value "<<Projection<<" Id: "<<(iNode)->Id()<<std::endl;
            // std::cout<<" Face: "<<FaceNormals<<" Tip: "<<TipNormals<<" Id: "<<(iNode)->Id()<<" NumberOfNormals "<<NumberOfNormals<<" SingleFaces "<<SingleFaces<<" Normals: "<<NormalsTriad[0]<<" "<<NormalsTriad[1]<<" "<<NormalsTriad[2]<<std::endl;
          }

        }
        else if( SingleFaces == 2 ){
          rNormal = 0.5 * (NormalsTriad[0] + NormalsTriad[1]);
          // std::cout<<" WARNING something was wrong in normal calculation Triad Not found" <<std::endl;
          // std::cout<<" Face: "<<FaceNormals<<" Tip: "<<TipNormals<<" Id: "<<(iNode)->Id()<<" NumberOfNormals "<<NumberOfNormals<<" SingleFaces "<<SingleFaces<<std::endl;
        }
        else{
          std::cout<<" WARNING something was wrong in normal calculation Triad Not found" <<std::endl;
          std::cout<<" Face: "<<FaceNormals<<" Tip: "<<TipNormals<<" Id: "<<(iNode)->Id()<<" NumberOfNormals "<<NumberOfNormals<<" SingleFaces "<<SingleFaces<<std::endl;
        }

        //Check Boundary
        if(norm_2(rNormal)>2){ //not an attractive solution but needed to ensure consistency
          rNormal =rNormal/norm_2(rNormal);
          rNormal*=1.2;
        }

        //Check inconsistencies
        // unsigned int inverted = 0;
        // for (unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neighbour faces
        // {
        //   const array_1d<double,3>& rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions
        //   if( inner_prod(rEntityNormal,rNormal) < 0 )
        //     ++inverted;
        // }

        // if( inverted == NumberOfNeighbourNormals )
        //   std::cout<<" DIRECTION NEGATIVE "<<rNormal<<" Id: "<<(iNode)->Id()<<std::endl;

        // if( norm_2(rNormal) < 1e-5 )
        //   std::cout<<" ZERO NORMAL "<<rNormal<<" Id: "<<(iNode)->Id()<<std::endl;


        //Modify direction of tip normals in 3D ----------- end  -----------

      }
      else{

        if(norm_2(rNormal)!=0)
          rNormal=rNormal/norm_2(rNormal);

        for(unsigned int i_norm=0;i_norm<NumberOfNeighbourNormals;++i_norm)//loop over node neigbour faces
        {
          if (FaceNormals[i_norm]!=1){ //not coincident normal

            array_1d<double,3> rEntityNormal = rNeighbours[Id][i_norm].GetValue(NORMAL); //conditions

            if(norm_2(rEntityNormal))
              rEntityNormal/=norm_2(rEntityNormal);

            ProjectionValue=inner_prod(rEntityNormal,rNormal);

            MeanCosinus+=ProjectionValue;
          }
        }


        if(TotalFaces!=0)
          MeanCosinus*=(1.0/TotalFaces);

        //acute angles are problematic, reduction of the modulus in that cases
        if (MeanCosinus<=0.3){
          MeanCosinus=0.8;
        }

        if(MeanCosinus>3)
          std::cout<<" WARNING WRONG MeanCosinus "<<MeanCosinus<<std::endl;

        if (MeanCosinus!=0 && MeanCosinus>1e-3) { //to ensure consistency
          MeanCosinus=1.0/MeanCosinus;
          rNormal*=MeanCosinus;               //only to put the correct length
        }
        else{
          std::cout<<" Mean Cosinus not consistent "<<MeanCosinus<<" ["<<(iNode)->Id()<<"] rNormal "<<rNormal[0]<<" "<<rNormal[1]<<" "<<rNormal[2]<<std::endl;
        }
      }

      //Now Normalize Normal and store the Shrink_Factor
      rShrinkFactor=norm_2(rNormal);

      if(rShrinkFactor!=0)
      {
        if( mEchoLevel > 1 )
          std::cout<<"[Id "<<(iNode)->Id()<<" shrink_factor "<<rShrinkFactor<<" Normal "<<rNormal[0]<<" "<<rNormal[1]<<" "<<rNormal[2]<<" MeanCosinus "<<MeanCosinus<<"] shrink "<<std::endl;
        rNormal/=rShrinkFactor;

      }
      else{

        if( mEchoLevel > 1 )
          std::cout<<"[Id "<<(iNode)->Id()<<" Normal "<<rNormal[0]<<" "<<rNormal[1]<<" "<<rNormal[2]<<" MeanCosinus "<<MeanCosinus<<"] no shrink "<<std::endl;

        noalias(rNormal) = ZeroVector(3);
        rShrinkFactor=1;

        ++not_assigned;
      }


    }


    if(mEchoLevel > 0)
      std::cout<<"   [NORMALS SHRINKAGE (BOUNDARY NODES:"<<NumberOfNodes<<") [SET:"<<NumberOfNodes-not_assigned<<" / NOT_SET:"<<not_assigned<<"] "<<std::endl;


    KRATOS_CATCH( "" )
  }


  // template<class TClassType>
  // void ComputeBoundaryShrinkage(ModelPart::NodesContainerType& rNodes, const std::vector<WeakPointerVector<TClassType> >& rNeighbours, const std::vector<int>& rIds, const unsigned int& dimension )
  // {

  //   KRATOS_TRY

  //   //********** Construct the normals in order to have a good direction and length for applying a contraction
  //   std::vector<double> storenorm;
  //   std::vector<double> tipnormal;

  //   ModelPart::NodesContainerType::iterator boundary_nodes_begin = rNodes.begin();
  //   int boundary_nodes = rNodes.size();

  //   int not_assigned = 0;

  //   int pn=0;
  //   // #pragma omp parallel for private(pn,storenorm,tipnormal)
  //   for (pn=0; pn<boundary_nodes; ++pn)
  //     {

  //       double cosmedio=0;
  //       double totalnorm=0;
  //       double directiontol=0;
  //       double numnorm=0;
  //       int indepnorm=0,acuteangle=0;

  //       array_1d<double,3>&  Normal=(boundary_nodes_begin + pn)->FastGetSolutionStepValue(NORMAL);
  //       double&       shrink_factor=(boundary_nodes_begin + pn)->FastGetSolutionStepValue(SHRINK_FACTOR);

  //       //initialize

  //       unsigned int normals_size = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]].size();
  //       if( normals_size != 0 )
  //       {
  //         Normal.clear();
  //         shrink_factor=0;
  //       }

  //       if( mEchoLevel > 1 ){
  //         std::cout<<" Id "<<rIds[(boundary_nodes_begin + pn)->Id()]<<" normals size "<<normals_size<<" normal "<<Normal<<" shrink "<<shrink_factor<<std::endl;
  //         for (unsigned int esnod=0; esnod<normals_size; ++esnod)//loop over node neighbor faces
  //         {
  //           std::cout<<" normal ["<<esnod<<"]["<<(boundary_nodes_begin + pn)->Id()<<"]: "<<rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL)<<std::endl;
  //         }
  //       }


  //       storenorm.resize(normals_size);
  //       std::fill(storenorm.begin(), storenorm.end(), 0 );

  //       tipnormal.resize(normals_size);
  //       std::fill(tipnormal.begin(), tipnormal.end(), 0 );


  //       //Check coincident faces-normals
  //       for (unsigned int esnod=0; esnod<normals_size; ++esnod)//loop over node neighbor faces
  //         {

  //           const array_1d<double,3>& AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

  //           if (storenorm[esnod]!=1){ //if is not marked as coincident

  //     	if (esnod+1<normals_size){

  //     	  double tooclose=0;
  //     	  for (unsigned int esn=esnod+1; esn<normals_size; ++esn)//loop over node neighbor faces
  //     	    {

  //     	      const array_1d<double,3>& NormalVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esn].GetValue(NORMAL); //conditions

  //     	      directiontol=inner_prod(AuxVector,NormalVector);

  //     	      //std::cout<<"  -- > direction tol"<<directiontol<<" AuxVector "<<AuxVector<<" NormalVector "<<NormalVector<<std::endl;

  //     	      if (directiontol>0.995 && storenorm[esn]==0){//to not have coincident normals
  //     		storenorm[esn]=1;
  //     	      }
  //     	      else{
  //     		if (directiontol>0.95 && directiontol<0.995 && storenorm[esn]==0){//to not have close normals
  //     		  tooclose+=1;
  //     		  storenorm[esn]  =2;
  //     		  storenorm[esnod]=2;
  //     		}
  //     		else{
  //     		  if(directiontol<-0.005){
  //     		    //std::cout<<" acute angle "<<directiontol<<std::endl;
  //     		    acuteangle      +=1;
  //     		    tipnormal[esn]  +=1;
  //     		    tipnormal[esnod]+=1;
  //     		  }
  //     		}
  //     	      }

  //     	    }//end for the esnod for

  //     	  for (unsigned int esn=esnod; esn<normals_size; ++esn)//loop over node neighbour faces
  //     	    {
  //     	      if(storenorm[esn]==2 && tooclose>0)
  //     		storenorm[esn]=(1.0/(tooclose+1));
  //     	    }


  //     	}

  //     	if(storenorm[esnod]!=0){ //!=1

  //     	  Normal+=AuxVector*storenorm[esnod]; //outer normal
  //     	  numnorm  +=storenorm[esnod];

  //     	}
  //     	else{

  //     	  Normal+=AuxVector; //outer normal
  //     	  numnorm  +=1;

  //     	}

  //     	totalnorm+=1;


  //           }

  //         }


  //       //std::cout<<" Normal "<<Normal<<" ID "<<(boundary_nodes_begin + pn)->Id()<<std::endl;

  //       //Modify direction of tip normals in 3D  --------- start ----------

  //       if(dimension==3 && numnorm>=3 && acuteangle){

  //         //std::cout<<" pn "<<pn<<" normal number: "<<numnorm<<" acute angles: "<<acuteangle<<std::endl;
  //         //std::cout<<"storenormal"<<storenorm<<std::endl;
  //         //std::cout<<"tipnormal"<<tipnormal<<std::endl;


  //         std::vector< array_1d<double,3> > N(3);
  //         indepnorm=0;

  //         if(numnorm==3){ //Definite solution for 3 planes

  //           for (unsigned int esnod=0;esnod<normals_size;++esnod)//loop over node neighbor faces
  //     	{
  //     	  if(tipnormal[esnod]>=1 && storenorm[esnod]==0){ //tip node correction

  //     	    const array_1d<double,3>& AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
  //     	    N[indepnorm]=AuxVector;
  //     	    indepnorm+=1;
  //     	  }
  //     	}

  //         }
  //         else{ //I must select only 3 planes

  //           double maxprojection=0;
  //           double projection=0;
  //           std::vector<int> pronormal(normals_size);
  //           pronormal.clear();

  //           //get the positive biggest projection between planes
  //           for (unsigned int esnod=0;esnod<normals_size;++esnod)//loop over node neighbor faces
  //     	{
  //     	  maxprojection=0;

  //     	  if(tipnormal[esnod]>=1 && storenorm[esnod]!=1){

  //     	    const array_1d<double,3>& AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

  //     	    for (unsigned int esn=0;esn<normals_size;++esn)//loop over node neighbor faces to check the most coplanar
  //     	      {
  //     		if(tipnormal[esn]>=1 && storenorm[esn]!=1 && esnod!=esn){


  //     		  const array_1d<double,3>& NormalVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esn].GetValue(NORMAL); //conditions
  //     		  projection=inner_prod(AuxVector,NormalVector);

  //     		  if(maxprojection<projection && projection>0.5){
  //     		    maxprojection=projection;
  //     		    pronormal[esnod]=esn+1; //set in pronormal one position added to the real one the most coplanar planes
  //     		  }
  //     		}
  //     	      }

  //     	  }

  //     	}

  //           // 	    std::cout<<" PROJECTIONS "<<pn<<" pronormal"<<pronormal<<std::endl;

  //           //get the most obtuse normals
  //           for (unsigned int esnod=0;esnod<normals_size;++esnod)//loop over node neighbor faces
  //     	{

  //     	  if(indepnorm<3){

  //     	    if(tipnormal[esnod]>=1 && storenorm[esnod]!=1){ //tip node correction

  //     	      array_1d<double,3> AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

  //     	      if(storenorm[esnod]!=0) //!=1
  //     		AuxVector*=storenorm[esnod];

  //     	      for (unsigned int esn=esnod+1;esn<normals_size;++esn)//loop over node neighbor faces to check the most coplanar
  //     		{
  //     		  if(tipnormal[esn]>=1 && storenorm[esn]!=1){ //tip node correction

  //     		    if(esnod+1==(unsigned int)pronormal[esn]){

  //     		      const array_1d<double,3>& NormalVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

  //     		      if(storenorm[esnod]!=0){ //!=1
  //     			AuxVector+=NormalVector*storenorm[esnod];
  //     		      }
  //     		      else{
  //     			AuxVector+=NormalVector;
  //     		      }
  //     		      AuxVector=AuxVector/norm_2(AuxVector);
  //     		      //std::cout<<" plane "<<esnod<<" esn "<<esn<<std::endl;
  //     		      storenorm[esn]=1; //to not use this plane again
  //     		    }
  //     		  }

  //     		}


  //     	      N[indepnorm]=AuxVector;
  //     	      indepnorm+=1;

  //     	    }

  //     	  }
  //     	  // 		else{

  //     	  // 		  std::cout<<" pn "<<pn<<" pronormal"<<pronormal<<std::endl;
  //     	  // 		}


  //     	}

  //         }

  //         //std::cout<<" indepnorm "<<indepnorm<<std::endl;

  //         if(indepnorm<3){

  //           for (unsigned int esnod=0;esnod<normals_size;++esnod)//loop over node neighbor faces
  //     	{
  //     	  if(indepnorm==3)
  //     	    break;

  //     	  if(tipnormal[esnod]>=1 && storenorm[esnod]!=0 && storenorm[esnod]!=1){ //tip node correction
  //     	    const array_1d<double,3>& AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
  //     	    N[indepnorm]=AuxVector;
  //     	    indepnorm+=1;
  //     	  }

  //     	  if(storenorm[esnod]==0 && tipnormal[esnod]==0 && indepnorm>0){ //if only two of the tip normals are acute one is not stored (corrected here 05/09/2011)
  //     	    const array_1d<double,3>& AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions
  //     	    N[indepnorm]=AuxVector;
  //     	    indepnorm+=1;
  //     	  }

  //     	}
  //         }



  //         if(indepnorm==3){
  //           array_1d<double,3> CrossProductN;
  //           MathUtils<double>::CrossProduct(CrossProductN,N[1],N[2]);
  //           double aux=inner_prod(N[0],CrossProductN);
  //           if(aux!=0){
  //     	MathUtils<double>::CrossProduct(CrossProductN,N[1],N[2]);
  //     	Normal = CrossProductN;
  //     	MathUtils<double>::CrossProduct(CrossProductN,N[2],N[0]);
  //     	Normal += CrossProductN;
  //     	MathUtils<double>::CrossProduct(CrossProductN,N[0],N[1]);
  //     	Normal += CrossProductN;
  //     	if( aux > 1e-15 )
  //     	  Normal /= aux;  //intersection of three planes
  //           }
  //           // else{
  //           //   std::cout<<" aux "<<aux<<" Normals :";
  //           // }

  //         }


  //         //Check Boundary
  //         if(norm_2(Normal)>2){ //a bad solution... but too ensure consistency
  //           //Normal.write(" TO LARGE NORMAL ");
  //           //std::cout<<" modulus 1 "<<Normal.modulus()<<std::endl;
  //           Normal=Normal/norm_2(Normal);
  //           Normal*=1.2;
  //           //Normal.write(" CORRRECTION ");
  //           //std::cout<<" modulus 2 "<<Normal.modulus()<<std::endl;
  //         }

  //         //Modify direction of tip normals in 3D ----------- end  -----------

  //       }
  //       else{

  //         //std::cout<<" Normal "<<Normal<<std::endl;

  //         if(norm_2(Normal)!=0)
  //           Normal=Normal/norm_2(Normal); //normalize normal *this/modulus();

  //         for(unsigned int esnod=0;esnod<normals_size;++esnod)//loop over node neigbour faces
  //           {
  //     	if (storenorm[esnod]!=1){

  //     	  array_1d<double,3> AuxVector = rNeighbours[rIds[(boundary_nodes_begin + pn)->Id()]][esnod].GetValue(NORMAL); //conditions

  //     	  if(norm_2(AuxVector))
  //     	    AuxVector/=norm_2(AuxVector);

  //     	  directiontol=inner_prod(AuxVector,Normal);

  //     	  cosmedio+=directiontol;

  //     	}
  //           }


  //         if (totalnorm!=0)
  //           cosmedio*=(1.0/totalnorm);

  //         //acute angles are problematic, reduction of the modulus in that cases
  //         if (cosmedio<=0.3){
  //           cosmedio=0.8;
  //           //std::cout<<" acute correction "<<std::endl;
  //         }

  //         if(cosmedio>3)
  //           std::cout<<" cosmedio "<<cosmedio<<std::endl;

  //         //if (cosmedio!=0) {
  //         if (cosmedio!=0 && cosmedio>1e-3) { //to ensure consistency
  //           cosmedio=1.0/cosmedio;
  //           Normal*=cosmedio;               //only to put the correct length
  //           //(boundary_nodes_begin + pn)->SetValue(Normal);
  //           //std::cout<<pn<<" cosmedio "<<cosmedio<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<std::endl;
  //         }
  //       }

  //       //std::cout<<(boundary_nodes_begin + pn)->Id()<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<std::endl;

  //       //Now Normalize Normal and store the Shrink_Factor
  //       shrink_factor=norm_2(Normal);

  //       if(shrink_factor!=0)
  //         {
  //           if( mEchoLevel > 2 )
  //     	std::cout<<"[Id "<<rIds[(boundary_nodes_begin + pn)->Id()]<<" shrink_factor "<<shrink_factor<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<" cosmedio "<<cosmedio<<"] shrink "<<std::endl;
  //           Normal/=shrink_factor;

  //         }
  //       else{

  //         if( mEchoLevel > 1 )
  //           std::cout<<"[Id "<<rIds[(boundary_nodes_begin + pn)->Id()]<<" Normal "<<Normal[0]<<" "<<Normal[1]<<" "<<Normal[2]<<" cosmedio "<<cosmedio<<"] no shrink "<<std::endl;

  //         Normal.clear();
  //         shrink_factor=1;

  //         //std::cout<<" ERROR: normal shrinkage calculation failed "<<std::endl;
  //         not_assigned +=1;
  //       }


  //     }


  //   if(mEchoLevel > 0)
  //     std::cout<<"   [NORMALS SHRINKAGE (BOUNDARY NODES:"<<boundary_nodes<<") [SET:"<<boundary_nodes-not_assigned<<" / NOT_SET:"<<not_assigned<<"] "<<std::endl;


  //   KRATOS_CATCH( "" )

  // }

  ///@}
  ///@name Private  Access
  ///@{
  ///@}
  ///@name Serialization
  ///@{
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class BoundaryNormalsCalculationUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{
///@}

}  // namespace Kratos.

#endif // KRATOS_BOUNDARY_NORMALS_CALCULATION_UTILITIES_H_INCLUDED  defined

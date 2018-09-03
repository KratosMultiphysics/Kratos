//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/EP_point_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  EPPointRigidContactPenalty2DCondition::EPPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : EPPointRigidContactPenalty3DCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  EPPointRigidContactPenalty2DCondition::EPPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : EPPointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  EPPointRigidContactPenalty2DCondition::EPPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : EPPointRigidContactPenalty3DCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  EPPointRigidContactPenalty2DCondition::EPPointRigidContactPenalty2DCondition( EPPointRigidContactPenalty2DCondition const& rOther )
  : EPPointRigidContactPenalty3DCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer EPPointRigidContactPenalty2DCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<EPPointRigidContactPenalty2DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }


   //************************************CLONE*******************************************
   //************************************************************************************

   Condition::Pointer EPPointRigidContactPenalty2DCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
   {
      EPPointRigidContactPenalty2DCondition NewCondition( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall);
      NewCondition.mCurrentInfo = this->mCurrentInfo;
      NewCondition.mSavedInfo   = this->mSavedInfo;

      // in the constructor of NewCondition I create a new friction law and here I clone the this->
      NewCondition.mpFrictionLaw = this->mpFrictionLaw->Clone();

      return Kratos::make_shared<EPPointRigidContactPenalty2DCondition>(NewCondition);
   }

  //************************************************************************************
  //************************************************************************************
   EPPointRigidContactPenalty2DCondition::~EPPointRigidContactPenalty2DCondition()
   {

   }



   //**************************** CalculateSomeSortOfArea ******************************
   //***********************************************************************************
  double EPPointRigidContactPenalty2DCondition::CalculateSomeSortOfArea()
  {

     double Area = 1;
     const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();
     if ( dimension != 2)
        return Area;


     WeakPointerVector<Element >& rNeighbourElements = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);


     std::vector< double > AreaVector;
     for ( unsigned int el = 0; el < rNeighbourElements.size() ; el++) {

        const Geometry< Node < 3 > > & rElemGeom = rNeighbourElements[el].GetGeometry();
        unsigned int nBoundary = 0;

        std::vector< unsigned int > BoundaryNodes;
        for ( unsigned int i = 0; i < rElemGeom.size(); i++) {

           if ( rElemGeom[i].Is(BOUNDARY) ) {
               array_1d<double, 3>  CN = rElemGeom[i].FastGetSolutionStepValue(CONTACT_NORMAL);

               if ( ( fabs(CN[0]) + fabs(CN[1]) ) > 0.01) {
              BoundaryNodes.push_back( i );
              nBoundary += 1;
              if ( nBoundary == 2)
                 break;
               }
           }
        }

        if ( nBoundary == 2)
        {
           array_1d< double, 3 > Vector1 = rElemGeom[ BoundaryNodes[1] ].Coordinates() - rElemGeom[ BoundaryNodes[0] ].Coordinates();
           double ThisArea = MathUtils<double>::Norm3( Vector1);
           AreaVector.push_back(ThisArea);
        }

     }

     if ( AreaVector.size() > 0) {
        Area = 0;
        for (unsigned int i = 0; i < AreaVector.size() ; i++)
           Area += AreaVector[i];
        Area /= 2.0;
     }


     return Area;
  }




  //***********************************************************************************
  //************************************************************************************

  double& EPPointRigidContactPenalty2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( this->GetProperties().Has(THICKNESS) ) {
       if ( dimension == 2 && GetProperties()[THICKNESS]>0 )
          rIntegrationWeight *= GetProperties()[THICKNESS];
    }

    return rIntegrationWeight;
  }


} // Namespace Kratos

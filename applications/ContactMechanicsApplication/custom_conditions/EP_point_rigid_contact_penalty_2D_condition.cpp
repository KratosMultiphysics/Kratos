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
    return Condition::Pointer(new EPPointRigidContactPenalty2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
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

      return Condition::Pointer( new EPPointRigidContactPenalty2DCondition( NewCondition)  ); 
   }

  //************************************************************************************
  //************************************************************************************
   EPPointRigidContactPenalty2DCondition::~EPPointRigidContactPenalty2DCondition()
   {

   }

   //************************************************************************************
   //************************************************************************************
   void EPPointRigidContactPenalty2DCondition::CalculateContactFactors(ConditionVariables &rVariables)
  {

      KRATOS_TRY

      //Compute the neighbour distance, then a stress-"like" may be computed.
      WeakPointerVector<Node<3> >& rN  = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);
      array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
      array_1d<double,3> Neighb_Point;

      double distance = 0;
      double counter = 0;

      for(unsigned int i = 0; i < rN.size(); i++)
      {
         if(rN[i].Is(BOUNDARY)){

            Neighb_Point[0] = rN[i].X();
            Neighb_Point[1] = rN[i].Y();
            Neighb_Point[2] = rN[i].Z();

            distance += norm_2(Contact_Point-Neighb_Point);

            counter ++;
  }
      }

      if( counter != 0 )
         distance /= counter;

      if( distance == 0 )
         distance = 1;           


      rVariables.ContributoryFactor = distance;

      //get contact properties and parameters
      double PenaltyParameter = 1;
      if( GetProperties().Has(PENALTY_PARAMETER) )
         PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

      WeakPointerVector<Element >& rE = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
      double ElasticModulus = 0;
      if( GetProperties().Has(YOUNG_MODULUS) )
         ElasticModulus = GetProperties()[YOUNG_MODULUS];
      else
         ElasticModulus = rE.front().GetProperties()[YOUNG_MODULUS];

      // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
      if (ElasticModulus <= 1.0e-5) {
         if ( mElasticYoungModulus < 0) {
            std::vector<double> ModulusVector;
            ProcessInfo SomeProcessInfo;
            for ( unsigned int i = 0; i < rE.size(); i++)
            {
               rE[i].CalculateOnIntegrationPoints(EQUIVALENT_YOUNG_MODULUS, ModulusVector, SomeProcessInfo);
               ElasticModulus += ModulusVector[0];
            }
            ElasticModulus /= double(rE.size());
            mElasticYoungModulus = ElasticModulus;
         } else {
            ElasticModulus = mElasticYoungModulus;
         }
      }


      rVariables.Penalty.Normal  = 1e4 * PenaltyParameter * ElasticModulus;


      double ContributoryArea = this->CalculateSomeSortOfArea();


      //std::cout << " this->Id() " << this->Id() << " , " << rVariables.Penalty.Normal << " area " << ContributoryArea << " , " << Contact_Point[0] << " , " << rVariables.Penalty.Normal * ContributoryArea << std::endl;

      rVariables.Penalty.Normal *= ContributoryArea;


      double PenaltyRatio = 1;
      if( GetProperties().Has(TANGENTIAL_PENALTY_RATIO) )
         PenaltyRatio = GetProperties()[TANGENTIAL_PENALTY_RATIO];

      rVariables.Penalty.Tangent = rVariables.Penalty.Normal * PenaltyRatio ;

      //std::cout<<" ContactPoint["<<this->Id()<<"]: penalty_n"<<rVariables.Penalty.Normal<<", ElasticModulus: "<<ElasticModulus<<", distance: "<<distance<<std::endl;

      //set contact normal
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         GetGeometry()[i].SetLock();

         array_1d<double, 3> &ContactNormal  = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_NORMAL);

         for(unsigned int i=0; i<3; i++)
            ContactNormal[i] = rVariables.Surface.Normal[i];	 

         GetGeometry()[i].UnSetLock();
      }

      KRATOS_CATCH( "" )

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

    if ( dimension == 2 && GetProperties()[THICKNESS]>0 ) 
      rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
  }


} // Namespace Kratos




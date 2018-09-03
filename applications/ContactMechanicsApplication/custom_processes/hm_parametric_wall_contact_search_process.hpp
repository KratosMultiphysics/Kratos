//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:               January 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_HM_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_HM_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED


// External includes

// System includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "custom_processes/parametric_wall_contact_search_process.hpp"

#include "custom_conditions/hydraulic_rigid_contact_penalty_3D_condition.hpp"
#include "custom_conditions/hydraulic_axisym_rigid_contact_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"

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


   ///@name Kratos Classes
   ///@{

   /// The base class for all processes in Kratos.
   /** The process is the base class for all processes and defines a simple interface for them.
     Execute method is used to execute the Process algorithms. While the parameters of this method
     can be very different from one Process to other there is no way to create enough overridden
     versions of it. For this reason this method takes no argument and all Process parameters must
     be passed at construction time. The reason is that each constructor can take different set of
     argument without any dependency to other processes or the base Process class.
    */
   class HMParametricWallContactSearchProcess
      : public ParametricWallContactSearchProcess
   {
      public:
         ///@name Type Definitions
         ///@{

         /// Pointer definition of Process
         KRATOS_CLASS_POINTER_DEFINITION( HMParametricWallContactSearchProcess );

         typedef ModelPart::NodeType                   NodeType;
         typedef ModelPart::ConditionType         ConditionType;
         typedef ModelPart::PropertiesType       PropertiesType;
         typedef ConditionType::GeometryType       GeometryType;
         typedef Point2D<ModelPart::NodeType>       Point2DType;
         typedef Point3D<ModelPart::NodeType>       Point3DType;
         //typedef FrictionLaw::pointer           FrictionLawType;

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         HMParametricWallContactSearchProcess(ModelPart& rMainModelPart): ParametricWallContactSearchProcess(rMainModelPart) {}


         HMParametricWallContactSearchProcess( ModelPart& rMainModelPart,
               std::string rSubModelPartName,
               SpatialBoundingBox::Pointer pParametricWall,
               Parameters CustomParameters)
            : ParametricWallContactSearchProcess( rMainModelPart, rSubModelPartName, pParametricWall, CustomParameters)
         {
            KRATOS_TRY

            mpConditionType->SetValue( HYDRAULIC, false);

            mpHydraulicConditionType = CreateConditionPrototypeHM( CustomParameters );
            mpHydraulicConditionType->SetValue(HYDRAULIC, true);

            KRATOS_CATCH(" ")

         }

         /// Destructor.
         virtual ~HMParametricWallContactSearchProcess() {}


         ///@}
         ///@name Operators
         ///@{

         /// This operator is provided to call the process as a function and simply calls the Execute method.
         void operator()()
         {
            Execute();
         }


         ///@}
         ///@name Operations
         ///@{


         ///@}
         ///@name Access
         ///@{


         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         virtual std::string Info() const override
         {
            return "HMParametricWallContactSearchProcess";
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "HMParametricWallContactSearchProcess";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
         }


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



         virtual void CreateContactConditions() override
         {
            KRATOS_TRY

            ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
            double Dimension = rCurrentProcessInfo[SPACE_DIMENSION];

            ModelPart::ConditionsContainerType ContactConditions;

            ModelPart& rContactModelPart = mrMainModelPart.GetSubModelPart(mContactModelPartName);

            if( mEchoLevel > 1 ){
               std::cout<<"    ["<<rContactModelPart.Name()<<" :: CONDITIONS [OLD:"<<rContactModelPart.NumberOfConditions();
            }

            unsigned int id = mrMainModelPart.Conditions().back().Id() + 1;

            ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

            // create contact condition for rigid and deformable bodies
            for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
            {
               if( (*nd)->Is(BOUNDARY) && (*nd)->Is(CONTACT) ){

                  ConditionType::Pointer pCondition;

                  if( (*nd)->Is(RIGID) ){  //rigid wall contacting with a rigid body

                     GeometryType::Pointer pGeometry;
                     if( Dimension == 2 )
                       pGeometry = Kratos::make_shared<Point2DType>(*nd);
                     else if( Dimension == 3 )
                       pGeometry = Kratos::make_shared<Point3DType>(*nd);

                     //pCondition= Kratos::make_shared<RigidBodyPointRigidContactCondition>(id, pGeometry, mpProperties, mpParametricWall);

                     ContactConditions.push_back(pCondition);

                  }
                  else{ //rigid wall contacting with a deformable body

                     Condition::NodesArrayType pConditionNode;
                     pConditionNode.push_back( (*nd) );

                     ConditionType::Pointer pConditionType = FindPointConditionHM(rContactModelPart, (*nd) , false);
                     pCondition = pConditionType->Clone(id, pConditionNode);
                     pCondition->Set(CONTACT);
                     pCondition->SetValue(HYDRAULIC, false);
                     ContactConditions.push_back(pCondition);

                     pConditionType = FindPointConditionHM(rContactModelPart, (*nd) , true);
                     pCondition = pConditionType->Clone(id+1, pConditionNode);
                     pCondition->Set(CONTACT);
                     pCondition->SetValue(HYDRAULIC, true);
                     ContactConditions.push_back(pCondition);
                  }

                  id +=2;
               }

            }


            rContactModelPart.Conditions().swap(ContactConditions);


            if( mEchoLevel > 1 ){
               std::cout<<" / NEW:"<<rContactModelPart.NumberOfConditions()<<"] "<<std::endl;
            }

            std::string ModelPartName;

            //Add contact conditions to computing domain
            for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin(); i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
            {
               if(i_mp->Is(SOLID) && i_mp->Is(ACTIVE))
                  ModelPartName = i_mp->Name();
            }

            AddContactConditions(rContactModelPart, mrMainModelPart.GetSubModelPart(ModelPartName));

            //Add contact conditions to  main domain( with AddCondition are already added )
            //AddContactConditions(rContactModelPart, mrMainModelPart);

            if( mEchoLevel >= 1 )
               std::cout<<"  [CONTACT CANDIDATES : "<<rContactModelPart.NumberOfConditions()<<"] ("<<mContactModelPartName<<") "<<std::endl;

            KRATOS_CATCH( "" )

         }


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

         ConditionType::Pointer  mpHydraulicConditionType;

         ///@}
         ///@name Private Operators
         ///@{

         //**************************************************************************
         //**************************************************************************

         ConditionType::Pointer CreateConditionPrototypeHM( Parameters& CustomParameters )
         {
            KRATOS_TRY




            ProcessInfo& rCurrentProcessInfo= mrMainModelPart.GetProcessInfo();
            double Dimension = rCurrentProcessInfo[SPACE_DIMENSION];

            // create geometry prototype for the contact conditions
            GeometryType::Pointer pGeometry;
            if( Dimension == 2 )
              pGeometry = Kratos::make_shared<Point2DType>(*((mrMainModelPart.Nodes().begin()).base()));
            else if( Dimension == 3 )
              pGeometry = Kratos::make_shared<Point3DType>(*((mrMainModelPart.Nodes().begin()).base()));


            unsigned int LastConditionId = 1;
            if( mrMainModelPart.NumberOfConditions() != 0 )
               LastConditionId = mrMainModelPart.Conditions().back().Id() + 1;

            std::string ConditionName = CustomParameters["hydraulic_condition_type"].GetString();

            if ( ConditionName == "HydraulicPointContactCondition2D1N" ) {
              return  Kratos::make_shared<HydraulicRigidContactPenalty3DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
            }
            else if ( ConditionName == "HydraulicAxisymPointContactCondition2D1N") {
              return  Kratos::make_shared<HydraulicAxisymRigidContactPenalty2DCondition>(LastConditionId, pGeometry, mpProperties, mpParametricWall);
            } else {
              std::cout << ConditionName << std::endl;
              KRATOS_ERROR << "the specified hydraulic contact condition does not exist " << std::endl;
            }

            return NULL;

            KRATOS_CATCH( "" )
         }


         ///@}
         ///@name Private Operations
         ///@{

         ConditionType::Pointer FindPointConditionHM(ModelPart & rModelPart, Node<3>::Pointer pPoint, bool rHydraulic)
         {
            KRATOS_TRY

            for(ModelPart::ConditionsContainerType::iterator i_cond =rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); ++i_cond)
            {
               if( i_cond->Is(CONTACT) && i_cond->GetGeometry().size() == 1 ){
                  if( i_cond->GetGeometry()[0].Id() == pPoint->Id() ){
                     if ( i_cond->GetValue(HYDRAULIC) == rHydraulic) {
                        return ( *(i_cond.base()) );
                     }
                  }
               }
            }

            if ( rHydraulic == false) {
               return mpConditionType;
            }
            return mpHydraulicConditionType;

            KRATOS_CATCH("")
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
         HMParametricWallContactSearchProcess& operator=(HMParametricWallContactSearchProcess const& rOther);

         /// Copy constructor.
         //Process(Process const& rOther);


         ///@}

   }; // Class Process

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   /// input stream function
   inline std::istream& operator >> (std::istream& rIStream,
         HMParametricWallContactSearchProcess& rThis);

   /// output stream function
   inline std::ostream& operator << (std::ostream& rOStream,
         const HMParametricWallContactSearchProcess& rThis)
   {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
   }
   ///@}


}  // namespace Kratos.

#endif // KRATOS_HM_PARAMETRIC_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined

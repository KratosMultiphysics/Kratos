//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_RIGID_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED


// External includes

// System includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.hpp"
#include "custom_conditions/axisym_point_rigid_contact_penalty_water_2D_condition.hpp"
#include "custom_conditions/beam_point_rigid_contact_penalty_3D_condition.hpp"
#include "custom_conditions/beam_point_rigid_contact_LM_3D_condition.hpp"
#include "custom_conditions/rigid_body_point_rigid_contact_condition.hpp"

#include "pfem_solid_mechanics_application.h"


namespace Kratos
{

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
   class RigidWallContactSearchProcess
      : public Process
   {
      public:
         ///@name Type Definitions
         ///@{

         /// Pointer definition of Process
         KRATOS_CLASS_POINTER_DEFINITION( RigidWallContactSearchProcess );

         typedef ModelPart::ConditionType         ConditionType;
         typedef ModelPart::PropertiesType       PropertiesType;
         typedef ConditionType::GeometryType       GeometryType;
         typedef Point2D<ModelPart::NodeType>       Point2DType;
         typedef Point3D<ModelPart::NodeType>       Point3DType;
         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         RigidWallContactSearchProcess(ModelPart& rModelPart): mrModelPart(rModelPart) {}


         RigidWallContactSearchProcess(SpatialBoundingBox::Pointer pRigidWall, ModelPart& rModelPart, int EchoLevel, bool waterCondition) 
            : mrModelPart(rModelPart)
         {
            mpRigidWall = pRigidWall;
            mEchoLevel = EchoLevel;
            mWaterCondition = waterCondition;
         } 

         /// Destructor.
         virtual ~RigidWallContactSearchProcess() {}    ///@{

         /// This operator is provided to call the process as a function and simply calls the Execute method.
         void operator()()
         {
            Execute();
         }


         ///@}
         ///@name Operations
         ///@{


         /// Execute method is used to execute the Process algorithms.
         virtual void Execute() {}

         /// this function is designed for being called at the beginning of the computations
         /// right after reading the model and the groups
         virtual void ExecuteInitialize()
         {
         }

         /// this function is designed for being execute once before the solution loop but after all of the
         /// solvers where built
         virtual void ExecuteBeforeSolutionLoop()
         {
         }


         /// this function will be executed at every time step BEFORE performing the solve phase
         virtual void ExecuteInitializeSolutionStep()
         {

            KRATOS_TRY

            std::cout << " FINALLY I GOT IT " << mWaterCondition <<  std::endl;
            ProcessInfo& CurrentProcessInfo= mrModelPart.GetProcessInfo();	  
            double Time      = CurrentProcessInfo[TIME];   

            //update center position
            mpRigidWall->UpdatePosition( Time );

            if (Time == 0)
               KRATOS_THROW_ERROR( std::logic_error, "detected time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )

                  ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

            //Create Rigid Contact Conditions
            int id = mrModelPart.Conditions().back().Id() + 1;

            mConditionsNumber =  mrModelPart.Conditions().size();

            //Check if elements are beams or two nodes and set BOUNDARY flag to nodes
            for(ModelPart::ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie!=mrModelPart.ElementsEnd(); ie++){

               if( ie->GetGeometry().size() == 2 ){
                  for( unsigned int i=0; i<2; i++ )
                  {
                     ie->GetGeometry()[i].Set(BOUNDARY,true);
                  }
               }

            }


            for(ModelPart::NodesContainerType::iterator nd = rNodes.begin(); nd != rNodes.end(); ++nd)
            {
               nd->Set(CONTACT,false);
            }


            //Check ModelPart meshes for Rigid Domains and set BOUNDARY flag to nodes
            ModelPart::MeshesContainerType& rMeshes = mrModelPart.GetMeshes();
            //bool RigidBodyPresent = false;
            for(unsigned int m=0; m<rMeshes.size(); m++){

               if( rMeshes[m].Is(RIGID) ){

                  //RigidBodyPresent = true;

                  for(ModelPart::ElementsContainerType::iterator ie = rMeshes[m].ElementsBegin(); ie!=rMeshes[m].ElementsEnd(); ie++){

                     for( unsigned int i=0; i<ie->GetGeometry().size(); i++ )
                     {
                        ie->GetGeometry()[i].Set(BOUNDARY,true);
                     }
                  }
               }

            }


#ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
#else
            int number_of_threads = 1;
#endif


            vector<unsigned int> nodes_partition;
            OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), nodes_partition);

            Vector WallVelocity =  mpRigidWall->Velocity();
            int  MovementLabel  =  mpRigidWall->GetMovementLabel();

#pragma omp parallel
            {
               int k = OpenMPUtils::ThisThread();
               ModelPart::NodesContainerType::iterator NodesBegin = rNodes.begin() + nodes_partition[k];
               ModelPart::NodesContainerType::iterator NodesEnd = rNodes.begin() + nodes_partition[k + 1];

               for(ModelPart::NodesContainerType::const_iterator nd = NodesBegin; nd != NodesEnd; nd++)
               {

                  //set point rigid wall condition : usually in non rigid_wall points
                  if( nd->SolutionStepsDataHas(RIGID_WALL) ){

                     if( nd->GetSolutionStepValue(RIGID_WALL) == MovementLabel ){

                        //nd->Set(STRUCTURE);
                        nd->Set(RIGID);

                        //set new coordinates (MOVE MESH)
                        nd->X() = nd->X0() + WallVelocity[0] * Time;
                        nd->Y() = nd->Y0() + WallVelocity[1] * Time;
                        nd->Z() = nd->Z0() + WallVelocity[2] * Time;

                        //std::cout<<" node "<<(nd)->Id()<<" Position ("<<(nd)->X()<<", "<<(nd)->Y()<<" "<<(nd)->Z()<<") "<<std::endl;
                     }
                  }

               }

               for(ModelPart::NodesContainerType::const_iterator nd = NodesBegin; nd != NodesEnd; nd++)
               {
                  if( nd->Is(BOUNDARY) ){

                     if( nd->IsNot(RIGID) )
                     {
                        //to perform contact with a tube radius must be set
                        double Radius = 0;
                        //if( !RigidBodyPresent ){
                        if( nd->IsNot(SLAVE) ){
                           Radius = nd->GetValue(MEAN_RADIUS);
                        }
                        else{
                           Radius = 0;
                        }

                        Vector Point(3);
                        Point[0] = nd->X();
                        Point[1] = nd->Y();
                        Point[2] = nd->Z();

                        if( mpRigidWall->IsInside(Point,Time,Radius) ){
                           nd->Set(CONTACT);
                        }
                        else{ //clean nodal contact forces
                           nd->FastGetSolutionStepValue(CONTACT_FORCE).clear();
                        }
                     }

                     }

                  }
               }


               // **************** Serial version of the same search:
               // Vector WallVelocity =  mpRigidWall->Velocity();
               // int  MovementLabel  =  mpRigidWall->GetMovementLabel();

               // //Check RIGID walls and search contacts
               // for(ModelPart::NodesContainerType::iterator nd = rNodes.begin(); nd != rNodes.end(); ++nd)
               // 	{
               // 	  //set point rigid wall condition : usually in non rigid_wall points
               // 	  if( nd->SolutionStepsDataHas(RIGID_WALL) ){

               // 	    if( nd->GetSolutionStepValue(RIGID_WALL) == MovementLabel ){

               // 	      //nd->Set(STRUCTURE);
               // 	      nd->Set(RIGID);

               // 	      //set new coordinates (MOVE MESH)
               // 	      nd->X() = nd->X0() + WallVelocity[0] * Time;
               // 	      nd->Y() = nd->Y0() + WallVelocity[1] * Time;
               // 	      nd->Z() = nd->Z0() + WallVelocity[2] * Time;

               // 	      //std::cout<<" node "<<(nd)->Id()<<" Position ("<<(nd)->X()<<", "<<(nd)->Y()<<" "<<(nd)->Z()<<") "<<std::endl;
               // 	    }
               // 	  }

               // 	  if( nd->Is(BOUNDARY) ){

               // 	    if( nd->IsNot(RIGID) )
               // 	      {
               // 		//to perform contact with a tube radius must be set
               // 		double Radius = 0;
               // 		if( !RigidBodyPresent ){
               // 		  Radius = nd->GetValue(MEAN_RADIUS);
               // 		}
               // 		else{
               // 		  Radius = 0;
               // 		}

               // 		Vector Point(3);
               // 		Point[0] = nd->X();
               // 		Point[1] = nd->Y();
               // 		Point[2] = nd->Z();

               // 		if( mpRigidWall->IsInside(Point,Time,Radius) ){
               // 		  nd->Set(CONTACT);
               // 		}
               // 	      }

               // 	  }

               // 	}

               //*********************

               //Check RIGID walls and search contacts
               for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
               {

                  // if( (*nd)->Is(BOUNDARY) )
                  //std::cout<<" Node "<<(*nd)->Id()<<" Is boundary "<<std::endl;

                  //if( (*nd)->Is(BOUNDARY) && (*nd)->IsNot(RIGID) ){
                  if( (*nd)->Is(BOUNDARY) ){

                     //if( (*nd)->IsNot(RIGID) && !RigidBodyPresent ){//rigid wall contacting with a deformable body 
                     if( (*nd)->IsNot(RIGID) && (*nd)->IsNot(SLAVE) ){//rigid wall contacting with a deformable body 

                        if( (*nd)->Is(CONTACT) ){

                           //std::cout<<" Node Selected "<<(*nd)->Id()<<std::endl;
                           //contact parameters in properties
                           int number_properties = mrModelPart.NumberOfProperties();
                           PropertiesType::Pointer p_properties = mrModelPart.pGetProperties(number_properties-1);

                           ConditionType::Pointer p_cond;

                           if( mpRigidWall->GetDimension() == 2 ){

                              GeometryType::Pointer p_geometry = GeometryType::Pointer(new Point2DType( (*nd) ));
                              if( mpRigidWall->Axisymmetric() == true ){
                                 if ( mWaterCondition)  {
                                    p_cond= ModelPart::ConditionType::Pointer(new AxisymPointRigidContactPenaltyWater2DCondition(id, p_geometry, p_properties, mpRigidWall) );
                                 }
                                 else {
                                    p_cond= ModelPart::ConditionType::Pointer(new AxisymPointRigidContactPenalty2DCondition(id, p_geometry, p_properties, mpRigidWall) );
                                 }
                                 //std::cout<<": Set Contact 2D axisymmetric condition "<<std::endl;
                              }
                              else{
                                 p_cond= ModelPart::ConditionType::Pointer(new PointRigidContactPenalty2DCondition(id, p_geometry, p_properties, mpRigidWall) );  
                                 //std::cout <<": Set Contact 2D condition "<<std::endl;
                              }

                           }
                           else if( mpRigidWall->GetDimension() == 3 ){

                              //Set Beam Element Properties
                              WeakPointerVector<Element>& rE = (*nd)->GetValue(NEIGHBOUR_ELEMENTS);

                              PropertiesType& Properties = rE[0].GetProperties();
                              double Radius = 0;
                              double Max_Radius = Properties[MEAN_RADIUS];

                              int element_id = 0;
                              for(unsigned int ie=1; ie<rE.size(); ie++)
                              {
                                 PropertiesType& Properties = rE[ie].GetProperties();      
                                 Radius = Properties[MEAN_RADIUS];
                                 if( Max_Radius < Radius )
                                    element_id =ie;
                              }


                              double PenaltyParameter = (*p_properties)[PENALTY_PARAMETER];
                              p_properties = rE[element_id].pGetProperties();      
                              p_properties->SetValue( PENALTY_PARAMETER, PenaltyParameter );

                              //std::cout<<" BEAM radius considered for contact "<<(*p_properties)[MEAN_RADIUS]<<std::endl;

                              GeometryType::Pointer p_geometry = GeometryType::Pointer(new Point3DType( (*nd) ));
                              //p_cond= ModelPart::ConditionType::Pointer(new PointRigidContactPenalty3DCondition(id, p_geometry, p_properties, mpRigidWall) ); 
                              p_cond= ModelPart::ConditionType::Pointer(new BeamPointRigidContactPenalty3DCondition(id, p_geometry, p_properties, mpRigidWall) ); 	       
                              //p_cond= ModelPart::ConditionType::Pointer(new BeamPointRigidContactLM3DCondition(id, p_geometry, p_properties, mpRigidWall) ); 

                              //std::cout<<" Node Selected for BEAM Contact "<<(*nd)->Id()<<": Set Contact 3D condition "<<std::endl;
                              //std::cout<<" with properties "<<*p_properties<<std::endl;
                           }

                           //pcond->SetValue(mpRigidWall); the boundingbox of the rigid wall must be passed to the condition

                           mrModelPart.Conditions().push_back(p_cond);

                           id +=1;

                        }
                     }
                     else{ //rigid wall contacting with a rigid body 


                        if( (*nd)->IsNot(RIGID) ){


                           if( (*nd)->Is(CONTACT) ){

                              //contact parameters in properties
                              int number_properties = mrModelPart.NumberOfProperties();
                              PropertiesType::Pointer p_properties = mrModelPart.pGetProperties(number_properties-1);

                              ConditionType::Pointer p_cond;


                              if( mpRigidWall->GetDimension() == 2 ){

                                 //rigid wall contacting with a 2D rigid body
                                 GeometryType::Pointer p_geometry = GeometryType::Pointer(new Point2DType( (*nd) ));
                                 std::cout<<" ERROR: POINT CONTACT CONDITION NOT BUILD "<<std::endl;
                                 //std::cout<<" Node Selected "<<(*nd)->Id()<<": Set Contact 2D condition "<<std::endl;

                              }
                              else if( mpRigidWall->GetDimension() == 3 ){

                                 //rigid wall contacting with a 3D rigid body
                                 GeometryType::Pointer p_geometry = GeometryType::Pointer(new Point3DType( (*nd) ));
                                 p_cond= ModelPart::ConditionType::Pointer(new RigidBodyPointRigidContactCondition(id, p_geometry, p_properties, mpRigidWall) ); 
                                 //std::cout<<" Node Selected "<<(*nd)->Id()<<": Set Contact 3D condition "<<std::endl;

                              }

                              mrModelPart.Conditions().push_back(p_cond);

                              id +=1;
                           }

                        }

                     }
                  }
                  }


                  if( mEchoLevel > 1 )
                     std::cout<<"  [ Rigid Contacts : "<<mrModelPart.Conditions().size() - mConditionsNumber<<" ]"<<std::endl;

                  KRATOS_CATCH( "" )

               }

               /// this function will be executed at every time step AFTER performing the solve phase
               virtual void ExecuteFinalizeSolutionStep()
               {
                  KRATOS_TRY	  

                  // To write correct displacements
                  ProcessInfo& CurrentProcessInfo= mrModelPart.GetProcessInfo();
                  double Time      = CurrentProcessInfo[TIME];
                  //double DeltaTime = CurrentProcessInfo[DELTA_TIME];

                  ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

                  int counter = 0;

                  for ( ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
                  {
                     if((*nd)->Is(RIGID) && (*nd)->FastGetSolutionStepValue(RIGID_WALL) == mpRigidWall->GetMovementLabel() ){
                        array_1d<double, 3 >& CurrentDisplacement  = (*nd)->FastGetSolutionStepValue(DISPLACEMENT);	    
                        CurrentDisplacement[0] = mpRigidWall->Velocity()[0] * Time;
                        CurrentDisplacement[1] = mpRigidWall->Velocity()[1] * Time;
                        CurrentDisplacement[2] = mpRigidWall->Velocity()[2] * Time;

                        (*nd)->Coordinates() = (*nd)->GetInitialPosition() + CurrentDisplacement;

                        counter++;
                     }
                  }

                  if( mEchoLevel > 1 )
                     std::cout<<"  [ Finalize Wall Contact : wall velocity:"<< mpRigidWall->Velocity()<<" wall nodes : "<<counter<<std::endl;

                  //Clean Rigid Contact Conditions
                  ModelPart::ConditionsContainerType NonRigidContactConditions;

                  unsigned int id=0;

                  //std::cout<<" [ NUMBER OF CONDITIONS before rigid contact update: "<<mrModelPart.Conditions().size()<<" ]"<<std::endl;

                  for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ic++)
                  {
                     if( id == mConditionsNumber )
                        break;

                     NonRigidContactConditions.push_back(*(ic.base()));  

                     id +=1;
                  }

                  mrModelPart.Conditions().swap( NonRigidContactConditions );

                  //std::cout<<"  [ NUMBER OF CONDITIONS after  rigid contact update: "<<mrModelPart.Conditions().size()<<" ]"<<std::endl;

                  //calculate elemental contribution
                  KRATOS_CATCH( "" )      
               }


               /// this function will be executed at every time step BEFORE  writing the output
               virtual void ExecuteBeforeOutputStep()
               {
               }


               /// this function will be executed at every time step AFTER writing the output
               virtual void ExecuteAfterOutputStep()
               {
               }


               /// this function is designed for being called at the end of the computations
               /// right after reading the model and the groups
               virtual void ExecuteFinalize()
               {
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

               /// Turn back information as a string.
               virtual std::string Info() const
               {
                  return "RigidWallContactSearchProcess";
               }

               /// Print information about this object.
               virtual void PrintInfo(std::ostream& rOStream) const
               {
                  rOStream << "RigidWallContactSearchProcess";
               }

               /// Print object's data.
               virtual void PrintData(std::ostream& rOStream) const
               {
               }


               ///@}
               ///@name Friends
               ///@{


               ///@}


               private:
               ///@name Static Member Variables
               ///@{

               ///@}
               ///@name Static Member Variables
               ///@{
               ModelPart&  mrModelPart;

               SpatialBoundingBox::Pointer mpRigidWall;

               unsigned int mConditionsNumber;

               int mEchoLevel;

               bool mWaterCondition;
               ///@}
               ///@name Un accessible methods
               ///@{

               /// Assignment operator.
               RigidWallContactSearchProcess& operator=(RigidWallContactSearchProcess const& rOther);

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
                  RigidWallContactSearchProcess& rThis);

            /// output stream function
            inline std::ostream& operator << (std::ostream& rOStream,
                  const RigidWallContactSearchProcess& rThis)
            {
               rThis.PrintInfo(rOStream);
               rOStream << std::endl;
               rThis.PrintData(rOStream);

               return rOStream;
            }
            ///@}


   }  // namespace Kratos.

#endif // KRATOS_RIGID_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined 



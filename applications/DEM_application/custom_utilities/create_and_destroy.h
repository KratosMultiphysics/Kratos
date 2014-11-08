//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_CREATE_AND_DESTROY )
#define  KRATOS_CREATE_AND_DESTROY

// /* External includes */

// System includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"

// Project includes
#include "includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/spheric_particle.h"
#include "../custom_elements/spheric_continuum_particle.h"
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../DEM_application_variables.h"


namespace Kratos {

    static double rand_normal(double mean, double stddev, double max_radius, double min_radius) {
        
        if(!stddev) return mean;
        
        double return_value;
        
        do {
        
            double x, y, r;

            do {
                x = 2.0 * rand() / RAND_MAX - 1;
                y = 2.0 * rand() / RAND_MAX - 1;
                r = x * x + y*y;
            } while (r == 0.0 || r > 1.0);

            double d = sqrt(-2.0 * log(r) / r);
            return_value = x*d * stddev + mean;
            
        } while (return_value < min_radius || return_value > max_radius);
        
        return return_value;
    }        
    
    
class ParticleCreatorDestructor
{
public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;
        unsigned int mMaxNodeId;
		
    KRATOS_CLASS_POINTER_DEFINITION(ParticleCreatorDestructor);

    /// Default constructor. 
    ParticleCreatorDestructor(): mGreatestParticleId(0)
    {
      mScaleFactor        = 1.0;
      mHighPoint[0]       = 10e18;
      mHighPoint[1]       = 10e18;
      mHighPoint[2]       = 10e18;
      mLowPoint[0]        = -10e18;
      mLowPoint[1]        = -10e18;
      mLowPoint[2]        = -10e18;
    }

    //Particle_Creator_Destructor() {};

    /// Destructor.

    virtual ~ParticleCreatorDestructor() {};

    
    int FindMaxNodeIdInModelPart(ModelPart& r_modelpart){
        
        int max_Id = 1; //GID accepts Id's >= 1
        //TODO: should we parallelize this?
        for (ModelPart::NodesContainerType::iterator node_it = r_modelpart.GetCommunicator().LocalMesh().NodesBegin(); node_it != r_modelpart.GetCommunicator().LocalMesh().NodesEnd(); node_it++){
	  if( (int)(node_it->Id()) > max_Id ) max_Id = node_it->Id();
	}              
        r_modelpart.GetCommunicator().MaxAll( max_Id );
        return max_Id;
    }
    
    void NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart, Node < 3 > ::Pointer& pnew_node, int aId, Node < 3 > ::Pointer & reference_node ,double radius, Properties& params, bool has_sphericity, bool has_rotation, bool initial)
    {
        
      double bx= reference_node->X();
      double cy= reference_node->Y();
      double dz= reference_node->Z();
      
      if (initial){
          pnew_node = reference_node;
          r_modelpart.AddNode(pnew_node);   // The same node is added to r_modelpart (the calculation model part)
          pnew_node->SetId(aId);
          pnew_node->FastGetSolutionStepValue(VELOCITY)[0]= 0.0;
          pnew_node->FastGetSolutionStepValue(VELOCITY)[1]= 0.0;
          pnew_node->FastGetSolutionStepValue(VELOCITY)[2]= 0.0;
          //(actually it should be the velocity of the inlet layer, which is different from the particles being inserted)
      }

      else {
          pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz);      //ACTUAL node creation and addition to model part
          pnew_node->FastGetSolutionStepValue(VELOCITY)                     = params[VELOCITY];
      }
            
      if(has_rotation){
          pnew_node->FastGetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO) = params[PARTICLE_ROTATION_DAMP_RATIO];
      }
      
      if(has_sphericity){
        pnew_node->FastGetSolutionStepValue(PARTICLE_SPHERICITY)        = params[PARTICLE_SPHERICITY];
      }      
      
      pnew_node->FastGetSolutionStepValue(RADIUS)                       = radius;
      array_1d<double, 3 > null_vector(3,0.0);                              
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY)             = null_vector;
      pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL)            = params[PARTICLE_MATERIAL];
     
      if (pnew_node->SolutionStepsDataHas(SOLID_FRACTION_PROJECTED)){
          pnew_node->GetSolutionStepValue(SOLID_FRACTION_PROJECTED)           = 0.0;
      }

      if (pnew_node->SolutionStepsDataHas(MESH_VELOCITY1)){
          pnew_node->GetSolutionStepValue(MESH_VELOCITY1_X)         = 0.0;
          pnew_node->GetSolutionStepValue(MESH_VELOCITY1_Y)         = 0.0;
          pnew_node->GetSolutionStepValue(MESH_VELOCITY1_Z)         = 0.0;
      }

      if (pnew_node->SolutionStepsDataHas(DRAG_REACTION)){
          pnew_node->GetSolutionStepValue(DRAG_REACTION_X)          = 0.0;
          pnew_node->GetSolutionStepValue(DRAG_REACTION_Y)          = 0.0;
          pnew_node->GetSolutionStepValue(DRAG_REACTION_Z)          = 0.0;
      }
      
      pnew_node->AddDof(VELOCITY_X,REACTION_X);
      pnew_node->AddDof(VELOCITY_Y,REACTION_Y);
      pnew_node->AddDof(VELOCITY_Z,REACTION_Z);
      pnew_node->AddDof(ANGULAR_VELOCITY_X,REACTION_X);
      pnew_node->AddDof(ANGULAR_VELOCITY_Y,REACTION_Y);
      pnew_node->AddDof(ANGULAR_VELOCITY_Z,REACTION_Z);
      
      pnew_node->pGetDof(VELOCITY_X)->FixDof();
      pnew_node->pGetDof(VELOCITY_Y)->FixDof();
      pnew_node->pGetDof(VELOCITY_Z)->FixDof();           
      pnew_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
      pnew_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
      pnew_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();
            
      pnew_node->Set(DEMFlags::FIXED_VEL_X,true);
      pnew_node->Set(DEMFlags::FIXED_VEL_Y,true);
      pnew_node->Set(DEMFlags::FIXED_VEL_Z,true);
      pnew_node->Set(DEMFlags::FIXED_ANG_VEL_X,true);
      pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Y,true);
      pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Z,true);

    }   

    void ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart, int r_Elem_Id, Node < 3 > ::Pointer reference_node, 
                                              Properties::Pointer r_params, const Element& r_reference_element, PropertiesProxy* p_fast_properties, bool has_sphericity, bool has_rotation, bool initial) {          
        
      Node < 3 > ::Pointer pnew_node;
      
      double radius = (*r_params)[RADIUS];   
      double max_radius = 1.5 * radius;
      if(initial){
          radius = max_radius;          
      } else {
        double std_deviation = (*r_params)[STANDARD_DEVIATION]; 

        
        double min_radius = 0.5 * radius;      

        radius = rand_normal(radius, std_deviation, max_radius, min_radius);                      
      }
            
      NodeCreatorWithPhysicalParameters(r_modelpart, pnew_node, r_Elem_Id, reference_node, radius, *r_params, has_sphericity, has_rotation, initial); 
      
      Geometry< Node < 3 > >::PointsArrayType nodelist;
      
      nodelist.push_back(pnew_node);                                              
      
      Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);

      //p_particle->GetProperties() = r_params;
      
      p_particle->Set(NEW_ENTITY);            
      pnew_node->Set(NEW_ENTITY);
      
      if(initial) {
          p_particle->Set(BLOCKED);
          pnew_node->Set(BLOCKED);  
      }
      
      p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());      
      
      Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*>(p_particle.get());       
          
      spheric_p_particle->SetFastProperties(p_fast_properties); 
            
      double density = spheric_p_particle->GetDensity();
      spheric_p_particle->SetRadius(radius);      
      double mass = 4.0 / 3.0 * KRATOS_M_PI * density * radius * radius * radius;
      spheric_p_particle->SetSqrtOfRealMass(sqrt(mass)); 
      
      if (has_rotation) spheric_p_particle->Set(DEMFlags::HAS_ROTATION,true);
      else              spheric_p_particle->Set(DEMFlags::HAS_ROTATION,false);
      
      p_particle->Initialize(); //////////////// STANDARD INITIALIZATION!!
      
      r_modelpart.Elements().push_back(p_particle);          
}    
    
    
    void NodeCreatorForClusters(ModelPart& r_modelpart, 
                                Node < 3 > ::Pointer& pnew_node,
                                int aId,
                                array_1d<double, 3 >& reference_coordinates,                                  
                                double radius, 
                                Properties& params)
    {

        double bx= reference_coordinates[0];
        double cy= reference_coordinates[1];
        double dz= reference_coordinates[2];


        pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz);      //ACTUAL node creation and addition to model part

        pnew_node->FastGetSolutionStepValue(RADIUS)                       = radius;
        array_1d<double, 3 > null_vector(3,0.0);
        pnew_node->FastGetSolutionStepValue(VELOCITY)                     = null_vector;
        pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY)             = null_vector;
        pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL)            = params[PARTICLE_MATERIAL];

        pnew_node->AddDof(VELOCITY_X,REACTION_X);
        pnew_node->AddDof(VELOCITY_Y,REACTION_Y);
        pnew_node->AddDof(VELOCITY_Z,REACTION_Z);
        pnew_node->AddDof(ANGULAR_VELOCITY_X,REACTION_X);
        pnew_node->AddDof(ANGULAR_VELOCITY_Y,REACTION_Y);
        pnew_node->AddDof(ANGULAR_VELOCITY_Z,REACTION_Z);

        pnew_node->pGetDof(VELOCITY_X)->FixDof();
        pnew_node->pGetDof(VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(VELOCITY_Z)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();

        pnew_node->Set(DEMFlags::FIXED_VEL_X,true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Y,true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Z,true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_X,true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Y,true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Z,true);
        pnew_node->Set(DEMFlags::BELONGS_TO_A_CLUSTER,true);

    }           
    
    Kratos::SphericParticle* ElementCreatorForClusters( ModelPart& r_modelpart, 
                                    int r_Elem_Id, 
                                    double radius,
                                    array_1d<double, 3 >& reference_coordinates, 
                                    double sqrt_of_cluster_mass,
                                    Properties::Pointer r_params, 
                                    const Element& r_reference_element,
                                    const int cluster_id) {          

        Node<3> ::Pointer pnew_node;

        NodeCreatorForClusters(r_modelpart, pnew_node, r_Elem_Id, reference_coordinates, radius, *r_params);

        Geometry< Node <3> >::PointsArrayType nodelist;
        nodelist.push_back(pnew_node);

        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);
        p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());

        Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*>(p_particle.get());

        spheric_p_particle->SetRadius(radius);
        spheric_p_particle->SetSqrtOfRealMass(sqrt_of_cluster_mass);

        spheric_p_particle->Set(DEMFlags::HAS_ROLLING_FRICTION,false);
        spheric_p_particle->Set(DEMFlags::BELONGS_TO_A_CLUSTER,true);
        spheric_p_particle->SetClusterId(cluster_id);
        
        r_modelpart.Elements().push_back(p_particle);
        return spheric_p_particle;
    }

    void CalculateSurroundingBoundingBox(ModelPart& r_balls_model_part, ModelPart& r_clusters_model_part, ModelPart& r_rigid_faces_model_part, double scale_factor, bool automatic)
    {
        KRATOS_TRY

        if (automatic){
            double ref_radius = 0.0;

            if (r_balls_model_part.NumberOfElements(0) == 0 && r_clusters_model_part.NumberOfElements(0) == 0 && r_rigid_faces_model_part.NumberOfElements(0) == 0){
                KRATOS_ERROR(std::logic_error, "The Bounding Box cannot be calculated automatically when there are no elements. Kratos stops." , "");
            }

            if (scale_factor < 0.0){
                KRATOS_ERROR(std::logic_error, "The enlargement factor for the automatic calculation of the bounding box must be a positive value." , "");
            }

            if (scale_factor < 1.0){
                std::cout << "\n WARNING" + std::string( 2, '\n' );
                std::cout << "The enlargement factor for the automatic calculation of the bounding box is less than 1!. \n";
                std::cout << "It is not guaranteed that the model fits inside of it." + std::string( 2, '\n' );
            }

            if (r_balls_model_part.NumberOfElements(0)){ // loop over spheric elements (balls)
                Configure::ElementsContainerType Elements = r_balls_model_part.GetCommunicator().LocalMesh().Elements();

                ref_radius               = (*(Elements.begin().base()))->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                array_1d<double, 3> coor = (*(Elements.begin().base()))->GetGeometry()(0)->Coordinates();

                mStrictLowPoint          = coor;
                mStrictHighPoint         = coor;

                for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin(); particle_pointer_it != Elements.end(); ++particle_pointer_it){
                    coor = (*(particle_pointer_it.base()))->GetGeometry()(0)->Coordinates();
                    double radius = particle_pointer_it->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                    ref_radius = (ref_radius < radius) ? radius : ref_radius;

                    for (std::size_t i = 0; i < 3; i++){
                        mStrictLowPoint[i]  = (mStrictLowPoint[i]  > coor[i]) ? coor[i] : mStrictLowPoint[i];
                        mStrictHighPoint[i] = (mStrictHighPoint[i] < coor[i]) ? coor[i] : mStrictHighPoint[i];
                    }
                }
            }

            if (r_clusters_model_part.NumberOfElements(0)){} // loop over clusters

            if (r_rigid_faces_model_part.NumberOfConditions(0)){ // loop over rigid faces
                ModelPart::ConditionsContainerType Conditions = r_rigid_faces_model_part.GetCommunicator().LocalMesh().Conditions();

                std::vector <array_1d<double, 3> > face_coor;
                std::size_t n_nodes = (*(Conditions.begin().base()))->GetGeometry().size();
                face_coor.resize(n_nodes);

                for (std::size_t i = 0; i < n_nodes; ++i){
                    face_coor[i] = (*(Conditions.begin().base()))->GetGeometry()(i)->Coordinates();
                }

                if (r_balls_model_part.NumberOfElements(0) == 0){ // initialize if not initialized already
                    mStrictLowPoint  = face_coor[0];
                    mStrictHighPoint = face_coor[0];
                }

                for (ModelPart::ConditionsContainerType::iterator particle_pointer_it = Conditions.begin(); particle_pointer_it != Conditions.end(); ++particle_pointer_it){
                    std::size_t n_nodes = (*(particle_pointer_it.base()))->GetGeometry().size();
                    face_coor.resize(n_nodes);

                    for (std::size_t i = 0; i < n_nodes; ++i){
                        face_coor[i] = (*(particle_pointer_it.base()))->GetGeometry()(i)->Coordinates();

                        for (std::size_t j = 0; j < 3; j++){
                            mStrictLowPoint[j]  = (mStrictLowPoint[j]  > face_coor[i][j]) ? face_coor[i][j] : mStrictLowPoint[j];
                            mStrictHighPoint[j] = (mStrictHighPoint[j] < face_coor[i][j]) ? face_coor[i][j] : mStrictHighPoint[j];
                        }
                    }
                }
            }

            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[0]);
            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[1]);
            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[2]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[0]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[1]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[2]);
            r_balls_model_part.GetCommunicator().MaxAll(ref_radius);

            array_1d<double, 3> midpoint = 0.5 * (mStrictHighPoint + mStrictLowPoint);
            mHighPoint                   = midpoint * (1 - scale_factor) + scale_factor * mStrictHighPoint;
            mLowPoint                    = midpoint * (1 - scale_factor) + scale_factor * mStrictLowPoint;

            for (std::size_t i = 0; i < 3; i++){
                mLowPoint[i]  -= 2 * ref_radius;
                mHighPoint[i] += 2 * ref_radius;
            }

            mStrictDiameter  = norm_2(mStrictHighPoint - mStrictLowPoint);
            mDiameter        = norm_2(mHighPoint - mLowPoint);
        } // if (automatic)

        else {
            for (int i = 0; i < 3; ++i){

                if (mHighPoint[i] < mLowPoint[i]){
                    KRATOS_ERROR(std::logic_error,  "Check limits of the Bounding Box, minimum coordinates exceed maximum coordinates." , "");
                }
            }

            mStrictHighPoint = mHighPoint; // mHighPoint and mLowPoint have been set as an input value
            mStrictLowPoint  = mLowPoint;
            mStrictDiameter  = norm_2(mStrictHighPoint - mStrictLowPoint);
            mDiameter        = norm_2(mHighPoint - mLowPoint);
        }
               
        KRATOS_CATCH("")
    }

    void DestroyParticles(ModelPart& r_model_part)
    {

        KRATOS_TRY

        //Type definitions
        typedef ModelPart::ElementsContainerType                          ElementsArrayType;
        typedef ElementsArrayType::iterator                               ElementsIterator;
        
        ElementsArrayType& rElements        = r_model_part.GetCommunicator().LocalMesh().Elements();                       
        ModelPart::NodesContainerType& rNodes       = r_model_part.GetCommunicator().LocalMesh().Nodes();

        //ElementsArrayType& rElements                          = r_model_part.Elements();
        //ModelPart::NodesContainerType& rNodes               = r_model_part.Nodes();

        ElementsArrayType temp_particles_container;
        ModelPart::NodesContainerType temp_nodes_container;

        temp_particles_container.reserve(rElements.size());
        temp_nodes_container.reserve(rNodes.size());
        
        temp_particles_container.swap(rElements);
        temp_nodes_container.swap(rNodes);

        //Add the ones not marked with TO_ERASE
        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = temp_particles_container.ptr_begin(); particle_pointer_it != temp_particles_container.ptr_end(); ++particle_pointer_it){	  
            
            if( !(*particle_pointer_it)->GetGeometry()(0)->Is(TO_ERASE) && !(*particle_pointer_it)->Is(TO_ERASE)) {
                (rElements).push_back(*particle_pointer_it); //adding the elements

                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++){ //GENERAL FOR ELEMENTS OF MORE THAN ONE NODE
                    ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                    (rNodes).push_back(pNode);
                }
            }

        }                           
        KRATOS_CATCH("")
    }
    
    
        void DestroyContactElements(ModelPart& r_model_part)
    {
        KRATOS_TRY

        //Type definitions
        typedef ModelPart::ElementsContainerType                          ElementsArrayType;
        typedef ElementsArrayType::iterator                               ElementsIterator;

        ElementsArrayType& rElements                          = r_model_part.GetCommunicator().LocalMesh().Elements();

        ElementsArrayType temp_elements_container;

        //Copy the elements and clear the element container
        //temp_elements_container.reserve(pElements->size());
        temp_elements_container.reserve(rElements.size());        
        temp_elements_container.swap(rElements);

        //Add the ones not marked with TO_ERASE
        for (Configure::ElementsContainerType::ptr_iterator element_pointer_it = temp_elements_container.ptr_begin(); element_pointer_it != temp_elements_container.ptr_end(); ++element_pointer_it){	              

            if( !(*element_pointer_it)->Is(TO_ERASE) ) {
                (rElements).push_back(*element_pointer_it); //adding the elements               
	    }	  
        }                           
        KRATOS_CATCH("")
    }


    void MarkDistantParticlesForErasing(ModelPart& r_model_part)
    {

      MarkParticlesForErasingGivenBoundingBox(r_model_part, mLowPoint, mHighPoint);

    }

    void MarkParticlesForErasingGivenScalarVariableValue(ModelPart& r_model_part, const Variable<double>& rVariable, double value, double tol)
    {

      KRATOS_TRY

      Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){

          const double& i_value     = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(rVariable);
          bool include=true;             // = (erase_flag < 0.5);

          include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

          if(include) 
              (*particle_pointer_it)->GetGeometry()(0)->Set(TO_ERASE);

      }

      KRATOS_CATCH("")

    }

    void MarkParticlesForErasingGivenVectorVariableModulus(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& rVariable, double value, double tol)
    {

      KRATOS_TRY

      Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){

          array_1d<double, 3 >& i_var = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(rVariable);
          double i_value              = sqrt(i_var[0] * i_var[0] + i_var[1] * i_var[1] + i_var[2] * i_var[2]);
          bool include=true;              //  = (erase_flag < 0.5);

          include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

          if(include) 
              (*particle_pointer_it)->GetGeometry()(0)->Set(TO_ERASE);

      }

      KRATOS_CATCH("")

    }

    void MarkParticlesForErasingGivenBoundingBox(ModelPart& r_model_part, array_1d<double, 3> low_point, array_1d<double, 3> high_point)
    {

      KRATOS_TRY
              
      ModelPart::NodesContainerType& rNodes       = r_model_part.GetCommunicator().LocalMesh().Nodes();
      Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
      
      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){
          
          if( (*particle_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER) ) continue;

          array_1d<double, 3 > coor = (*particle_pointer_it)->GetGeometry()(0)->Coordinates();
          bool include=true;            
          
          for (unsigned int i = 0; i < 3; i++){
              include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
          }
          
          if(!include) {
              (*particle_pointer_it)->GetGeometry()(0)->Set(TO_ERASE); 
              (*particle_pointer_it)->Set(TO_ERASE); 
          }

      }
      
      for (ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin();
              node_pointer_it != rNodes.ptr_end(); ++node_pointer_it){

          if( (*node_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER) ) continue;
          array_1d<double, 3 > coor = (*node_pointer_it)->Coordinates();
          bool include=true;            
          
          for (unsigned int i = 0; i < 3; i++){
              include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
          }
          
          if(!include) {
              (*node_pointer_it)->Set(TO_ERASE);
          }
          
      }
      
      KRATOS_CATCH("")
    }
    
    
    void MarkContactElementsForErasing(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {
        ///ONLY FOR CONTINUUM!!!

      KRATOS_TRY
                                                           
      Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){
          
          if( (*particle_pointer_it)->GetGeometry()(0)->Is(TO_ERASE) )   {

              Element* p_element = particle_pointer_it->get();
              SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*>( p_element );    
              std::vector<Particle_Contact_Element*>& node_to_neigh_element_pointer = p_continuum_spheric_particle->mBondElements;
              unsigned int number_of_contact_elements = node_to_neigh_element_pointer.size();
              for(unsigned int i = 0; i < number_of_contact_elements; i++)
              {
                  Particle_Contact_Element* p_to_contact_element = node_to_neigh_element_pointer[i];
                  p_to_contact_element->Set(TO_ERASE);
              }
                            
          }       

      }

      KRATOS_CATCH("")
    }
    
    void DestroyParticlesOutsideBoundingBox(ModelPart& r_model_part)
    {
       MarkDistantParticlesForErasing(r_model_part);
       DestroyParticles( r_model_part);
    }
    
    void DestroyContactElementsOutsideBoundingBox(ModelPart& r_model_part, ModelPart& mcontacts_model_part)
    {
        MarkContactElementsForErasing(r_model_part, mcontacts_model_part);
        DestroyContactElements(mcontacts_model_part);
    }
    
    

    ///@}
    ///@name Access
    ///@{

    array_1d<double, 3> GetHighNode()
    {
        return (mHighPoint);
    }

    array_1d<double, 3> GetLowNode()
    {
        return (mLowPoint);
    }

    array_1d<double, 3> GetStrictHighNode()
    {
        return (mHighPoint);
    }

    array_1d<double, 3> GetStrictLowNode()
    {
        return (mLowPoint);
    }

    double GetDiameter()
    {
        return (mDiameter);
    }

    double GetStrictDiameter()
    {
        return (mStrictDiameter);
    }

    void SetHighNode(array_1d<double, 3> node)
    {
        mHighPoint = node;
    }

    void SetLowNode(array_1d<double, 3> node)
    {
        mLowPoint = node;
    }
    
    unsigned int GetCurrentMaxNodeId()
    {
        return mMaxNodeId;
    }
    
    void SetMaxNodeId(unsigned int id)
    {
        mMaxNodeId = id;
    }
    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.

    virtual std::string Info() const
    {
        return "";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member rVariables
    ///@{


    ///@}
    ///@name Protected member rVariables
    ///@{ template<class T, std::size_t dim>


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

    ///@name Static Member rVariables
    ///@{


    ///@}
    ///@name Member rVariables
    ///@{
    array_1d<double, 3 > mHighPoint;
    array_1d<double, 3 > mLowPoint;
    array_1d<double, 3 > mStrictHighPoint;
    array_1d<double, 3 > mStrictLowPoint;
    double mDiameter;
    double mStrictDiameter;
    double mScaleFactor;
    int mGreatestParticleId;

    inline void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size)
    {
        unsigned int buffer_size = node_it->GetBufferSize();
        for (unsigned int step = 0; step < buffer_size; step++)
        {
            //getting the data of the solution step
            double* step_data = (node_it)->SolutionStepData().Data(step);
            //copying this data in the position of the vector we are interested in
            for (int j = 0; j < step_data_size; j++)
            {
                step_data[j] = 0.0;
            }
        }
    }

    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable)
    {
        array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable);
        noalias(Aux_var) = ZeroVector(3);
    }

    inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable)
    {
        /* ///WARNING M: aixo activar-ho també
        double& Aux_var = (*particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
        Aux_var = 0.0;
         */
    }

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


    /// Assignment operator.
    ParticleCreatorDestructor & operator=(ParticleCreatorDestructor const& rOther);

    ///@}

}; // Class ParticleCreatorDestructor

class DEM_Inlet
{
public:        
    
    typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;
    typedef WeakPointerVector<Element> ParticleWeakVectorType; 
    typedef ModelPart::ElementsContainerType                          ElementsArrayType;
    
    
    Vector PartialParticleToInsert; //array of doubles, must be resized in the constructor to the number of meshes
    ModelPart& InletModelPart; //The model part used to insert elements
    bool mFirstTime;
    boost::numeric::ublas::vector<bool> mLayerRemoved;
    bool mBallsModelPartHasSphericity;
    bool mBallsModelPartHasRotation;
    std::vector<PropertiesProxy>                 mFastProperties;
    
    /// Constructor:               
    
    DEM_Inlet(ModelPart& inlet_modelpart): InletModelPart(inlet_modelpart)
    {                
        
        PartialParticleToInsert.resize(inlet_modelpart.NumberOfMeshes(),false);
        mLayerRemoved.resize(inlet_modelpart.NumberOfMeshes(),false);
        
        int mesh_iterator_number=0;   
        
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin();
                                               mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
        {
         PartialParticleToInsert[mesh_iterator_number]  = 0.0;
         mLayerRemoved[mesh_iterator_number]  = false;
         
         mesh_iterator_number++;
        }                
        
        mFirstTime=true;
        mBallsModelPartHasSphericity=false;  
        mBallsModelPartHasRotation  =false;
        
        
    }
            
    /// Destructor.
    virtual ~DEM_Inlet() {};
        
    void InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const std::string& ElementNameString){
        
        unsigned int& max_Id=creator.mMaxNodeId;       
        
        CreatePropertiesProxies(mFastProperties, InletModelPart);      
               
        VariablesList r_modelpart_nodal_variables_list = r_modelpart.GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTICLE_SPHERICITY) )  mBallsModelPartHasSphericity = true;        
        if(r_modelpart.GetProcessInfo()[ROTATION_OPTION] ) {
            mBallsModelPartHasRotation   = true;
            InletModelPart.GetProcessInfo()[ROTATION_OPTION] = true;            
        }
        else {
            InletModelPart.GetProcessInfo()[ROTATION_OPTION] = false;             
        }
        
        const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
        
        int mesh_number=0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = InletModelPart.GetMeshes().begin()+1;
                                               mesh_it != InletModelPart.GetMeshes().end();    ++mesh_it)
        {            
            mesh_number++;
            int mesh_size=mesh_it->NumberOfNodes();
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
            
            int general_properties_id = InletModelPart.GetProperties(mesh_number).Id();  
            PropertiesProxy* p_fast_properties = NULL;
            
            for (unsigned int i = 0; i < mFastProperties.size(); i++){
                int fast_properties_id = mFastProperties[i].GetId(); 
                if( fast_properties_id == general_properties_id ){  
                    p_fast_properties = &(mFastProperties[i]);
                    break;
                }
            }   
            
            for (int i = 0; i < mesh_size; i++){                
                creator.ElementCreatorWithPhysicalParameters(r_modelpart, max_Id+1, all_nodes[i], InletModelPart.pGetProperties(mesh_number), r_reference_element, p_fast_properties, mBallsModelPartHasSphericity, mBallsModelPartHasRotation,true); 
		max_Id++;
            }      
            
        } //for mesh_it                                               
    } //InitializeDEM_Inlet
        
    void DettachElements(ModelPart& r_modelpart, unsigned int& max_Id){    
                            
        vector<unsigned int> ElementPartition;
        OpenMPUtils::CreatePartition( OpenMPUtils::GetNumThreads(), r_modelpart.GetCommunicator().LocalMesh().Elements().size(), ElementPartition);
        typedef ElementsArrayType::iterator ElementIterator;       

        #pragma omp parallel for
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            ElementIterator it_begin = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k];
            ElementIterator it_end   = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator elem_it = it_begin; elem_it != it_end; ++elem_it){
         //   for (ElementsArrayType::iterator elem_it = r_modelpart.ElementsBegin(); elem_it != r_modelpart.ElementsEnd(); ++elem_it){
                            
                if(elem_it->IsNot(NEW_ENTITY)) continue;

                Kratos::SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*elem_it);

                Node < 3 > ::Pointer& node_it = elem_it->GetGeometry()(0);                  

                bool still_touching=false;

                for (unsigned int i = 0; i < spheric_particle.mNeighbourElements.size(); i++) {
                    SphericParticle* neighbour_iterator = spheric_particle.mNeighbourElements[i];
                    Node < 3 > ::Pointer& neighbour_node = neighbour_iterator->GetGeometry()(0); 
                    if( (node_it->IsNot(BLOCKED) && neighbour_node->Is(BLOCKED)) || (neighbour_node->IsNot(BLOCKED) && node_it->Is(BLOCKED)) ) {
                        still_touching=true;
                        break;
                    }              
                }

                if(!still_touching){ 
                  if(node_it->IsNot(BLOCKED)){//The ball must be freed
                    node_it->Set(DEMFlags::FIXED_VEL_X,false);
                    node_it->Set(DEMFlags::FIXED_VEL_Y,false);
                    node_it->Set(DEMFlags::FIXED_VEL_Z,false);
                    node_it->Set(DEMFlags::FIXED_ANG_VEL_X,false);
                    node_it->Set(DEMFlags::FIXED_ANG_VEL_Y,false);
                    node_it->Set(DEMFlags::FIXED_ANG_VEL_Z,false);
                    elem_it->Set(NEW_ENTITY,0);
                    node_it->Set(NEW_ENTITY,0);
                  }
                  else{
                  //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted.             
                    node_it->Set(ACTIVE,false);
                    elem_it->Set(ACTIVE,false);
                  }
                }

            } //loop nodes
        } //loop threads
        
        mFirstTime=false;
        
    } //DettachElements
    
    void CreateElementsFromInletMesh( ModelPart& r_modelpart, ModelPart& inlet_modelpart, ParticleCreatorDestructor& creator, const std::string& ElementNameString){
                        
        //unsigned int max_Id=0; 
        unsigned int& max_Id=creator.mMaxNodeId; 
        
        DettachElements(r_modelpart, max_Id);                
                
        int mesh_number=0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin()+1;
                                               mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
        {            
            mesh_number++;

            if(r_modelpart.GetProcessInfo()[TIME] <  InletModelPart.GetProperties(mesh_number)[INLET_START_TIME] ) continue;
            
            int mesh_size=mesh_it->NumberOfNodes();
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
                        
            if(r_modelpart.GetProcessInfo()[TIME] >  InletModelPart.GetProperties(mesh_number)[INLET_STOP_TIME]  ) {
                if (mLayerRemoved[mesh_number]) continue;
                for (int i = 0; i < mesh_size; i++){                   
		     all_nodes[i]->Set(TO_ERASE);		   
                }
                mLayerRemoved[mesh_number] = true;
                continue;
            }                        
                        
            double num_part_surface_time=  InletModelPart.GetProperties(mesh_number)[INLET_NUMBER_OF_PARTICLES]; 
            double delta_t=  r_modelpart.GetProcessInfo()[DELTA_TIME]; // FLUID DELTA_T CAN BE USED ALSO, it will depend on how often we call this function
            double surface=  1.0;//inlet_surface; // this should probably be projected to velocity vector
            
            //calculate number of particles to insert from input data
            double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + PartialParticleToInsert[mesh_number-1];            
            int number_of_particles_to_insert = floor(double_number_of_particles_to_insert);
            PartialParticleToInsert[mesh_number-1] = double_number_of_particles_to_insert - number_of_particles_to_insert;
            
            if (number_of_particles_to_insert) {
              //randomizing mesh
               srand( time(NULL)*r_modelpart.GetProcessInfo()[TIME_STEPS] );
               ModelPart::NodesContainerType::ContainerType inserting_nodes(number_of_particles_to_insert);               
               ModelPart::NodesContainerType::ContainerType valid_nodes = mesh_it->NodesArray();
               int valid_nodes_length=0;
               
               for (int i = 0; i < mesh_size; i++){
                   if( all_nodes[i]->IsNot(ACTIVE) ) { 
		     valid_nodes[valid_nodes_length]=all_nodes[i];   
		     valid_nodes_length++; 
		   } // (push_back) //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
               }

               if (valid_nodes_length < number_of_particles_to_insert) {
                   number_of_particles_to_insert = valid_nodes_length;
                   //std::cout<<"The number of DEM particles has been reduced to match the available number of nodes of the DEM Inlet mesh"<<std::endl<<std::flush;
               }
               
               for (int i = 0; i < number_of_particles_to_insert; i++) {
		   int pos = rand() % valid_nodes_length;
		   //int pos = i;
                   inserting_nodes[i] = valid_nodes[pos]; //This only works for pos as real position in the vector if 
                   //we use ModelPart::NodesContainerType::ContainerType 
                   //instead of ModelPart::NodesContainerType
                   valid_nodes[pos] = valid_nodes[valid_nodes_length - 1];
                   valid_nodes_length = valid_nodes_length - 1;
               }
               
                PropertiesProxy* p_fast_properties = NULL;
                int general_properties_id = InletModelPart.GetProperties(mesh_number).Id();  
                for (unsigned int i = 0; i < mFastProperties.size(); i++){
                    int fast_properties_id = mFastProperties[i].GetId(); 
                    if( fast_properties_id == general_properties_id ){  
                        p_fast_properties = &(mFastProperties[i]);
                        break;
                    }
                }     
                                             
               const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
               
               for (int i = 0; i < number_of_particles_to_insert; i++) {
                   creator.ElementCreatorWithPhysicalParameters(r_modelpart, max_Id+1, inserting_nodes[i], InletModelPart.pGetProperties(mesh_number), r_reference_element, p_fast_properties, mBallsModelPartHasSphericity, mBallsModelPartHasRotation, false);
                   inserting_nodes[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                   max_Id++;
               }                                                                   
               
           } //if (number_of_particles_to_insert)
            
        } // for mesh_it
    }    //CreateElementsFromInletMesh       
            
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{




/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}


} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined



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

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"

// Project includes
#include "includes/define.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/spheric_particle.h"
#include "../custom_elements/spheric_swimming_particle.h"
#
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../DEM_application.h"

#
//SALVA_ENDING
//const double prox_tol = 0.00000000001;
namespace Kratos
{
class ParticleCreatorDestructor
{
public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;
		


    KRATOS_CLASS_POINTER_DEFINITION(ParticleCreatorDestructor);

    /// Default constructor. 
    ParticleCreatorDestructor()
    {
      mScaleFactor  = 1.0;
      mHighPoint[0] = 10e18;
      mHighPoint[1] = 10e18;
      mHighPoint[2] = 10e18;
      mLowPoint[0]  = -10e18;
      mLowPoint[1]  = -10e18;
      mLowPoint[2]  = -10e18;
    }

    //Particle_Creator_Destructor() {};

    /// Destructor.

    virtual ~ParticleCreatorDestructor() {};

    
    void NodeCreator(ModelPart& r_modelpart, Node < 3 > ::Pointer& pnew_node, int aId, double bx, double cy, double dz) {
              
      pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz, 0.0);
      pnew_node->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(RADIUS) = 0.05;
      pnew_node->FastGetSolutionStepValue(PARTICLE_DENSITY) = 1000;
      pnew_node->FastGetSolutionStepValue(YOUNG_MODULUS) = 10000;
      pnew_node->FastGetSolutionStepValue(POISSON_RATIO) = 0.5;
      pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = 1;
          
    }
    void NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart, Node < 3 > ::Pointer& pnew_node, int aId, double bx, double cy, double dz ) {
              
      /*pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz, 0.0);
      pnew_node->FastGetSolutionStepValue(VELOCITY_X) = velX;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Y) = velY;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Z) = velZ;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(RADIUS) = radius;
      pnew_node->FastGetSolutionStepValue(PARTICLE_DENSITY) = density;
      pnew_node->FastGetSolutionStepValue(YOUNG_MODULUS) = 10000;
      pnew_node->FastGetSolutionStepValue(POISSON_RATIO) = 0.25;
      pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = color;*/
          
    }


    void ElementCreator(ModelPart& r_modelpart, int r_Elem_Id, int r_Id, double r_x, double r_y, double r_z)
    {

      Node < 3 > ::Pointer pnew_node;
      
      NodeCreator(r_modelpart, pnew_node, r_Id, r_x, r_y, r_z);
      
      Geometry< Node < 3 > >::PointsArrayType nodelist;
      
      nodelist.push_back(pnew_node);

      Element::Pointer p_swimming_particle = Element::Pointer(new SphericSwimmingParticle(r_Elem_Id, nodelist)); 
      
      r_modelpart.Elements().push_back(p_swimming_particle); 
           
    }

    void PrintingTest() {}

    //SALVA
    //MA
    void ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart, int r_Elem_Id, Node < 3 > ::Pointer reference_node, Properties& params) {          
        
      Node < 3 > ::Pointer pnew_node;
      
      NodeCreatorWithPhysicalParameters(r_modelpart, pnew_node, r_Elem_Id, reference_node->X(), reference_node->Y(), reference_node->Z() //THE ID IS WRONG
                                        /*radius, sphericity, density,  
                                        rest_coeff, fric_angle,
                                        roll_fric, rot_damp_ratio, color,
                                        velX,  velY, velZ  */); 
      
      Geometry< Node < 3 > >::PointsArrayType nodelist;
      
      nodelist.push_back(pnew_node);                                        

      Element::Pointer p_swimming_particle = Element::Pointer(new SphericSwimmingParticle(r_Elem_Id, nodelist)); 
      
      p_swimming_particle->Set(NEW_ENTITY);
      
      pnew_node->Set(NEW_ENTITY);
      
      r_modelpart.Elements().push_back(p_swimming_particle);          
}    
    //MA

    void CalculateSurroundingBoundingBox( ModelPart& r_model_part, double scale_factor)
    {
        KRATOS_TRY

        //Type definitions
        Configure::ElementsContainerType::Pointer pElements = r_model_part.pElements();
        Configure::ElementsContainerType Elements           = r_model_part.Elements();
       
        double ref_radius         = (*(Elements.begin().base()))->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        array_1d<double, 3 > coor = (*(Elements.begin().base()))->GetGeometry()(0)->Coordinates();
        mLowPoint                 = coor;
        mHighPoint                = coor;
        mStrictLowPoint           = coor;
        mStrictHighPoint          = coor;

        for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin(); particle_pointer_it != Elements.end(); ++particle_pointer_it){
            coor = (*(particle_pointer_it.base()))->GetGeometry()(0)->Coordinates();

            for (std::size_t i = 0; i < 3; i++){
                mStrictLowPoint[i]  = (mStrictLowPoint[i] > coor[i]) ? coor[i] : mStrictLowPoint[i];
                mStrictHighPoint[i] = (mStrictHighPoint[i] < coor[i]) ? coor[i] : mStrictHighPoint[i];
            }

        }

        array_1d<double, 3 > midpoint = 0.5 * (mStrictHighPoint + mStrictLowPoint);
        mHighPoint                    = midpoint * (1 - scale_factor) + scale_factor * mStrictHighPoint;
        mLowPoint                     = midpoint * (1 - scale_factor) + scale_factor * mStrictLowPoint;

        for (std::size_t i = 0; i < 3; i++){
            mLowPoint[i]  -= 2 * ref_radius;
            mHighPoint[i] += 2 * ref_radius;
        }

        mStrictDiameter = norm_2(mStrictHighPoint - mStrictLowPoint);
        mDiameter       = norm_2(mHighPoint - mLowPoint);
        
        KRATOS_CATCH("")
         
    }

    void DestroyParticles(ModelPart& r_model_part)
    {

        KRATOS_TRY

        //Type definitions
        Configure::ElementsContainerType::Pointer pElements = r_model_part.pElements();
        ModelPart::NodesContainerType::Pointer pNodes       = r_model_part.pNodes();

        Configure::ElementsContainerType& rElements         = r_model_part.Elements();
        ModelPart::NodesContainerType& rNodes               = r_model_part.Nodes();

        Configure::ElementsContainerType temp_particles_container;
        ModelPart::NodesContainerType temp_nodes_container;

        //Copy the elements and clear the element container
        temp_particles_container.reserve(pElements->size());
        temp_nodes_container.reserve(pNodes->size());
        
        temp_particles_container.swap(rElements);
        temp_nodes_container.swap(rNodes);

        //Add the ones inside the bounding box
        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = temp_particles_container.ptr_begin();
                particle_pointer_it != temp_particles_container.ptr_end(); ++particle_pointer_it){

            bool erase_flag = (0.5 < ((*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(ERASE_FLAG)));

            if (!erase_flag){
               (rElements).push_back(*particle_pointer_it); //adding the elements
               
               for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++){ //GENERAL FOR ELEMENTS OF MORE THAN ONE NODE
                   ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                   (rNodes).push_back(pNode);
               }

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

      Configure::ElementsContainerType& rElements = r_model_part.Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){

          const double& i_value     = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(rVariable);
          double& erase_flag        = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(ERASE_FLAG);
          bool include              = (erase_flag < 0.5);

          include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

          erase_flag = include ? 0.0 : 1.0;

      }

      KRATOS_CATCH("")

    }

    void MarkParticlesForErasingGivenVectorVariableModulus(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& rVariable, double value, double tol)
    {

      KRATOS_TRY

      Configure::ElementsContainerType& rElements = r_model_part.Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){

          array_1d<double, 3 >& i_var = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(rVariable);
          double i_value              = sqrt(i_var[0] * i_var[0] + i_var[1] * i_var[1] + i_var[2] * i_var[2]);
          double& erase_flag          = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(ERASE_FLAG);
          bool include                = (erase_flag < 0.5);

          include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

          erase_flag = include ? 0.0 : 1.0;

      }

      KRATOS_CATCH("")

    }

    void MarkParticlesForErasingGivenBoundingBox(ModelPart& r_model_part, array_1d<double, 3> low_point, array_1d<double, 3> high_point)
    {

      KRATOS_TRY

      Configure::ElementsContainerType& rElements = r_model_part.Elements();

      for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
              particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it){

          double& erase_flag        = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(ERASE_FLAG);
          array_1d<double, 3 > coor = (*particle_pointer_it)->GetGeometry()(0)->Coordinates();
          bool include              = (erase_flag < 0.5);

          for (unsigned int i = 0; i < 3; i++){
              include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
          }

          erase_flag = include ? 0.0 : 1.0;

      }

      KRATOS_CATCH("")
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

    array_1d<double, 3>& GetStrictHighNode()
    {
        return (mHighPoint);
    }

    array_1d<double, 3>& GetStrictLowNode()
    {
        return (mLowPoint);
    }

    double& GetDiameter()
    {
        return (mDiameter);
    }

    double& GetStrictDiameter()
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
        array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
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
  //MA
class DEM_Inlet
{
public:        
    
    Vector PartialParticleToInsert; //array of doubles, must be resized in the constructor to the number of meshes
    ModelPart& InletModelPart; //The model part used to insert elements
    
    /// Constructor:               
    
    DEM_Inlet(ModelPart& inlet_modelpart): InletModelPart(inlet_modelpart)
    {                
        
        PartialParticleToInsert.resize(inlet_modelpart.NumberOfMeshes(),false);
        
        int mesh_iterator_number=0;   
        
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin();
                                               mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
        {
         PartialParticleToInsert[mesh_iterator_number]  = 0.0;
         
         mesh_iterator_number++;
        }                
        
    }
            
    /// Destructor.
    virtual ~DEM_Inlet() {};
    
    void InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator){
        int r_Elem_Id=100; ///This Id must be set!!!
        int mesh_number=0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = InletModelPart.GetMeshes().begin()+1;
                                               mesh_it != InletModelPart.GetMeshes().end();    ++mesh_it)
        {
            mesh_number++;
            int mesh_size=mesh_it->NumberOfNodes();                                                                      
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
            //for (ModelPart::NodesContainerType::iterator node_it = all_nodes.begin();
                                                        //node_it != all_nodes.end(); node_it++){
            //for (ModelPart::NodesContainerType::iterator node_it = mesh_it->NodesBegin(); node_it != mesh_it->NodesEnd(); node_it++){ 
            for (int i = 0; i < mesh_size; i++){                
                creator.ElementCreatorWithPhysicalParameters(r_modelpart, r_Elem_Id, all_nodes[i],InletModelPart.GetProperties(mesh_number));
                
            }
        }
                                               
    }
        
    
    void CreateElementsFromInletMesh( ModelPart& r_modelpart, ModelPart& inlet_modelpart, ParticleCreatorDestructor& creator ){
        
  //      if(r_modelpart.GetProcessInfo()[TIME]< inlet_start  || r_modelpart.GetProcessInfo()[TIME]< inlet_start) return;
        
        int mesh_iterator_number=0;
        int mesh_number=0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin()+1;
                                               mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
        {
            mesh_number++;
            int mesh_size=mesh_it->NumberOfNodes();
            
            double num_part_surface_time=  1; //num_parts; ////////////////////////////////
            double delta_t=  r_modelpart.GetProcessInfo()[DELTA_TIME]; // FLUID DELTA_T CAN BE USED ALSO
            double surface=  1.0;//inlet_surface; // this should probably be projected to velocity vector
            
            //calculate number of particles to insert from input data
            double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + PartialParticleToInsert[mesh_iterator_number];            
            int number_of_particles_to_insert = floor(double_number_of_particles_to_insert);
            PartialParticleToInsert[mesh_iterator_number] = double_number_of_particles_to_insert - number_of_particles_to_insert;

            if (number_of_particles_to_insert) {
              //randomizing mesh
               srand(time(NULL));
               ModelPart::NodesContainerType::ContainerType inserting_nodes(number_of_particles_to_insert);
               ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
               ModelPart::NodesContainerType::ContainerType valid_nodes = mesh_it->NodesArray();
               int valid_nodes_length=0;
               
               for (int i = 0; i < mesh_size; i++){
                   if( all_nodes[i]->IsNot(ACTIVE) ) { valid_nodes[valid_nodes_length]=all_nodes[i];   valid_nodes_length++;  }
               }

               if (valid_nodes_length < number_of_particles_to_insert) {
                   number_of_particles_to_insert = valid_nodes_length;
                   std::cout<<"The number of DEM particles has been reduced to match the available number of nodes of the DEM Inlet mesh"<<std::endl<<std::flush;
               }
               
               for (int i = 0; i < number_of_particles_to_insert; i++) {
                   int pos = rand() % valid_nodes_length;
                   inserting_nodes[i] = valid_nodes[pos]; //This only works for pos as real position in the vector if 
                   //we use ModelPart::NodesContainerType::ContainerType 
                   //instead of ModelPart::NodesContainerType
                   valid_nodes[pos] = valid_nodes[valid_nodes_length - 1];
                   valid_nodes_length = valid_nodes_length - 1;
               }

               for (int i = 0; i < number_of_particles_to_insert; i++) {
                   creator.ElementCreatorWithPhysicalParameters(r_modelpart, /*r_Elem_Id*/1, inserting_nodes[i],InletModelPart.GetProperties(mesh_number));                                                                                                          
               }               
           } //if (number_of_particles_to_insert)
           mesh_iterator_number++;
        }                
    }    
    //MA
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



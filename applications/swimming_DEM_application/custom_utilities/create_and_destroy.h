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
#include "utilities/timer.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"
//SALVA_BEGINNING
// Project includes
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "../../DEM_application/custom_elements/spheric_swimming_particle.h"
#
#include "includes/define.h"
//#include "custom_elements/spheric_swimming_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../../DEM_application/DEM_application.h"

#include "../../DEM_application/custom_elements/spheric_particle.h"
#
//SALVA_ENDING
//const double prox_tol = 0.00000000001;
namespace Kratos
{

class Particle_Creator_Destructor
{
public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;


    KRATOS_CLASS_POINTER_DEFINITION(Particle_Creator_Destructor);

  
    Particle_Creator_Destructor() {};
    /// Destructor.

    virtual ~Particle_Creator_Destructor() {};

    /// Default constructor.
    
    //SALVA
    
    /*void NodeCreator(ModelPart& r_modelpart, Node < 3 > ::Pointer& pnew_node) {
              
      pnew_node = r_modelpart.CreateNewNode(10000, 1.0, 1.5, 0.0); //CASE TUBE09
      //CASE TUBE03 pnew_node = r_modelpart.CreateNewNode(10000, 0.2, 1.0, 0.0);
      pnew_node->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
      pnew_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY_X) = 0.0;
      pnew_node->FastGetSolutionStepValue(RADIUS) = 0.05;
      pnew_node->FastGetSolutionStepValue(PARTICLE_DENSITY) = 1000;
      //pnew_node->FastGetSolutionStepValue(NODAL_MASS) = 1.0; //SALVA_FALTABA?
      pnew_node->FastGetSolutionStepValue(YOUNG_MODULUS) = 1000;
      pnew_node->FastGetSolutionStepValue(POISSON_RATIO) = 0.5;
      pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = 1;
          
    }*/
    
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
      pnew_node->FastGetSolutionStepValue(YOUNG_MODULUS) = 100000;
      pnew_node->FastGetSolutionStepValue(POISSON_RATIO) = 0.5;
      pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = 1;
          
    }
    
    /*void ElementCreator(ModelPart& r_modelpart) {

      KRATOS_WATCH("Hola")

      Node < 3 > ::Pointer pnew_node;
      
      NodeCreator(r_modelpart, pnew_node);
      
      Geometry< Node < 3 > >::PointsArrayType nodelist;
      
      nodelist.push_back(pnew_node);

      Element::Pointer p_swimming_particle = Element::Pointer(new SphericSwimmingParticle(11000, nodelist)); //POOYAN
      
      r_modelpart.Elements().push_back(p_swimming_particle); //POOYAN
      
      //SphericSwimmingParticle swimming_particle(11000, nodelist); //SALVA
            
      //r_modelpart.Elements().push_back(swimming_particle); //SALVA
          
}*/

    void ElementCreator(ModelPart& r_modelpart, int r_Elem_Id, int r_Id, double r_x, double r_y, double r_z) {

      Node < 3 > ::Pointer pnew_node;
      
      /*int Id=10000;
      double x= 1.0;
      double y= 1.5;*/
      
      NodeCreator(r_modelpart, pnew_node, r_Id, r_x, r_y, r_z);
      
      Geometry< Node < 3 > >::PointsArrayType nodelist;
      
      nodelist.push_back(pnew_node);

      Element::Pointer p_swimming_particle = Element::Pointer(new SphericSwimmingParticle(r_Elem_Id, nodelist)); //POOYAN
      
      r_modelpart.Elements().push_back(p_swimming_particle); //POOYAN
           
}

    void PrintingTest() {}

    //SALVA
    
    void CalculateSurroundingBoundingBox( ModelPart& r_model_part, double scale_factor)
    {


        KRATOS_TRY

        //Type definitions
        Configure::ElementsContainerType::Pointer pElements         = r_model_part.pElements();
        Configure::ElementsContainerType Elements                   = r_model_part.Elements();
       

        double ref_radius = (*(Elements.begin().base()))->GetValue(RADIUS);
        array_1d<double, 3 > coor = (*(Elements.begin().base()))->GetGeometry()(0)->Coordinates();
        mLowPoint = coor;
        mHighPoint = coor;

        for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin();
                particle_pointer_it != Elements.end(); ++particle_pointer_it)
        {

            coor = (*(particle_pointer_it.base()))->GetGeometry()(0)->Coordinates();
            for (std::size_t i = 0; i < 3; i++)
            {
                mLowPoint[i] = (mLowPoint[i] > coor[i]) ? coor[i] : mLowPoint[i];
                mHighPoint[i] = (mHighPoint[i] < coor[i]) ? coor[i] : mHighPoint[i];          
            }
        }
        array_1d<double, 3 > midpoint = 0.5 * (mHighPoint + mLowPoint);
        mHighPoint = midpoint * (1 - scale_factor) + scale_factor * mHighPoint;
        mLowPoint = midpoint * (1 - scale_factor) + scale_factor * mLowPoint;
        for (std::size_t i = 0; i < 3; i++)
        {
            mLowPoint[i] -= 2 * ref_radius;
            mHighPoint[i] += 2 * ref_radius;
        }
        Particle_Creator_Destructor::GetHighNode() = mHighPoint;
        Particle_Creator_Destructor::GetLowNode() = mLowPoint;
        
        //KRATOS_WATCH(mHighPoint)
        //KRATOS_WATCH(mLowPoint)

        KRATOS_CATCH("")
         
    }

    void DestroyDistantParticles(ModelPart& r_model_part)
    {

        KRATOS_TRY
        //Type definitions
        Configure::ElementsContainerType::Pointer pElements      = r_model_part.pElements();
        ModelPart::NodesContainerType::Pointer pNodes = r_model_part.pNodes();

        Configure::ElementsContainerType& rElements      = r_model_part.Elements();
        ModelPart::NodesContainerType& rNodes = r_model_part.Nodes();

        Configure::ElementsContainerType temp_particles_container;
        ModelPart::NodesContainerType temp_nodes_container;

        //Copy the elements and clear the element container
// 	KRATOS_WATCH(pElements->size());
// 	KRATOS_WATCH(pNodes->size());
// 	KRATOS_WATCH(temp_particles_container);
// 	KRATOS_WATCH(rNodes.size());
// 	KRATOS_WATCH(rElements.size());
        temp_particles_container.reserve(pElements->size());
        temp_nodes_container.reserve(pNodes->size());
        //KRATOS_WATCH(rElements);
	//KRATOS_WATCH(rNodes);
	temp_nodes_container.swap(rNodes);
        temp_particles_container.swap(rElements);
	
        
        //Add the ones inside the bounding box
        //Add the ones inside the bounding box
        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = temp_particles_container.ptr_begin();
                particle_pointer_it != temp_particles_container.ptr_end(); ++particle_pointer_it)
        {	
	  array_1d<double, 3 > coor = ( *particle_pointer_it )->GetGeometry()(0)->Coordinates();

            //KRATOS_WATCH(coor)
            bool include = true;

            for (std::size_t i = 0; i < 3; i++)
            {
	      include = include && (coor[i] >= mLowPoint[i]) && (coor[i] <= mHighPoint[i]);


            }


            if (include)
            {
	      (rElements).push_back(*particle_pointer_it); //adding the elements
   

               for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++) //GENERAL FOR ELEMENTS OF MORE THAN ONE NODE
               {
                   ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
		   (rNodes).push_back( pNode );
               }

            

	    
            }

            else
            {

            //KRATOS_WATCH((*(*particle_pointer_it)).Id()) KRATOS_WATCH(coor[1])  KRATOS_WATCH(mLowPoint[1]) KRATOS_WATCH(mHighPoint[1])

            }

            //rNodes.Sort(); //this makes the calculation go so slowly
            //rNodes.Unique();

            //KRATOS_WATCH(rElements.size())
            //KRATOS_WATCH((r_model_part.Elements()).size())


        }
	//KRATOS_WATCH(r_model_part);
        KRATOS_CATCH("")
       
    }

    void DestroyDistantParticlesGivenBBox( ModelPart& r_model_part, array_1d<double, 3 > low_point,
                                          array_1d<double, 3 > high_point)
    {
    
        KRATOS_TRY

        mLowPoint = low_point;
        mHighPoint = high_point;

        //Type definitions
        Configure::ElementsContainerType::Pointer pElements      = r_model_part.pElements();
        Configure::ElementsContainerType Elements      = r_model_part.Elements();
        Configure::ElementsContainerType temp_particles_container;

        //Copy the elements and clear the element container
        temp_particles_container.reserve(pElements->size());
        temp_particles_container.swap(Elements);

        //Add the ones inside the bounding box
        for (Configure::ElementsContainerType::iterator particle_pointer_it = temp_particles_container.begin();
                particle_pointer_it != temp_particles_container.end(); ++particle_pointer_it)
        {

            array_1d<double, 3 > coor = (*(particle_pointer_it.base()))->GetGeometry()(0)->Coordinates();

            bool include = true;

            for (std::size_t i = 0; i < 3; i++)
            {
                include = include && (coor[i] >= mLowPoint[i]) && (coor[i] <= mHighPoint[i]);
            }

            if (include)
            {
               (r_model_part.Elements()).push_back(*particle_pointer_it); //adding the elements
            }
        }

        KRATOS_CATCH("")
  
    }


    ///@}
    ///@name Access
    ///@{

    array_1d<double, 3 > & GetHighNode()
    {
        return (mHighPoint);
    };

    array_1d<double, 3 > & GetLowNode()
    {
        return (mLowPoint);
    };


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
    Particle_Creator_Destructor & operator=(Particle_Creator_Destructor const& rOther);


    ///@}

}; // Class Particle_Creator_Destructor

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



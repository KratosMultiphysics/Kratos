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
// Project includes
#include "includes/define.h"
#include "custom_elements/discrete_element.h"
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

    void CalculateSurroundingBoundingBox( ModelPart& r_model_part, double scale_factor)
    {


        KRATOS_TRY

        //Type definitions
        Configure::ElementsContainerType::Pointer pElements = r_model_part.pElements();
        Configure::ElementsContainerType Elements           = r_model_part.Elements();
       
        double ref_radius         = (*(Elements.begin().base()))->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
        array_1d<double, 3 > coor = (*(Elements.begin().base()))->GetGeometry()(0)->Coordinates();
        mLowPoint                 = coor;
        mHighPoint                = coor;
        mStrictLowPoint           = coor;
        mStrictHighPoint          = coor;

        for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin();
                particle_pointer_it != Elements.end(); ++particle_pointer_it){

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
        temp_particles_container.reserve(pElements->size());
        temp_nodes_container.reserve(pNodes->size());
        
        temp_particles_container.swap(rElements);
        temp_nodes_container.swap(rNodes);

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

        KRATOS_CATCH("")
       
    }
      
    
    

    void DestroyDistantParticlesGivenBBox( ModelPart& r_model_part, array_1d<double, 3 > low_point,
                                          array_1d<double, 3 > high_point)
    {
    /*
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
  */
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

    array_1d<double, 3 > & GetStrictHighNode()
    {
        return (mHighPoint);
    };

    array_1d<double, 3 > & GetStrictLowNode()
    {
        return (mLowPoint);
    };

    double & GetDiameter()
    {
        return (mDiameter);
    };

    double & GetStrictDiameter()
    {
        return (mStrictDiameter);
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
    array_1d<double, 3 > mStrictHighPoint;
    array_1d<double, 3 > mStrictLowPoint;
    double mDiameter;
    double mStrictDiameter;

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



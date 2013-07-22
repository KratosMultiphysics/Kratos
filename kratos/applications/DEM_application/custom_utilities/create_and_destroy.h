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

    /// Destructor.
    virtual ~ParticleCreatorDestructor(){}

    void CalculateSurroundingBoundingBox(ModelPart& r_model_part, double scale_factor)
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



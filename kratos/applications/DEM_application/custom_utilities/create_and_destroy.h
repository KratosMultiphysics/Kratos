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
#include "custom_utilities/particle_configure.h"
//const double prox_tol = 0.00000000001;
namespace Kratos {
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

    /** Detail class definition.
     */

    template <std::size_t TDim,
    class TParticle,
    class TParticlePointer,
    class TParticleVector,
    class TParticleWeakVector,
    class TParticlePointerVector,
    class TParticleIterator,
    class TParticleWeakIterator,
    class TParticlePointerIterator,
    class TDistanceVector,
    class TDistanceIterator
    >

    class Particle_Creator_Destructor {
    public:

        ///@name Type Definitions
        ///@{
        typedef TParticle Particle;
        typedef TParticlePointer ParticlePointer;
        typedef TParticleVector ParticleVector;
        typedef TParticleWeakVector ParticleWeakVector;
        typedef TParticlePointerVector ParticlePointerVector;
        typedef TParticleIterator ParticleIterator;
        typedef TParticleWeakIterator ParticleWeakIterator;
        typedef TParticlePointerIterator ParticlePointerIterator;
        typedef TDistanceVector DistanceVector;
        typedef TDistanceIterator DistanceIterator;

        /// Pointer definition of Particle_Creator_Destructor
        KRATOS_CLASS_POINTER_DEFINITION(Particle_Creator_Destructor);

        ///@}
        ///@name Life Cycle
        ///@{
         Particle_Creator_Destructor() {};
        /// Destructor.

        virtual ~Particle_Creator_Destructor() {};
      
        /// Default constructor.

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

    void CalculateSurroundingBoundingBox(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part, double scale_factor) {
        KRATOS_TRY
        double ref_radius = (*(vector_of_particle_pointers.begin().base()))->GetRadius();
        array_1d<double, 3 > coor = (*(vector_of_particle_pointers.begin().base()))->GetPosition();
        mLowPoint = coor;
        mHighPoint = coor;
        for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
            particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it) {
            coor = (*(particle_pointer_it.base()))->GetPosition();
            for (std::size_t i = 0; i < 3; i++) {
                mLowPoint[i] = (mLowPoint[i] > coor[i]) ? coor[i] : mLowPoint[i];
                mHighPoint[i] = (mHighPoint[i] < coor[i]) ? coor[i] : mHighPoint[i];
            }
        }
        array_1d<double, 3 > midpoint = 0.5 * (mHighPoint + mLowPoint);
        mHighPoint = midpoint * (1 - scale_factor) + scale_factor * mHighPoint;
        mLowPoint = midpoint * (1 - scale_factor) + scale_factor * mLowPoint;
        for (std::size_t i = 0; i < 3; i++){
            mLowPoint[i] -= 2 * ref_radius;
            mHighPoint[i] += 2 * ref_radius;
        }
        Particle_Creator_Destructor::GetHighNode() = mHighPoint;
        Particle_Creator_Destructor::GetLowNode() = mLowPoint;
        KRATOS_WATCH("Bounding box for the model, bigger than strict by a factor of:");
        KRATOS_WATCH(scale_factor);
        KRATOS_WATCH(ref_radius);
        KRATOS_WATCH(mLowPoint);
        KRATOS_WATCH(mHighPoint);
        KRATOS_CATCH("")
    }

    void DestroyDistantParticles(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part) {
        KRATOS_TRY
        ModelPart::NodesContainerType temp_nodes_container;
        ParticlePointerVector temp_particles_container;
        temp_nodes_container.reserve(model_part.Nodes().size());
        temp_particles_container.reserve(vector_of_particle_pointers.size());
        temp_nodes_container.swap(model_part.Nodes());
        temp_particles_container.swap(vector_of_particle_pointers);
        for (ParticlePointerIterator particle_pointer_it = temp_particles_container.begin();
                particle_pointer_it != temp_particles_container.end(); ++particle_pointer_it) {
            array_1d<double, 3 > coor = (*(particle_pointer_it.base()))->GetPosition();
            bool include = true;
            for (std::size_t i = 0; i < 3; i++) {
                include = include && (coor[i] > mLowPoint[i]) && (coor[i] < mHighPoint[i]);
            }
            if (include) {
                vector_of_particle_pointers.push_back(*(particle_pointer_it.base()));
                model_part.Nodes().push_back((*(particle_pointer_it.base()))->GetPointerToCenterNode());
            }
        }
        KRATOS_CATCH("")
    }

    void DestroyDistantParticlesGivenBBox(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part, array_1d<double, 3 > low_point,
    array_1d<double, 3 > high_point) {
        KRATOS_TRY
        mLowPoint = low_point;
        mHighPoint = high_point;
        ModelPart::NodesContainerType temp_nodes_container;
        ParticlePointerVector temp_particles_container;
        temp_nodes_container.reserve(model_part.Nodes().size());
        temp_particles_container.reserve(vector_of_particle_pointers.size());
        temp_nodes_container.swap(model_part.Nodes());
        temp_particles_container.swap(vector_of_particle_pointers);
        for (ParticlePointerIterator particle_pointer_it = temp_particles_container.begin();
                particle_pointer_it != temp_particles_container.end(); ++particle_pointer_it) {
            array_1d<double, 3 > coor = (*(particle_pointer_it.base()))->GetPosition();
            bool include = true;
            for (std::size_t i = 0; i < 3; i++) {
                include = (coor[i] > mLowPoint[i]) && (coor[i] < mHighPoint[i]);
            }
            if (include) {
                vector_of_particle_pointers.push_back(*(particle_pointer_it.base()));
                (model_part.Nodes()).push_back(particle_pointer_it->GetPointerToCenterNode());
            }
        }
        KRATOS_CATCH("")
}


        ///@}
        ///@name Access
        ///@{

        array_1d<double, 3 > & GetHighNode() {
            return (mHighPoint);
        };

        array_1d<double, 3 > & GetLowNode() {
            return (mLowPoint);
        };


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
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

        inline void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
            unsigned int buffer_size = node_it->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++) {
                //getting the data of the solution step
                double* step_data = (node_it)->SolutionStepData().Data(step);
                //copying this data in the position of the vector we are interested in
                for (int j = 0; j < step_data_size; j++) {
                    step_data[j] = 0.0;
                }
            }
        }

        inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
            array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            noalias(Aux_var) = ZeroVector(3);
        }

        inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {
            double& Aux_var = (particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
            Aux_var = 0.0;
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



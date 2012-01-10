//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_NEIGHBOURS_CALCULATOR )
#define  KRATOS_NEIGHBOURS_CALCULATOR
#define KD_TREE 0
#define STATIC_BINS 1
#define DYNAMIC_BINS 0

// /* External includes */

// System includes

// Project includes

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/particle_configure.h"

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

    class Neighbours_Calculator {
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
        typedef std::vector<array_1d<double, 3 > > TangDisplacementsVectorType;
        typedef std::vector<array_1d<double, 3 > >::iterator TangDisplacementsIteratorType;
        //**************************************************************************************************************************************************************
        // Bucket types
        typedef Bucket < TDim, Particle, ParticlePointerVector> BucketType;
        //        typedef Bucket < TDim, Particle, ParticlePointerVector, ParticlePointer, ParticlePointerIterator, DistanceIterator > BucketType;
        //        typedef Bins < TDim, Particle, ParticleVector, ParticlePointer, ParticleIterator, DistanceIterator > StaticBins;
        //        typedef BinsDynamic < TDim, Particle, ParticleVector, ParticlePointer, ParticleIterator, DistanceIterator > DynamicBins;
        typedef ParticleConfigure < TParticle > ConfigureType;
        //**************************************************************************************************************************************************************

        //**************************************************************************************************************************************************************
        // DynamicBins;
#if KD_TREE == 1
        typedef Tree< KDTreePartition<BucketType> > tree;
#endif
        //    typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
        //    typedef Tree< StaticBins > tree;                                  //Binstree;
        //        typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
        //    typedef typename KdtreeBins::Partitions SubPartitions;
        //    typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
        //    typedef Bins< TDim, PointType, stdPointVector> stdBins;
        //        typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/
#if DYNAMIC_BINS == 1
        typedef BinsObjectDynamic <ConfigureType> bins;
#endif
#if STATIC_BINS == 1
        typedef BinsObjectStatic <ConfigureType> bins;
#endif
        /// Pointer definition of Neighbour_calculator
        KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        Neighbours_Calculator() {
        };
        /// Destructor.

        virtual ~Neighbours_Calculator() {
        };

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        static void Search_Neighbours(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part, double max_radius, double prox_tol) {
            KRATOS_TRY
            ParticlePointerVector aux_list_of_particles;
            boost::timer kdtree_construction;


            //**************************************************************************************************************************************************************

            //**************************************************************************************************************************************************************
            //create a spatial database with the list of new particles
            for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
                    particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it) {
                ParticlePointer pParticle = *(particle_pointer_it.base());
                //putting the nodes of the destination_model part in an auxiliary list
                aux_list_of_particles.push_back(pParticle);
            }
            //**************************************************************************************************************************************************************

            unsigned int MaximumNumberOfResults = 100;
//            int step_data_size = model_part.GetNodalSolutionStepDataSize();
//            double radius = 2.0 * max_radius;
            ParticlePointerVector Results(MaximumNumberOfResults);
            DistanceVector ResultsDistances(MaximumNumberOfResults);
//            unsigned int bucket_size = 20;
#if KD_TREE == 1
            tree nodes_tree(aux_list_of_particles.begin(), aux_list_of_particles.end(), bucket_size);
#endif
#if KD_TREE == 0
            bins particle_bin(aux_list_of_particles.begin(), aux_list_of_particles.end());
#endif
            boost::timer search_time;
            //**************************************************************************************************************************************************************

            //loop over all of the particles in the list to perform search
            //**************************************************************************************************************************************************************
            for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
                    particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it) {
                //find all of the new particles within the radius
                //looks which of the new particles is inside the radius around the working particle
#if KD_TREE == 1
                (*particle_pointer_it)->GetNumberOfNeighbours() = nodes_tree.SearchInRadius(**particle_pointer_it, radius, Results.begin(), ResultsDistances.begin(), MaximumNumberOfResults);
#endif
#if KD_TREE == 0
                typename ConfigureType::ResultIteratorType results_begin = Results.begin();
                (*particle_pointer_it)->GetNumberOfNeighbours() = particle_bin.SearchObjects(*(particle_pointer_it.base()), results_begin, MaximumNumberOfResults) - 1;
#endif
                int n_neighbours = (*particle_pointer_it)->GetNumberOfNeighbours();
                int neighbour_counter = -1;
                (*particle_pointer_it)->GetNeighbours().clear();
                for (ParticlePointerIterator neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it) {
                    if ((*particle_pointer_it)->GetPointerToCenterNode() != (*neighbour_it)->GetPointerToCenterNode()) {
                        (*particle_pointer_it)->GetNeighbours().push_back(*neighbour_it);
                    }
                    ++neighbour_counter;
                }
                //            ParticlePointerIterator neighbour_iter = Results.begin();
                //            neighbour_counter = 0;
                //            (*particle_pointer_it)->GetDistancesToNeighbours().clear();
                //            for (DistanceIterator distance_it = ResultsDistances.begin();
                //                distance_it != ResultsDistances.end() && neighbour_counter != int(n_neighbours); ++distance_it){
                //                ++neighbour_counter;
                //                ++neighbour_iter;
                //                if (*distance_it > prox_tol){
                //                    (*particle_pointer_it)->GetDistancesToNeighbours().push_back(*distance_it);
                //                    }
                //                }
                //            KRATOS_WATCH((*particle_pointer_it)->Id());
                //            KRATOS_WATCH((*particle_pointer_it)->GetNumberOfNeighbours());
            }
            //**************************************************************************************************************************************************************

            //            std::cout << "search_time time " << search_time.elapsed() << std::endl;
            KRATOS_CATCH("")
        }// Search_Neighbours


//         static void Search_Contacting_Neighbours(ParticlePointerVector& vector_of_particle_pointers){
//            KRATOS_TRY
//            array_1d<double, 3 > zero_vec;
//            zero_vec[0] = 0.0;
//            zero_vec[1] = 0.0;
//            zero_vec[2] = 0.0;
//            for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
//            particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it) {
//                //typename ConfigureType::ResultIteratorType results_begin = Results.begin();
//                //(*particle_pointer_it)->GetNumberOfNeighbours() = particle_bin.SearchObjects(*(particle_pointer_it.base()), results_begin, MaximumNumberOfResults);
//                int n_neighbours = (*particle_pointer_it)->GetNumberOfNeighbours();
//                ParticleWeakVector neighbours = (*particle_pointer_it)->GetNeighbours();
//                ParticleWeakVector& contact_neighbours = (*particle_pointer_it)->GetContactingNeighbours();
//                ParticleWeakVector old_contact_neighbours = (*particle_pointer_it)->GetContactingNeighbours();
//                TangDisplacementsVectorType tangential_displacements;
//                (*particle_pointer_it)->GetTangentialDisplacementOfNeighbours(tangential_displacements);
//                (*particle_pointer_it)->GetContactingNeighbours().clear();
//                tangential_displacements.clear();
//                double prox_tol = (*particle_pointer_it)->GetProxTol();
//                int neighbour_counter = 0;
//                for (ParticleWeakIterator neighbour_it = neighbours.begin();
//                neighbour_it != neighbours.end() && neighbour_counter != n_neighbours; ++neighbour_it) {
//                    ++neighbour_counter;
////                    if (*particle_pointer_it != *neighbour_it) {
////                        (*particle_pointer_it)->GetNeighbours().push_back(*neighbour_it);
////                    }
//                    array_1d<double, 3 > other_to_me_vect = (*particle_pointer_it)->GetPosition() - (*neighbour_it)->GetPosition();
//                    double radius = (*particle_pointer_it)->GetRadius();
//                    double other_radius = (*neighbour_it)->GetRadius();
//                    double distance_2 = other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2];
//                    double radius_sum = radius + other_radius;
//                    double aux_1 = distance_2 - radius_sum * radius_sum;
//                    bool found = false;
//                    TangDisplacementsIteratorType i_tang = tangential_displacements.begin();
//                    if (aux_1 < prox_tol){
//                        contact_neighbours.push_back(*(neighbour_it.base()));
//                        for (ParticleWeakIterator i_old_contact = old_contact_neighbours.begin(); i_old_contact != old_contact_neighbours.end(); i_old_contact++){
//                            if(*(neighbour_it.base()) == *(i_old_contact.base())){
//                                tangential_displacements.push_back(*(i_tang.base()));
//                                found = true;
//                                break;
//                            }
//                            ++i_tang;
//                        }
//                        if(!found){
//                            tangential_displacements.push_back(zero_vec);
//                            }
//                    }
//                }
//            }
//            KRATOS_CATCH("")
//         }// Search_Contacting_Neighbours
        ///@}
        ///@name Access
        ///@{


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
        Neighbours_Calculator & operator=(Neighbours_Calculator const& rOther);


        ///@}

    }; // Class Neighbours_calculator

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

#endif // KRATOS_NEIGHBOURS_CALCULATOR  defined 



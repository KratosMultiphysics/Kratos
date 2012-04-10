//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change
//something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_NEIGHBOURS_CALCULATOR )
#define  KRATOS_NEIGHBOURS_CALCULATOR

//M: we are using static bins for objects...

#include "includes/define.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/particle_configure.h"

namespace Kratos {
 
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
        typedef ParticleConfigure < TParticle > ConfigureType;
     
        typedef BinsObjectStatic <ConfigureType> bins;  //static Bins..?

        /// Pointer definition of Neighbour_calculator
        KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);

        Neighbours_Calculator() {
        };
        /// Destructor.

        virtual ~Neighbours_Calculator() {
        };

        static void Search_Neighbours(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part, double max_radius, double prox_tol)
        {
            KRATOS_TRY
            ParticlePointerVector aux_list_of_particles;
            boost::timer kdtree_construction;
            //**************************************************************************************************************************************************************

            //**************************************************************************************************************************************************************
            //create a spatial database with the list of new particles
            for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
                 particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it)

                {
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
            bins particle_bin(aux_list_of_particles.begin(), aux_list_of_particles.end());
            boost::timer search_time;
            //**************************************************************************************************************************************************************

            //loop over all of the particles in the list to perform search
            for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
                particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it)
            {
                //find all of the new particles within the radius
                //looks which of the new particles is inside the radius around the working particle

                typename ConfigureType::ResultIteratorType results_begin = Results.begin();
                (*particle_pointer_it)->GetNumberOfNeighbours() = particle_bin.SearchObjects(*(particle_pointer_it.base()), results_begin, MaximumNumberOfResults) - 1;

       // M: ARA VE LU DEL CFENG.         
  /////////////////////////////////////////////////////////////////////////////////////

                // SAVING THE OLD NEIGHBOURS, FORCES, FAILURE TYPES AND NUMBER OF NEIGHBOURS.

                ParticleWeakVector TempNeighbours;
                TempNeighbours.swap((*particle_pointer_it)->GetNeighbours());

                std::vector< array_1d<double,3> > TempContactForce;
                TempContactForce.swap((*particle_pointer_it)->GetContactForces());

                std::vector< int > TempContactFailureId;
                TempContactFailureId.swap((*particle_pointer_it)->GetContactFailureId());

                std::vector< double > TempInitialDelta;
                TempInitialDelta.swap((*particle_pointer_it)->GetContactInitialDelta());
                //M:in general we don't search here but we need it for the neigbours that the search calculator doesnt find.

                if(TempNeighbours.size() != TempContactForce.size())
                {
                    KRATOS_WATCH("Neighbour size and Contact force size do not match!");
                    KRATOS_WATCH(TempNeighbours.size());
                    KRATOS_WATCH(TempContactForce.size());
                }

                int n_neighbours = (*particle_pointer_it)->GetNumberOfNeighbours();

                // CLEARING AND INITIALITZING.

                int neighbour_counter = -1;
                (*particle_pointer_it)->GetNeighbours().clear();
                (*particle_pointer_it)->GetContactForces().clear();
		(*particle_pointer_it)->GetContactFailureId().clear();
                (*particle_pointer_it)->GetContactInitialDelta().clear();

                // GETTING NEW NEIGHBOURS

                for (ParticlePointerIterator neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
                {
                    if ((*particle_pointer_it)->GetPointerToCenterNode() != (*neighbour_it)->GetPointerToCenterNode())  // Excluding the particle as neighbour itself
                        {
                            (*particle_pointer_it)->GetNeighbours().push_back(*neighbour_it);

                            // LOOP TO EXTEND THE VECTORS AND SET A 0.0 VALUE EACH TIME

                            size_t Notemp = (*particle_pointer_it)->GetContactForces().size(); //Notemp its an index that will be amplified for every step of the loop as the GetContactForces().size() is getting larger

                            (*particle_pointer_it)->GetContactForces().resize(Notemp + 1);
                            (*particle_pointer_it)->GetContactForces()[Notemp][0] = 0.0;
                            (*particle_pointer_it)->GetContactForces()[Notemp][1] = 0.0;
                            (*particle_pointer_it)->GetContactForces()[Notemp][2] = 0.0;

                            (*particle_pointer_it)->GetContactFailureId().resize(Notemp + 1);
                            (*particle_pointer_it)->GetContactFailureId()[Notemp] = 1;

                            (*particle_pointer_it)->GetContactInitialDelta().resize(Notemp + 1);
                            (*particle_pointer_it)->GetContactInitialDelta()[Notemp] = 0.0;

                            // LOOP OVER THE OLD NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN EXISTING ONE AND COPYING THE OLD DATA

                            int OldNeighbourCounter = 0;
                            for (ParticleIterator old_neighbour = TempNeighbours.begin(); old_neighbour != TempNeighbours.end(); old_neighbour++)
                            {
                                if ( (old_neighbour.base())->expired() == false )
                                {
                                    if ((*neighbour_it)->GetPointerToCenterNode() == old_neighbour->GetPointerToCenterNode())  // IF IT'S AN EXISTING NEIGHBOUR
                                    {
                                        (*particle_pointer_it)->GetContactForces()[Notemp][0] = TempContactForce[OldNeighbourCounter][0];
                                        (*particle_pointer_it)->GetContactForces()[Notemp][1] = TempContactForce[OldNeighbourCounter][1];
                                        (*particle_pointer_it)->GetContactForces()[Notemp][2] = TempContactForce[OldNeighbourCounter][2];

                                        (*particle_pointer_it)->GetContactFailureId()[Notemp] = TempContactFailureId[OldNeighbourCounter];

                                        break;
                                    }
                                }

                                OldNeighbourCounter++;
                            } //loop old neighbours

                            //M: CORRECT IT!
                            //definir deltas si o no:
                            bool delta_OPTION = true;

                            if(delta_OPTION) {
                                // LOOP OVER THE INITIAL NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN INITIAL ONE AND THEN COPYING THE DELTA DATA
                                int InitialNeighboursCounter = 0;

                                for (ParticleIterator ini_neighbour = (*particle_pointer_it)->GetInitialNeighbours().begin(); ini_neighbour != (*particle_pointer_it)->GetInitialNeighbours().end(); ini_neighbour++)
                                {

                                  
                                  
                                    if ( (ini_neighbour.base())->expired() == false )
                                    {
                                        if ((*neighbour_it)->GetPointerToCenterNode() == ini_neighbour->GetPointerToCenterNode())  // IF IT'S AN INITIAL NEIGHBOUR
                                        {
                                            //KRATOS_WATCH( (*particle_pointer_it)->GetInitialDelta()[InitialNeighboursCounter] )

                                            (*particle_pointer_it)->GetContactInitialDelta()[Notemp] =  (*particle_pointer_it)->GetInitialDelta()[InitialNeighboursCounter];

                                            break;

                                        }
                                    }
                                InitialNeighboursCounter++;
                                } // for initial neighbours
                            }//deltaOPTION
                        }
                    ++neighbour_counter;
                } // for each neighbour, neighbour_it.

               
                    //ADDING NOT FOUND NEIGHBOURS (the ones with negative identation still in tensile contact are not detected, but they are on the old neighbours list).

                    int TempNeighbourCounter = 0;
                
                    for (ParticleIterator temp_neighbour = TempNeighbours.begin(); temp_neighbour != TempNeighbours.end(); temp_neighbour++)
                    {

                        if (TempContactFailureId[TempNeighbourCounter] == 0) // if they are not detached.
                        {
                            if ( (temp_neighbour.base())->expired() == false )
                            {

                                bool AlreadyAdded = false; //identifying if they are already found ot not.
                       
                                for(ParticleIterator new_neighbour  = (*particle_pointer_it)->GetNeighbours().begin();
                                         new_neighbour != (*particle_pointer_it)->GetNeighbours().end(); new_neighbour++)
                                {
                                    if (new_neighbour->GetPointerToCenterNode() == (temp_neighbour)->GetPointerToCenterNode())
                                    {
                                        AlreadyAdded = true; //for the ones already found in the new search.
                                        break;
                                    }
                                }

                                if (AlreadyAdded == false) //for the ones not included!
                                {

                                    (*particle_pointer_it)->GetNeighbours().push_back( TempNeighbours(TempNeighbourCounter) ); //adding the not found neighbours.

                                    size_t Notemp = (*particle_pointer_it)->GetContactForces().size();
                                    (*particle_pointer_it)->GetContactForces().resize(Notemp + 1);  // adding one more space for every missing neighbour.
                                    (*particle_pointer_it)->GetContactForces()[Notemp][0] = TempContactForce[TempNeighbourCounter][0]; //copying properties.
                                    (*particle_pointer_it)->GetContactForces()[Notemp][1] = TempContactForce[TempNeighbourCounter][1];
                                    (*particle_pointer_it)->GetContactForces()[Notemp][2] = TempContactForce[TempNeighbourCounter][2];

                                    (*particle_pointer_it)->GetContactFailureId().resize(Notemp + 1);
                                    (*particle_pointer_it)->GetContactFailureId()[Notemp] = TempContactFailureId[TempNeighbourCounter];

                                    //M:correct it!
                                    bool delta_OPTION=false;
                                    if(delta_OPTION) {
                                    (*particle_pointer_it)->GetContactInitialDelta().resize(Notemp + 1);
                                    (*particle_pointer_it)->GetContactInitialDelta()[Notemp] = TempInitialDelta[TempNeighbourCounter];
                                }

                                    (*particle_pointer_it)->GetNumberOfNeighbours()++;

                               }
                            }
                        }

                        TempNeighbourCounter++;
                     }
                 }//Loop for evey particle as a base.
               KRATOS_CATCH("")
        }// Search_Neighbours


        virtual std::string Info() const {
            return "neighbour_calculator";
        }

       

        virtual void PrintInfo(std::ostream& rOStream) const {
        }

      

        virtual void PrintData(std::ostream& rOStream) const {
        }



    protected:
    

    private:
    

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

     
        Neighbours_Calculator & operator=(Neighbours_Calculator const& rOther);


    }; // Class Neighbours_calculator


} // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_CALCULATOR  defined 



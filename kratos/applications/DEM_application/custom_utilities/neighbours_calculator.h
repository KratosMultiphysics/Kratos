//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
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
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/discrete_particle_configure.h"   //M: al altre hi tens lo del charlie
#include "containers/weak_pointer_vector.h"
#include "containers/pointer_vector.h"
#include "containers/pointer_vector_set.h"

<<<<<<< .mine
namespace Kratos {
    template<                    //pot no compilar en windows aquest tipus d'assignacio per template.
    std::size_t TDim,
    class TParticle
    >
=======
namespace Kratos
{
>>>>>>> .r5067

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

<<<<<<< .mine
            typedef DiscreteParticleConfigure < 3 > ConfigureType;

            
            typedef TParticle                                       Particle;                  // es el objecte
            typedef typename Particle::Pointer                      ParticlePointer;           // es punter al objecte
            typedef ModelPart::ElementsContainerType::ContainerType ParticleVector;            // es un vector d'objectes.
            typedef ParticleVector::iterator                        ParticleIterator;          // es un iterador d'objectes

            typedef ModelPart::ElementsContainerType                ParticlePointerVector;     // es un vector de iteradors
            typedef ParticlePointerVector::iterator                 ParticlePointerIterator;   // es un iterador de punters
            


            typedef  ConfigureType::PointType                       PointType;
            typedef  ConfigureType::DistanceIteratorType            DistanceIteratorType;
            typedef  ConfigureType::ContainerType                   ContainerType;
            typedef  ConfigureType::PointerType                     PointerType;
            typedef  ConfigureType::IteratorType                    IteratorType;   // iterador de punteros.
            typedef  ConfigureType::ResultContainerType             ResultContainerType;
            typedef  ConfigureType::ResultPointerType               ResultPointerType;
            typedef  ConfigureType::ResultIteratorType              ResultIteratorType;
            typedef  ConfigureType::ContactPairType                 ContactPairType;
            typedef  ConfigureType::ContainerContactType            ContainerContactType;
            typedef  ConfigureType::IteratorContactType             IteratorContactType;
            typedef  ConfigureType::PointerContactType              PointerContactType;
            typedef  ConfigureType::PointerTypeIterator             PointerTypeIterator;

            typedef WeakPointerVector<Element>                      ParticleWeakVector;
            typedef typename ParticleWeakVector::iterator           ParticleWeakIterator;

/*
            //typedef PointerVector<Particle>                    ParticleVector;
            typedef WeakPointerVector<Particle>                ParticleWeakVector;
            typedef PointerVector<ParticlePointer>             ParticlePointerVector;
            typedef typename ParticleVector::iterator          ParticleIterator;
            typedef typename ParticleWeakVector::iterator      ParticleWeakIterator;
            //typedef typename ParticlePointerVector::iterator   ParticlePointerIterator;

 * 
 */
            typedef std::vector<double>                        DistanceVector;
            typename DistanceVector::iterator                  DistanceIterator;

            

        typedef std::vector<array_1d<double, 3 > > TangDisplacementsVectorType;
        typedef TangDisplacementsVectorType::iterator TangDisplacementsIteratorType;
        //**************************************************************************************************************************************************************
        // Bucket types
        typedef Bucket < TDim, Particle, ParticlePointerVector> BucketType;
       
     
        typedef BinsObjectDynamic <ConfigureType> bins;  //static Bins..?
=======
class Neighbours_Calculator
{
public:
>>>>>>> .r5067

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

    typedef BinsObjectDynamic <ConfigureType> bins;  //static Bins..?

    /// Pointer definition of Neighbour_calculator
    KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);

<<<<<<< .mine
        static void Search_Neighbours(ModelPart& r_model_part)
        {

            KRATOS_TRY
            //typedef ModelPart::ElementsContainerType LocalParticleVector;
=======
    Neighbours_Calculator()
    {
    };
    /// Destructor.
>>>>>>> .r5067

<<<<<<< .mine
            ProcessInfo& CurrentProcessInfo     = r_model_part.GetProcessInfo();
            ContainerType& pElements            = r_model_part.ElementsArray();
            const double radius_extend          = CurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
            const bool delta_OPTION             = CurrentProcessInfo[DELTA_OPTION];
           // bool delta_OPTION = true;
=======
    virtual ~Neighbours_Calculator()
    {
    };
>>>>>>> .r5067

<<<<<<< .mine
=======
    static void Search_Neighbours(ParticlePointerVector& vector_of_particle_pointers, ModelPart& model_part, double search_radius)
    {
        KRATOS_TRY
        ParticlePointerVector aux_list_of_particles;
        boost::timer kdtree_construction;
        //**************************************************************************************************************************************************************
>>>>>>> .r5067

<<<<<<< .mine
             

            boost::timer kdtree_construction;
            
             
            unsigned int MaximumNumberOfResults = 100;
=======
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
>>>>>>> .r5067
//            int step_data_size = model_part.GetNodalSolutionStepDataSize();
//            double radius = 2.0 * max_radius;
<<<<<<< .mine
            ResultContainerType Results(MaximumNumberOfResults);
            DistanceVector ResultsDistances(MaximumNumberOfResults);
=======
        ParticlePointerVector Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);
>>>>>>> .r5067
//            unsigned int bucket_size = 20;
<<<<<<< .mine
            bins particle_bin(pElements.begin(),  pElements.end());
            boost::timer search_time;
            //**************************************************************************************************************************************************************
=======
        bins particle_bin(aux_list_of_particles.begin(), aux_list_of_particles.end());
        boost::timer search_time;
        //**************************************************************************************************************************************************************
>>>>>>> .r5067

<<<<<<< .mine
            
            //loop over all of the particles in the list to perform search
            for (IteratorType particle_pointer_it = pElements.begin();
                particle_pointer_it != pElements.end(); ++particle_pointer_it)
=======
        //loop over all of the particles in the list to perform search
        for (ParticlePointerIterator particle_pointer_it = vector_of_particle_pointers.begin();
                particle_pointer_it != vector_of_particle_pointers.end(); ++particle_pointer_it)
>>>>>>> .r5067
<<<<<<< .mine
            {

            
                Element::GeometryType& geom         =(*particle_pointer_it)->GetGeometry();
                double search_radius = (1+radius_extend) *geom(0)->GetSolutionStepValue(RADIUS);
                 
            
                //find all of the new particles within the radius
                //looks which of the new particles is inside the radius around the working particle
=======
        {
            //find all of the new particles within the radius
            //looks which of the new particles is inside the radius around the working particle
>>>>>>> .r5067

<<<<<<< .mine
                ResultIteratorType results_begin = Results.begin();
=======
            typename ConfigureType::ResultIteratorType results_begin = Results.begin();
            typename ConfigureType::DistanceIteratorType result_distances_begin = ResultsDistances.begin();
            (*particle_pointer_it)->GetNumberOfNeighbours() = particle_bin.SearchObjectsInRadius(*(particle_pointer_it.base()), search_radius, results_begin, result_distances_begin, MaximumNumberOfResults) - 1;
>>>>>>> .r5067

<<<<<<< .mine
                DistanceIteratorType result_distances_begin = ResultsDistances.begin();
                ///WARNING: particle_pointer_it  funcionava també amb .base()

                (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS) = particle_bin.SearchObjectsInRadius(*(particle_pointer_it), search_radius, results_begin, result_distances_begin, MaximumNumberOfResults) - 1;

                // SAVING THE OLD NEIGHBOURS, FORCES, FAILURE TYPES AND NUMBER OF NEIGHBOURS.
=======
            // M: ARA VE LU DEL CFENG.
            /////////////////////////////////////////////////////////////////////////////////////
>>>>>>> .r5067
            
                
<<<<<<< .mine
                ParticleWeakVector TempNeighbours;
                TempNeighbours.swap((*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS));
=======
            // SAVING THE OLD NEIGHBOURS, FORCES, FAILURE TYPES AND NUMBER OF NEIGHBOURS.
>>>>>>> .r5067

<<<<<<< .mine
                vector< array_1d<double,3> > TempContactForce;
                TempContactForce.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES));
=======
            ParticleWeakVector TempNeighbours;
            TempNeighbours.swap((*particle_pointer_it)->GetNeighbours());
>>>>>>> .r5067

<<<<<<< .mine
                vector< int > TempContactFailureId;
                TempContactFailureId.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID));
=======
            std::vector< array_1d<double,3> > TempContactForce;
            TempContactForce.swap((*particle_pointer_it)->GetContactForces());
>>>>>>> .r5067

<<<<<<< .mine
                vector< double > TempInitialDelta;
                TempInitialDelta.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA));
                //M:in general we don't search here but we need it for the neigbours that the search calculator doesnt find.
=======
            std::vector< int > TempContactFailureId;
            TempContactFailureId.swap((*particle_pointer_it)->GetContactFailureId());
>>>>>>> .r5067

            std::vector< double > TempInitialDelta;
            TempInitialDelta.swap((*particle_pointer_it)->GetContactInitialDelta());
            //M:in general we don't search here but we need it for the neigbours that the search calculator doesnt find.

<<<<<<< .mine
                int n_neighbours = (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS);   //M: number of neighbours es guarda automaticament???
                // com obting el NUM of NEIGHBOURS?? l'he de guardar o no????
                              
=======
            if(TempNeighbours.size() != TempContactForce.size())
            {
                KRATOS_WATCH("Neighbour size and Contact force size do not match!");
                KRATOS_WATCH(TempNeighbours.size());
                KRATOS_WATCH(TempContactForce.size());
            }

>>>>>>> .r5067
            int n_neighbours = (*particle_pointer_it)->GetNumberOfNeighbours();

<<<<<<< .mine
                int neighbour_counter = -1;
                (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).clear();
		(*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA).clear();
=======
            // CLEARING AND INITIALITZING.
>>>>>>> .r5067

<<<<<<< .mine
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(n_neighbours);  //M: COTELA: COM INICIALITZO AMB 3 ESPAIS VECTOR DE 3 DOUBLES.
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(n_neighbours);
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA).resize(n_neighbours);

                // GETTING NEW NEIGHBOURS
=======
            int neighbour_counter = -1;
            (*particle_pointer_it)->GetNeighbours().clear();
            (*particle_pointer_it)->GetContactForces().clear();
            (*particle_pointer_it)->GetContactFailureId().clear();
            (*particle_pointer_it)->GetContactInitialDelta().clear();
>>>>>>> .r5067
<<<<<<< .mine
                
                for (ResultIteratorType neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
=======

            // GETTING NEW NEIGHBOURS

            for (ParticlePointerIterator neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
            {
                if ((*particle_pointer_it)->GetPointerToCenterNode() != (*neighbour_it)->GetPointerToCenterNode())  // Excluding the particle as neighbour itself
>>>>>>> .r5067
                {
<<<<<<< .mine
                   // if ((*particle_pointer_it)->Id() != (*neighbour_it)->Id() )  //MIQUEL FES LA COMPROBACIÓ A POSTERIORI DE SI ES TROBEN O NO.
=======
                    (*particle_pointer_it)->GetNeighbours().push_back(*neighbour_it);
>>>>>>> .r5067

<<<<<<< .mine
                    {
                           (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);

                            // LOOP TO EXTEND THE VECTORS AND SET A 0.0 VALUE EACH TIME
=======
                    // LOOP TO EXTEND THE VECTORS AND SET A 0.0 VALUE EACH TIME
>>>>>>> .r5067

<<<<<<< .mine
                            size_t Notemp = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).size(); //Notemp its an index that will be amplified for every step of the loop as the NEIGHBOUR_ELEMENTS.size() is getting larger
                            
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp] = ZeroVector(3);
                                                     
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp] = 1;
                           
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA)[Notemp] = 0.0;
=======
                    size_t Notemp = (*particle_pointer_it)->GetContactForces().size(); //Notemp its an index that will be amplified for every step of the loop as the GetContactForces().size() is getting larger
>>>>>>> .r5067

<<<<<<< .mine
=======
                    (*particle_pointer_it)->GetContactForces().resize(Notemp + 1);
                    (*particle_pointer_it)->GetContactForces()[Notemp][0] = 0.0;
                    (*particle_pointer_it)->GetContactForces()[Notemp][1] = 0.0;
                    (*particle_pointer_it)->GetContactForces()[Notemp][2] = 0.0;

                    (*particle_pointer_it)->GetContactFailureId().resize(Notemp + 1);
                    (*particle_pointer_it)->GetContactFailureId()[Notemp] = 1;

                    (*particle_pointer_it)->GetContactInitialDelta().resize(Notemp + 1);
                    (*particle_pointer_it)->GetContactInitialDelta()[Notemp] = 0.0;

>>>>>>> .r5067
                    // LOOP OVER THE OLD NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN EXISTING ONE AND COPYING THE OLD DATA
                    
<<<<<<< .mine
                            int OldNeighbourCounter = 0;
                            for (ParticleWeakIterator old_neighbour = TempNeighbours.begin(); old_neighbour != TempNeighbours.end(); old_neighbour++)
=======
                    int OldNeighbourCounter = 0;
                    for (ParticleIterator old_neighbour = TempNeighbours.begin(); old_neighbour != TempNeighbours.end(); old_neighbour++)
                    {
                        if ( (old_neighbour.base())->expired() == false )
                        {
                            if ((*neighbour_it)->GetPointerToCenterNode() == old_neighbour->GetPointerToCenterNode())  // IF IT'S AN EXISTING NEIGHBOUR
>>>>>>> .r5067
                            {
<<<<<<< .mine
                                if ( (old_neighbour.base())->expired() == false )
                                {
                                    //if ((*neighbour_it)->Id() == old_neighbour->Id() )  // MIQUEL COMPROBA SI EL TROBES O NO
                                    {
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][0] = TempContactForce[OldNeighbourCounter][0];
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][1] = TempContactForce[OldNeighbourCounter][1];
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][2] = TempContactForce[OldNeighbourCounter][2];
=======
                                (*particle_pointer_it)->GetContactForces()[Notemp][0] = TempContactForce[OldNeighbourCounter][0];
                                (*particle_pointer_it)->GetContactForces()[Notemp][1] = TempContactForce[OldNeighbourCounter][1];
                                (*particle_pointer_it)->GetContactForces()[Notemp][2] = TempContactForce[OldNeighbourCounter][2];
>>>>>>> .r5067

<<<<<<< .mine
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp] = TempContactFailureId[OldNeighbourCounter];
=======
                                (*particle_pointer_it)->GetContactFailureId()[Notemp] = TempContactFailureId[OldNeighbourCounter];
>>>>>>> .r5067

                                break;
                            }
                        }

                        OldNeighbourCounter++;
                    } //loop old neighbours

<<<<<<< .mine
                          
=======
                    //M: CORRECT IT!
                    //definir deltas si o no:
                    bool delta_OPTION = true;

>>>>>>> .r5067
                    if(delta_OPTION)
                    {
                        // LOOP OVER THE INITIAL NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN INITIAL ONE AND THEN COPYING THE DELTA DATA
                        int InitialNeighboursCounter = 0;
<<<<<<< .mine
                            
                                for (ParticleWeakIterator ini_neighbour = ((*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS)).begin(); ini_neighbour != ((*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS)).end(); ini_neighbour++)
=======

                        for (ParticleIterator ini_neighbour = (*particle_pointer_it)->GetInitialNeighbours().begin(); ini_neighbour != (*particle_pointer_it)->GetInitialNeighbours().end(); ini_neighbour++)
                        {



                            if ( (ini_neighbour.base())->expired() == false )
                            {
                                if ((*neighbour_it)->GetPointerToCenterNode() == ini_neighbour->GetPointerToCenterNode())  // IF IT'S AN INITIAL NEIGHBOUR
>>>>>>> .r5067
                                {
<<<<<<< .mine
                                    if ( (ini_neighbour.base())->expired() == false )
                                    {
                                        if ((*neighbour_it)->Id() == ini_neighbour->Id())  // IF IT'S AN INITIAL NEIGHBOUR //POOOOOYAN, ESTARA BÉ AIXO???? TINDRAN DIFERENTS IDS??? INI NEIGH TINDRA ID??
                                        {
                                            //KRATOS_WATCH( (*particle_pointer_it)->GetInitialDelta()[InitialNeighboursCounter] )
=======
                                    //KRATOS_WATCH( (*particle_pointer_it)->GetInitialDelta()[InitialNeighboursCounter] )

                                    (*particle_pointer_it)->GetContactInitialDelta()[Notemp] =  (*particle_pointer_it)->GetInitialDelta()[InitialNeighboursCounter];
>>>>>>> .r5067

<<<<<<< .mine
                                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA)[Notemp] =  (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA)[InitialNeighboursCounter];
=======
                                    break;
>>>>>>> .r5067

                                }
                            }
                            InitialNeighboursCounter++;
                        } // for initial neighbours
                    }//deltaOPTION
                }
                ++neighbour_counter;
            } // for each neighbour, neighbour_it.

<<<<<<< .mine
                                        }
                                    }
                                    
                                InitialNeighboursCounter++;
                                } // for initial neighbours
                            }//deltaOPTION
                        }
                    ++neighbour_counter;
                } // for each neighbour, neighbour_it.
=======
>>>>>>> .r5067
            
            //ADDING NOT FOUND NEIGHBOURS (the ones with negative identation still in tensile contact are not detected, but they are on the old neighbours list).

<<<<<<< .mine
                    int TempNeighbourCounter = 0;
                
                    for (ParticleWeakIterator temp_neighbour = TempNeighbours.begin(); temp_neighbour != TempNeighbours.end(); temp_neighbour++)
=======
            int TempNeighbourCounter = 0;

            for (ParticleIterator temp_neighbour = TempNeighbours.begin(); temp_neighbour != TempNeighbours.end(); temp_neighbour++)
            {

                if (TempContactFailureId[TempNeighbourCounter] == 0) // if they are not detached.
                {
                    if ( (temp_neighbour.base())->expired() == false )
>>>>>>> .r5067
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

<<<<<<< .mine
                                bool AlreadyAdded = false; //identifying if they are already found ot not.
                       
                                for(ParticleWeakIterator new_neighbour  = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).begin();
                                         new_neighbour != (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).end(); new_neighbour++)
                                {
                                    if (new_neighbour->Id() == (temp_neighbour)->Id())
                                    {
                                        AlreadyAdded = true; //for the ones already found in the new search.
                                        break;
                                    }
                                }
=======
                        if (AlreadyAdded == false) //for the ones not included!
                        {
>>>>>>> .r5067

                            (*particle_pointer_it)->GetNeighbours().push_back( TempNeighbours(TempNeighbourCounter) ); //adding the not found neighbours.

<<<<<<< .mine
                                    (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back( TempNeighbours(TempNeighbourCounter) ); //adding the not found neighbours.
=======
                            size_t Notemp = (*particle_pointer_it)->GetContactForces().size();
                            (*particle_pointer_it)->GetContactForces().resize(Notemp + 1);  // adding one more space for every missing neighbour.
                            (*particle_pointer_it)->GetContactForces()[Notemp][0] = TempContactForce[TempNeighbourCounter][0]; //copying properties.
                            (*particle_pointer_it)->GetContactForces()[Notemp][1] = TempContactForce[TempNeighbourCounter][1];
                            (*particle_pointer_it)->GetContactForces()[Notemp][2] = TempContactForce[TempNeighbourCounter][2];
>>>>>>> .r5067

<<<<<<< .mine
                                    size_t Notemp = (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).size();
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(Notemp + 1);  // adding one more space for every missing neighbour.
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][0] = TempContactForce[TempNeighbourCounter][0]; //copying properties.
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][1] = TempContactForce[TempNeighbourCounter][1];
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][2] = TempContactForce[TempNeighbourCounter][2];
=======
                            (*particle_pointer_it)->GetContactFailureId().resize(Notemp + 1);
                            (*particle_pointer_it)->GetContactFailureId()[Notemp] = TempContactFailureId[TempNeighbourCounter];
>>>>>>> .r5067

<<<<<<< .mine
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(Notemp + 1);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp] = TempContactFailureId[TempNeighbourCounter];
=======
                            //M:correct it!
                            bool delta_OPTION=false;
                            if(delta_OPTION)
                            {
                                (*particle_pointer_it)->GetContactInitialDelta().resize(Notemp + 1);
                                (*particle_pointer_it)->GetContactInitialDelta()[Notemp] = TempInitialDelta[TempNeighbourCounter];
                            }
>>>>>>> .r5067

<<<<<<< .mine
                                    if(delta_OPTION) {
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA).resize(Notemp + 1);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_INITIAL_DELTA)[Notemp] = TempInitialDelta[TempNeighbourCounter];
                                    }
=======
                            (*particle_pointer_it)->GetNumberOfNeighbours()++;
>>>>>>> .r5067

<<<<<<< .mine
                                    (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS)++;
                               }
                            }
=======
>>>>>>> .r5067
                        }
                    }
                }

<<<<<<< .mine
                        TempNeighbourCounter++;
                    }
                 }//Loop for evey particle as a base.
              

          KRATOS_CATCH("")
        }// Search_Neighbours
=======
                TempNeighbourCounter++;
            }
        }//Loop for evey particle as a base.
        KRATOS_CATCH("")
    }// Search_Neighbours
>>>>>>> .r5067


    virtual std::string Info() const
    {
        return "neighbour_calculator";
    }



    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }



    virtual void PrintData(std::ostream& rOStream) const
    {
    }



<<<<<<< .mine
        inline void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
            /*
            unsigned int buffer_size = node_it->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++) {
                //getting the data of the solution step
                double* step_data = (node_it)->SolutionStepData().Data(step);
                //copying this data in the position of the vector we are interested in
                for (int j = 0; j < step_data_size; j++) {
                    step_data[j] = 0.0;
                }
=======
protected:


private:


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
>>>>>>> .r5067
            }
             */
        }
    }

<<<<<<< .mine
        inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
            /*
            array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            noalias(Aux_var) = ZeroVector(3);
             * */
        }
=======
    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable)
    {
        array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
        noalias(Aux_var) = ZeroVector(3);
    }
>>>>>>> .r5067

<<<<<<< .mine
        inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {
            /*
            double& Aux_var = (particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
            Aux_var = 0.0;
             * */
        }
=======
    inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable)
    {
        double& Aux_var = (particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
        Aux_var = 0.0;
    }
>>>>>>> .r5067


    Neighbours_Calculator & operator=(Neighbours_Calculator const& rOther);


}; // Class Neighbours_calculator


} // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_CALCULATOR  defined 



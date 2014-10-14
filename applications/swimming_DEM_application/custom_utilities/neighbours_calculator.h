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

#define _OPENMPI 1

//M: we are using static bins for objects...

#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/weak_pointer_vector.h"
#include "containers/pointer_vector.h"
#include "containers/pointer_vector_set.h"

#include "custom_utilities/discrete_particle_configure.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "utilities/openmp_utils.h"

namespace Kratos {

    template< //pot no compilar en windows aquest tipus d'assignacio per template.
//     std::size_t TDim, FORA!
    class TParticle
    >

    class Neighbours_Calculator {
    public:
        typedef DiscreteParticleConfigure < 3 > ConfigureType;

        typedef TParticle Particle; // es el objecte
        typedef typename Particle::Pointer ParticlePointer; // es punter al objecte
        typedef ModelPart::ElementsContainerType::ContainerType ParticleVector; // es un vector d'objectes.
        typedef ParticleVector::iterator ParticleIterator; // es un iterador d'objectes

        typedef ModelPart::ElementsContainerType ParticlePointerVector; // es un vector de iteradors
        typedef ParticlePointerVector::iterator ParticlePointerIterator; // es un iterador de punters

        typedef ConfigureType::PointType PointType;
        typedef ConfigureType::DistanceIteratorType DistanceIteratorType;
        typedef ConfigureType::ContainerType ContainerType;
        typedef ConfigureType::PointerType PointerType;
        typedef ConfigureType::IteratorType IteratorType; // iterador de punteros.
        typedef ConfigureType::ResultContainerType ResultContainerType;
        typedef ConfigureType::ResultPointerType ResultPointerType;
        typedef ConfigureType::ResultIteratorType ResultIteratorType;
        typedef ConfigureType::ContactPairType ContactPairType;
        typedef ConfigureType::ContainerContactType ContainerContactType;
        typedef ConfigureType::IteratorContactType IteratorContactType;
        typedef ConfigureType::PointerContactType PointerContactType;
        typedef ConfigureType::PointerTypeIterator PointerTypeIterator;

        typedef WeakPointerVector<Element> ParticleWeakVector;
        typedef typename ParticleWeakVector::iterator ParticleWeakIterator;
        typedef ParticleWeakVector::ptr_iterator ParticleWeakIteratorType_ptr;

        typedef std::vector<double> DistanceVector;
        typename DistanceVector::iterator DistanceIterator;

        typedef std::vector<array_1d<double, 3 > > TangDisplacementsVectorType;
        typedef TangDisplacementsVectorType::iterator TangDisplacementsIteratorType;

        // Bucket types
//         typedef Bucket < TDim, Particle, ParticlePointerVector> BucketType; No ho borro pero no veig per que hem de tindre aixo aqui si no fa res mes que molestar
        typedef BinsObjectDynamic <ConfigureType> Bins;


        /// Pointer definition of Neighbour_calculator
        //KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);  R: necesitu?

        virtual ~Neighbours_Calculator() {
        };

        //Aquesta da igual si es estatica
        static void Parallel_partitioning(ModelPart& r_model_part, bool extension_option)
        {
            /* Redefine this if you are using parallelism */
            /* Perform the repartition of the model */
        }

        //Aquestas va molt malament que sigin estaticas
        virtual void Add_To_Modelpart(ModelPart& r_model_part, ResultIteratorType neighbour_it)
        {
            /* Must be redefined */
        }

        virtual void Clean_Modelpart(ModelPart& r_model_part)
        {
            /* Must be redefined */
        }

        virtual void Sort_Modelpart(ModelPart& r_model_part)
        {
            /* Must be redefined */
        }

        virtual ContainerType& Get_Elements(ModelPart& r_model_part)
        {
            /* Must be redefined */
            return r_model_part.ElementsArray();
        }

        virtual void SearchNeighbours(ModelPart& r_model_part,
                                      ContainerType& pIteratorElements,
                                      int NumberOfElements,
                                      int MaximumNumberOfResults,
                                      std::vector<std::size_t> &NumberOfResults,
                                      std::vector<std::vector<PointerType> > &Results,
                                      std::vector<std::vector<double> > &ResultsDistances,
                                      std::vector<double> &Radius
        )
        {
            Bins particle_bin(pIteratorElements.begin(), pIteratorElements.end());

            particle_bin.SearchObjectsInRadiusInner(pIteratorElements.begin(),NumberOfElements,Radius,Results,ResultsDistances,NumberOfResults,MaximumNumberOfResults);
        }

        void Search_Ini_Neighbours(ModelPart& r_model_part, bool extension_option)
        {
            KRATOS_TRY

            ContainerType& pIteratorElements = Get_Elements(r_model_part);
            ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

            double radius_extend = 0.0;
            if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];

            boost::timer kdtree_construction;

            unsigned int MaximumNumberOfResults = 1000;
            unsigned int ResultIterator = 0;

            boost::timer search_time;

            //**********************************************************************************************************************************************//

            int NumberOfElements = pIteratorElements.end() - pIteratorElements.begin();

            Clean_Modelpart(r_model_part);

            std::vector<std::size_t>               NumberOfResults(NumberOfElements);
            std::vector<std::vector<PointerType> > Results(NumberOfElements, std::vector<PointerType>(MaximumNumberOfResults));
            std::vector<std::vector<double> >      ResultsDistances(NumberOfElements, std::vector<double>(MaximumNumberOfResults));
            std::vector<double>                    Radius(NumberOfElements);

            //Radius vector fill
            for (IteratorType particle_pointer_it = pIteratorElements.begin(); particle_pointer_it != pIteratorElements.end(); ++particle_pointer_it)
            {
                Radius[particle_pointer_it - pIteratorElements.begin()] = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS); //if this is changed, then compobation before adding neighbours must change also.
            }

            SearchNeighbours(r_model_part,pIteratorElements,NumberOfElements,MaximumNumberOfResults,NumberOfResults,Results,ResultsDistances,Radius);

            //Aqui ja tenim tots els resultats de tots el elements i fem el cambi de buffers de mpi a els iteradors normals de kratos
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif

            vector<unsigned int> element_partition;
            OpenMPUtils::CreatePartition(number_of_threads, pIteratorElements.size(), element_partition);

            #pragma omp parallel for private(ResultIterator)
            for(int k=0; k<number_of_threads; k++)
            {
                IteratorType it_begin=pIteratorElements.begin()+element_partition[k];
                IteratorType it_end=pIteratorElements.begin()+element_partition[k+1];

                ResultIterator = element_partition[k];

                for (IteratorType particle_pointer_it = it_begin;
                        particle_pointer_it != it_end; ++particle_pointer_it, ++ResultIterator)
                {
                    (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = NumberOfResults[ResultIterator];

                    // CLEARING AND INITIALITZING.

                    vector< int > TempIds;
                    TempIds.swap((*particle_pointer_it)->GetValue(NEIGHBOURS_IDS));

                    (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).clear();
                    (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).clear();

                    // STORING THE NEIGHBOURS AND PROPERTIES

                    unsigned int neighbour_counter = 0;

                    for (ResultIteratorType neighbour_it = Results[ResultIterator].begin(); neighbour_counter < NumberOfResults[ResultIterator]; ++neighbour_it)
                    {
                        Add_To_Modelpart(r_model_part,neighbour_it);

                        double particle_radius  = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                        double neigh_radius     = (*neighbour_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);

                        array_1d<double,3> other_to_me_vect      = (*particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();
                        double distance                          = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                                       other_to_me_vect[1] * other_to_me_vect[1] +
                                                                       other_to_me_vect[2] * other_to_me_vect[2]);

                        double neighbour_search_radius           = (1.0 + radius_extend) * neigh_radius;

                        if( (distance - particle_radius)  <= neighbour_search_radius )
                        {

                             (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);
                             size_t size = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).size();

                             (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).resize(size);
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(size);
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(size);
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(size);
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(size);
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(size);


                             (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS)[size-1] = (*neighbour_it)->Id();

                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][0] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][1] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][2] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][0] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][1] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][2] = 0.0;
                             (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[size-1] = 1;
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[size-1] = 1;
                             (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[size-1] = 0.0;

                        }
                        ++neighbour_counter;

                    }// For each neighbour, neighbour_it.
                }// Loop for evey particle as a base.
            }// OpenMP loop

            Sort_Modelpart(r_model_part);

            KRATOS_CATCH("")
        }// Search_Ini_Neighbours

        void Search_Neighbours(ModelPart& r_model_part, bool extension_option)
        {
            KRATOS_TRY

            //KRATOS_WATCH("SEARCH_NEIGHBOURS")

            ContainerType& pIteratorElements = Get_Elements(r_model_part);
            ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

            double radius_extend = 0.0;
            if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];

            boost::timer kdtree_construction;

            unsigned int MaximumNumberOfResults = 1000;
            unsigned int ResultIterator = 0;

            boost::timer search_time;
            //**************************************************************************************************************************************************************

            int NumberOfElements = pIteratorElements.end() - pIteratorElements.begin();

            Clean_Modelpart(r_model_part);

            std::vector<std::size_t>               NumberOfResults(NumberOfElements);
            std::vector<std::vector<PointerType> > Results(NumberOfElements, std::vector<PointerType>(MaximumNumberOfResults));
            std::vector<std::vector<double> >      ResultsDistances(NumberOfElements, std::vector<double>(MaximumNumberOfResults));
            std::vector<double>                    Radius(NumberOfElements);

            ///Radius vector fill
            for (IteratorType particle_pointer_it = pIteratorElements.begin(); particle_pointer_it != pIteratorElements.end(); ++particle_pointer_it)
            {
                //####################################################################################################################
                double new_extension = 2.0; ///WARNING: PROVISIONALLY SET AS 2.0. SHOULD BE CALCULATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //####################################################################################################################

                Radius[particle_pointer_it - pIteratorElements.begin()] = new_extension*((1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS));
            }

            ///Aqui es fa la cerca
            Timer::Start("SEARCH");
            SearchNeighbours(r_model_part,pIteratorElements,NumberOfElements,MaximumNumberOfResults,NumberOfResults,Results,ResultsDistances,Radius);
            Timer::Stop("SEARCH");

            Timer::Start("PROCESS");
            ///Aqui ja tenim tots els resultats de tots el elements i fem el cambi de buffers de mpi a els iteradors normals de kratos
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif

            vector<unsigned int> element_partition;
            OpenMPUtils::CreatePartition(number_of_threads, pIteratorElements.size(), element_partition);

            #pragma omp parallel for private(ResultIterator)
            for(int k=0; k<number_of_threads; k++)
            {

                IteratorType it_begin=pIteratorElements.begin()+element_partition[k];
                IteratorType it_end=pIteratorElements.begin()+element_partition[k+1];

                ResultIterator = element_partition[k];

                for (IteratorType particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it, ++ResultIterator)
                {
                    (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = NumberOfResults[ResultIterator];

                    ParticleWeakVector TempNeighbours;

                    TempNeighbours.swap((*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS));

                    vector< int > TempIds;
                    TempIds.swap((*particle_pointer_it)->GetValue(NEIGHBOURS_IDS));

                    vector< array_1d<double, 3 > > TempContactForce;
                    TempContactForce.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES));

                    vector< int > TempContactFailureId;
                    TempContactFailureId.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID));

                    vector< double > TempContactDelta;
                    TempContactDelta.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA));

                    vector< double > TempRotateSpringFailType;
                    TempRotateSpringFailType.swap((*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE));

                    vector< array_1d<double, 3 > > TempRotateSpringMoment;
                    TempRotateSpringMoment.swap((*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT));

                    //CLEARING AND INITIALITZING

                    (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).clear();
                    (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).clear();

                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).clear();

                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).clear();
                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).clear();

                    /*
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(n_neighbours);
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(n_neighbours);
                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(n_neighbours);

                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(n_neighbours);
                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(n_neighbours);
                    */

                    // GETTING NEW NEIGHBOURS

                    double radius = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                    int InitialNeighbourSize = (*particle_pointer_it)->GetValue(INI_NEIGHBOURS_IDS).size();
                    int OldNeighbourSize = TempIds.size();
                    unsigned int neighbour_counter = 0;

                    for (ResultIteratorType neighbour_it = Results[ResultIterator].begin(); neighbour_counter < NumberOfResults[ResultIterator]; ++neighbour_it)
                    {
                        double initial_delta = 0.0;
                        int failure_id       = 1;
                        //int ini_failure      = 0;
                        double indentation   = 0.0;

                        bool already_added = false;

                        Timer::Start("ADD");
                        //Add Elements to ghost mesh //TODO: Add elements to local mesh (NYI)
                        Add_To_Modelpart(r_model_part,neighbour_it);
                        Timer::Stop("ADD");

                        for (int IniNeighbourCounter = 0; IniNeighbourCounter != InitialNeighbourSize; IniNeighbourCounter++)
                        {
                            if ( static_cast<int>((*neighbour_it)->Id()) == (*particle_pointer_it)->GetValue(INI_NEIGHBOURS_IDS)[IniNeighbourCounter]) // is this and initial neighbour?
                            {
                                initial_delta = (*particle_pointer_it)->GetValue(PARTICLE_INITIAL_DELTA)[IniNeighbourCounter];
                                //ini_failure = (*particle_pointer_it)->GetValue(PARTICLE_INITIAL_FAILURE_ID)[IniNeighbourCounter];
                            }
                        } //getting initial deltas and initial failure values.

                        for (int OldNeighbourCounter = 0; OldNeighbourCounter != OldNeighbourSize; OldNeighbourCounter++)
                        {
                            if (static_cast<int>((*neighbour_it)->Id()) == TempIds[OldNeighbourCounter]) // is this and old neighbour?
                            {
                                failure_id                          = TempContactFailureId[OldNeighbourCounter];
                                double other_radius                 = (*neighbour_it)->GetGeometry()[0].GetSolutionStepValue(RADIUS);
                                array_1d<double,3> other_to_me_vect = (*particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();
                                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                                      other_to_me_vect[1] * other_to_me_vect[1] +
                                                                      other_to_me_vect[2] * other_to_me_vect[2]);
                                double radius_sum                   = radius + other_radius;
                                double indentation                  = radius_sum - distance - initial_delta;

                                if ( indentation >= -1.0e-6 || (indentation < 1.0e-6 && failure_id == 0 ) )  //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.
                                {
                                    already_added = true;

                                    (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);

                                    size_t size = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).size();

                                    (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).resize(size);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(size);
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(size);
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(size);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(size);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(size);
                                    (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS)[size-1] = (*neighbour_it)->Id();
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][0] = TempContactForce[OldNeighbourCounter][0];
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][1] = TempContactForce[OldNeighbourCounter][1];
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][2] = TempContactForce[OldNeighbourCounter][2];
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][0] = TempRotateSpringMoment[OldNeighbourCounter][0];
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][1] = TempRotateSpringMoment[OldNeighbourCounter][1];
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][2] = TempRotateSpringMoment[OldNeighbourCounter][2];
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[size-1] = TempRotateSpringFailType[OldNeighbourCounter];
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[size-1] = failure_id;
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[size-1] = initial_delta;
                                }
                            }
                        }

                        //Not old ones

                        double other_radius = (*neighbour_it)->GetGeometry()[0].GetSolutionStepValue(RADIUS);

                        array_1d<double,3> other_to_me_vect = (*particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();

                        double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                              other_to_me_vect[1] * other_to_me_vect[1] +
                                              other_to_me_vect[2] * other_to_me_vect[2]);

                        double radius_sum = radius + other_radius;
                        indentation = radius_sum - distance - initial_delta;

                        if ( (already_added == false) && (indentation > 0.0))
                        {
                            (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);

                            size_t size = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).size();

                            (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).resize(size);
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(size);
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(size);
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(size);
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(size);
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(size);
                            (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS)[size-1] = (*neighbour_it)->Id();
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][0] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][1] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[size-1][2] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][0] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][1] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[size-1][2] = 0.0;
                            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[size-1] = 1;
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[size-1] = 1;
                            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[size-1] = initial_delta;
                        }

                        ++neighbour_counter;

                    } // for each neighbour, neighbour_it.

                }//Loop for evey particle as a base.
            }
            Timer::Start("SORT");
            Sort_Modelpart(r_model_part);
            Timer::Stop("SORT");

            Timer::Stop("PROCESS");

            Timer::PrintTimingInformation();

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
            /*
            array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            noalias(Aux_var) = ZeroVector(3);
             * */
        }

        inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {
            /*
            double& Aux_var = (particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
            Aux_var = 0.0;
             * */
        }

        Neighbours_Calculator & operator=(Neighbours_Calculator const& rOther);

    }; // Class Neighbours_calculator


} // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_CALCULATOR  defined 



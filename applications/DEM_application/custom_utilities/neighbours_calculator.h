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

// Project Includes
#include "spatial_containers/spatial_search.h"

namespace Kratos {

    template< class TParticle>
    class Neighbours_Calculator {
    public:
        
        /// SpatialSearch
        typedef TParticle Particle; // es el objecte
        typedef typename Particle::Pointer ParticlePointer; // es punter al objecte
        typedef ModelPart::ElementsContainerType::ContainerType ParticleVector; // es un vector d'objectes.
        typedef ParticleVector::iterator ParticleIterator; // es un iterador d'objectes
        
        typedef SpatialSearch::Pointer                          SearchPointer;
        
        typedef SpatialSearch::ElementsContainerType::ContainerType ElementsContainerType;
        
        typedef ElementsContainerType::value_type               PointerType;
        typedef ElementsContainerType::iterator                 IteratorType;

        typedef SpatialSearch::ResultElementsContainerType      ResultElementsContainerType;
        typedef SpatialSearch::VectorResultElementsContainerType    VectorResultElementsContainerType;
        
        typedef SpatialSearch::RadiusArrayType                  RadiusArrayType;
        typedef SpatialSearch::DistanceType                     DistanceType;
        typedef SpatialSearch::VectorDistanceType               VectorDistanceType;

        typedef WeakPointerVector<Element>                      ParticleWeakVector;

        typedef std::vector<array_1d<double, 3 > >              TangDisplacementsVectorType;
        typedef TangDisplacementsVectorType::iterator           TangDisplacementsIteratorType;

        Neighbours_Calculator(ModelPart& rModelPart) 
        : mrModelPart(rModelPart),
          mpIteratorElements(mrModelPart.GetCommunicator().LocalMesh().ElementsArray()),
          mrCurrentProcessInfo(mrModelPart.GetProcessInfo())
        {
            
        }
        
        ~Neighbours_Calculator() 
        {
        }
        
        void Initialize()
        {
            int NumberOfElements = mpIteratorElements.end() - mpIteratorElements.begin();
          
            mResults.resize(NumberOfElements);
            mResultsDistances.resize(NumberOfElements);
            mRadius.resize(NumberOfElements);
        }
        
        void Search_Ini_Neighbours(ModelPart& mrModelPartNO, bool extension_option, SearchPointer searchScheme) 
        {   
            KRATOS_TRY

            double radius_extend = 0.0;
            if (extension_option) radius_extend = mrCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
            
            //**********************************************************************************************************************************************//
            
            int NumberOfElements = mpIteratorElements.end() - mpIteratorElements.begin();

            //Radius vector fill
	    
            for (IteratorType particle_pointer_it = mpIteratorElements.begin(); particle_pointer_it != mpIteratorElements.end(); ++particle_pointer_it)
            {   
                mRadius[particle_pointer_it - mpIteratorElements.begin()] = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS); //if this is changed, then compobation before adding neighbours must change also.
            }

            searchScheme->SearchElementsInRadiusExclusive(mrModelPart,mpIteratorElements,mRadius,mResults,mResultsDistances);

            int number_of_threads = OpenMPUtils::GetNumThreads();

            vector<unsigned int> element_partition;
            OpenMPUtils::CreatePartition(number_of_threads, mpIteratorElements.size(), element_partition);
            
            std::cout << mrModelPart.GetCommunicator().GhostMesh().ElementsArray().size() << std::endl;
            
            int initial = 0;

            #pragma omp parallel for
            for(int k=0; k<number_of_threads; k++)
            {
                IteratorType it_begin=mpIteratorElements.begin()+element_partition[k];
                IteratorType it_end=mpIteratorElements.begin()+element_partition[k+1];
//                 
                int ResultCounter = element_partition[k];
                
                for (IteratorType particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it, ++ResultCounter)
                {     
                    (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = mResults[ResultCounter].size();

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
    
                    for (IteratorType neighbour_it = mResults[ResultCounter].begin(); neighbour_counter < mResults[ResultCounter].size(); ++neighbour_it)
                    {                  
                        double particle_radius  = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                        double neigh_radius     = (*neighbour_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                        
                        array_1d<double,3> other_to_me_vect      = (*particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();
                        double distance                          = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                                       other_to_me_vect[1] * other_to_me_vect[1] +
                                                                       other_to_me_vect[2] * other_to_me_vect[2]);

                        double neighbour_search_radius           = (1.0 + radius_extend) * neigh_radius;
                                    initial++;                                    
                        if( true || (distance - particle_radius)  <= neighbour_search_radius ) 
                        {
                             (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);
//                              initial++;
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
            
            std::cout << "Initial neighbours found: " << initial << std::endl;
//             abort();

            KRATOS_CATCH("")
        }// Search_Ini_Neighbours

        void Search_Neighbours(ModelPart& mrModelPartNO, bool extension_option, SearchPointer searchScheme) 
        {   
            KRATOS_TRY

            double radius_extend = 0.0;
            if (extension_option) radius_extend = mrCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];

            unsigned int ResultCounter = 0;

            int NumberOfElements = mpIteratorElements.end() - mpIteratorElements.begin();

            double new_extension = mrCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION];

            ///Radius vector fill
            for (IteratorType particle_pointer_it = mpIteratorElements.begin(); particle_pointer_it != mpIteratorElements.end(); ++particle_pointer_it)
            {
                //####################################################################################################################
                double new_extension = mrCurrentProcessInfo[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION]; ///WARNING: PROVISIONALLY SET MANUALLY SHOULD BE CALCULATED!!!!!!!!! ITS OK FOR CONCRETE
                //####################################################################################################################
              
                mRadius[particle_pointer_it - mpIteratorElements.begin()] = new_extension*((1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS));
            }
            
            ///Aqui es fa la cerca
            searchScheme->SearchElementsInRadiusExclusive(mrModelPart,mpIteratorElements,mRadius,mResults,mResultsDistances);

            ///Aqui ja tenim tots els resultats de tots el elements i fem el cambi de buffers de mpi a els iteradors normals de kratos
            int number_of_threads = OpenMPUtils::GetNumThreads();

            vector<unsigned int> element_partition;
            OpenMPUtils::CreatePartition(number_of_threads, mpIteratorElements.size(), element_partition);

            #pragma omp parallel for private(ResultCounter)
            for(int k=0; k<number_of_threads; k++)
            {

                IteratorType it_begin=mpIteratorElements.begin()+element_partition[k];
                IteratorType it_end=mpIteratorElements.begin()+element_partition[k+1];
                
                ResultCounter = element_partition[k];
                
                for (IteratorType particle_pointer_it = it_begin; particle_pointer_it != it_end; ++particle_pointer_it, ++ResultCounter)
                {                   
                    (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = mResults[ResultCounter].size();
                    
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
      
                    // GETTING NEW NEIGHBOURS

                    double radius = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                    int InitialNeighbourSize = (*particle_pointer_it)->GetValue(INI_NEIGHBOURS_IDS).size();
                    int OldNeighbourSize = TempIds.size();
                    unsigned int neighbour_counter = 0;

                    for (IteratorType neighbour_it = mResults[ResultCounter].begin(); neighbour_counter < mResults[ResultCounter].size(); ++neighbour_it)
                    { 
                        double initial_delta = 0.0;
                        int failure_id       = 1;
                        //int ini_failure      = 0;
                        double indentation   = 0.0;
                        
                        bool already_added = false;
                                
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
            }// OpenMP loop
            
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
            
        ModelPart& mrModelPart;
        ElementsContainerType& mpIteratorElements;
        ProcessInfo& mrCurrentProcessInfo;
      
        VectorResultElementsContainerType   mResults;
        VectorDistanceType                  mResultsDistances;
        RadiusArrayType                     mRadius;

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



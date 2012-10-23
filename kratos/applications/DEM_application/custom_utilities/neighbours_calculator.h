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
#include "containers/weak_pointer_vector.h"
#include "containers/pointer_vector.h"
#include "containers/pointer_vector_set.h"

#include "custom_utilities/discrete_particle_configure.h"
//#include "../applications/DEM_application/DEM_application.h"
//#include "../applications/DEM_application/custom_utilities/discrete_particle_configure.h"   //M: al altre hi tens lo del charlie



namespace Kratos {

    template< //pot no compilar en windows aquest tipus d'assignacio per template.
    std::size_t TDim,
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
        //**************************************************************************************************************************************************************
        // Bucket types
        typedef Bucket < TDim, Particle, ParticlePointerVector> BucketType;


        typedef BinsObjectDynamic <ConfigureType> bins; //static Bins..?


        /// Pointer definition of Neighbour_calculator
        //KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);  R: necesitu?

        virtual ~Neighbours_Calculator() {
        };


        static void Search_Ini_Neighbours(ContainerType& pElements, ProcessInfo& rCurrentProcessInfo, bool extension_option)

        {

        KRATOS_TRY

        double radius_extend = 0.0;
        if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];

      
        boost::timer kdtree_construction;

        unsigned int MaximumNumberOfResults = 100;

        ResultContainerType Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);

        bins particle_bin(pElements.begin(), pElements.end());
        boost::timer search_time;

        //**************************************************************************************************************************************************************

        ResultIteratorType results_begin;
        DistanceIteratorType result_distances_begin;

        for (IteratorType particle_pointer_it = pElements.begin();
                    particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
            
            


            Element::GeometryType& geom = (*particle_pointer_it)->GetGeometry();

            double search_radius = (1.0 + radius_extend) * geom(0)->GetSolutionStepValue(RADIUS);
            
            //find all of the new particles within the radius
            //looks which of the new particles is inside the radius around the working particle

            results_begin = Results.begin();
            result_distances_begin = ResultsDistances.begin();

            (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = particle_bin.SearchObjectsInRadiusInner(*(particle_pointer_it),
            search_radius, results_begin, result_distances_begin, MaximumNumberOfResults);

            // CLEARING AND INITIALITZING.

            vector< int > TempIds;
            TempIds.swap((*particle_pointer_it)->GetValue(NEIGHBOURS_IDS)); //afegir a getting

            (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).clear();
            (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).clear();
            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).clear();
            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).clear();
            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).clear();
            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).clear();
            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).clear();


            int n_neighbours = (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS); // the number is correct.




            (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS).resize(n_neighbours);

            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(n_neighbours);

            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(n_neighbours);

            (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(n_neighbours);

            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(n_neighbours);

            (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(n_neighbours);


            // STORING THE NEIGHBOURS AND PROPERTIES

            int neighbour_counter = 0;

            for (ResultIteratorType neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
            {

                (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);
    
                (*particle_pointer_it)->GetValue(NEIGHBOURS_IDS)[neighbour_counter] = (*neighbour_it)->Id();

                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[neighbour_counter][0] = 0.0;
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[neighbour_counter][1] = 0.0;
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[neighbour_counter][2] = 0.0;

                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[neighbour_counter][0] = 0.0;
                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[neighbour_counter][1] = 0.0;
                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[neighbour_counter][2] = 0.0;

                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[neighbour_counter] = 1;
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[neighbour_counter] = 1;
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[neighbour_counter] = 0.0;


                //HE FET QUE AQUÍ NO EM MIRO RES... I EL SET INITIAL NEIGHBOURS AGAFA AIXO I HO MODIFICA CLAR. HO PRODRIA FER AQUI

                ++neighbour_counter;

            } // for each neighbour, neighbour_it.


        }//Loop for evey particle as a base.

          KRATOS_CATCH("")
        }// Search_Neighbours




    

         static void Search_Neighbours(ContainerType& pElements, ProcessInfo& rCurrentProcessInfo, bool extension_option)

        {


            
         KRATOS_TRY

        double radius_extend = 0.0;
        if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];

        const int case_OPTION = rCurrentProcessInfo[CASE_OPTION];
        bool delta_OPTION = false;
        bool continuum_simulation_OPTION = false;

        /* as we don't accept bool variables we need the CASE_OPTION
         * this has to be improved...
         */

        switch (case_OPTION) {
            case 0:
                delta_OPTION = false;
                continuum_simulation_OPTION = false;
                break;
            case 1:
                delta_OPTION = true;
                continuum_simulation_OPTION = false;
                break;
            case 2:
                delta_OPTION = true;
                continuum_simulation_OPTION = true;
                break;
            case 3:
                delta_OPTION = false;
                continuum_simulation_OPTION = true;
                break;
            default:
                delta_OPTION = false;
                continuum_simulation_OPTION = false;
        }

        boost::timer kdtree_construction;

        unsigned int MaximumNumberOfResults = 100;

        ResultContainerType Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);

        bins particle_bin(pElements.begin(), pElements.end());
        boost::timer search_time;
        //**************************************************************************************************************************************************************
       
        ResultIteratorType results_begin;
        DistanceIteratorType result_distances_begin; 

        //loop over all of the particles in the list to perform search

        for (IteratorType particle_pointer_it = pElements.begin();
                    particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
  
            Element::GeometryType& geom = (*particle_pointer_it)->GetGeometry();

            double search_radius = (1.0 + radius_extend) * geom(0)->GetSolutionStepValue(RADIUS);
            
            //NEW EXTENDED SEARCH METHOD
            //####################################################################################################################
            double new_extension = 2.0; ///WARNING: PROVISIONALLY SET AS 2.0. SHOULD BE CALCULATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            //####################################################################################################################
        
            search_radius = search_radius*new_extension;

            //find all of the new particles within the radius
            //looks which of the new particles is inside the radius around the working particle

            results_begin = Results.begin();
            result_distances_begin = ResultsDistances.begin();
           
            (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS) = particle_bin.SearchObjectsInRadiusInner(*(particle_pointer_it),
            search_radius, results_begin, result_distances_begin, MaximumNumberOfResults);

             int n_neighbours = (*particle_pointer_it)->GetValue(POTENTIAL_NEIGHBOURS); 

            // SAVING THE OLD NEIGHBOURS, FORCES, FAILURE TYPES AND NUMBER OF NEIGHBOURS.

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

            // CLEARING AND INITIALITZING.

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

            double radius = (*particle_pointer_it)->GetGeometry()[0].GetSolutionStepValue(RADIUS);

            int InitialNeighbourSize = (*particle_pointer_it)->GetValue(INI_NEIGHBOURS_IDS).size();

            int OldNeighbourSize = TempIds.size();

            int neighbour_counter = 0;

            for (ResultIteratorType neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
            {

                double initial_delta    = 0.0;
                int failure_id          = 1;
                int ini_failure         = 0;
                double indentation      = 0.0;

                bool already_added      = false;

                for (int IniNeighbourCounter = 0; IniNeighbourCounter != InitialNeighbourSize; IniNeighbourCounter++)
                {
                    if ( static_cast<int>((*neighbour_it)->Id()) == (*particle_pointer_it)->GetValue(INI_NEIGHBOURS_IDS)[IniNeighbourCounter]) // is this and initial neighbour?
                    {

                        initial_delta = (*particle_pointer_it)->GetValue(PARTICLE_INITIAL_DELTA)[IniNeighbourCounter];
                        ini_failure = (*particle_pointer_it)->GetValue(PARTICLE_INITIAL_FAILURE_ID)[IniNeighbourCounter]; 
                        //UUUUU sha de canviar al cpp que es vagi canviant

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

                //Not old ones:

                double other_radius                 = (*neighbour_it)->GetGeometry()[0].GetSolutionStepValue(RADIUS);

                array_1d<double,3> other_to_me_vect = (*particle_pointer_it)->GetGeometry()(0)->Coordinates() - (*neighbour_it)->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                              other_to_me_vect[1] * other_to_me_vect[1] +
                                                              other_to_me_vect[2] * other_to_me_vect[2]);
                double radius_sum                   = radius + other_radius;
                indentation                  = radius_sum - distance - initial_delta;

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




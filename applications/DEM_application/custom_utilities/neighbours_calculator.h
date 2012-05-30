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

//#include "DEM_application.h"
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
        /*
                    //typedef PointerVector<Particle>                    ParticleVector;
                    typedef WeakPointerVector<Particle>                ParticleWeakVector;
                    typedef PointerVector<ParticlePointer>             ParticlePointerVector;
                    typedef typename ParticleVector::iterator          ParticleIterator;
                    typedef typename ParticleWeakVector::iterator      ParticleWeakIterator;
                    //typedef typename ParticlePointerVector::iterator   ParticlePointerIterator;

         *
         */
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

        static void Search_Neighbours(ContainerType& pElements, ProcessInfo& rCurrentProcessInfo, bool extension_option) {

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


            //            int step_data_size = model_part.GetNodalSolutionStepDataSize();
            //            double radius = 2.0 * max_radius;

            ResultContainerType Results(MaximumNumberOfResults);
            DistanceVector ResultsDistances(MaximumNumberOfResults);

            //            unsigned int bucket_size = 20;

            bins particle_bin(pElements.begin(), pElements.end());
            boost::timer search_time;
            //**************************************************************************************************************************************************************

         //   KRATOS_WATCH(pElements.size())
            ResultIteratorType results_begin; //= Results.begin();
            DistanceIteratorType result_distances_begin; //
            //loop over all of the particles in the list to perform search
            for (IteratorType particle_pointer_it = pElements.begin();
                    particle_pointer_it != pElements.end(); ++particle_pointer_it)
 {
               
                Element::GeometryType& geom = (*particle_pointer_it)->GetGeometry();
                double search_radius = (1.0 + radius_extend) * geom(0)->GetSolutionStepValue(RADIUS);

                //find all of the new particles within the radius
                //looks which of the new particles is inside the radius around the working particle

                results_begin = Results.begin();
                result_distances_begin = ResultsDistances.begin();

                ///WARNING: particle_pointer_it  funcionava també amb .base()

                ///WARNING = To be change
                (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS) = particle_bin.SearchObjectsInRadius(*(particle_pointer_it),
                search_radius, results_begin, result_distances_begin, MaximumNumberOfResults) - 1;
              


                ///WARNING:

                // This function, SearchObjectsInRadius has some problems, it finds the particle itself as a neighbours but, the result
                //stored in GetValue(NUMBER_OF_NEIGHBOURS) = is the correct one, all the neighbours but not itself.


                //M:KRATOS WATCHES... DEBBUGUING
                /*
                                               int num_neighbours = (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS);

                                               const int& particle_id = (*particle_pointer_it)->Id();
                                               KRATOS_WATCH("NOVA BUSQUEDA")
                                               KRATOS_WATCH(particle_id)
                                               KRATOS_WATCH(search_radius)
                                               KRATOS_WATCH((*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS))

                                                for (ResultIteratorType iterador_out = Results.begin(); iterador_out!= Results.begin() + (num_neighbours+1); iterador_out++)

                                                {

                                                KRATOS_WATCH((*iterador_out)->Id())
                                                KRATOS_WATCH((*iterador_out)->GetGeometry()(0)->Coordinates() )

                                                }
                 */


                // SAVING THE OLD NEIGHBOURS, FORCES, FAILURE TYPES AND NUMBER OF NEIGHBOURS.

                ParticleWeakVector TempNeighbours;


                TempNeighbours.swap((*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS));


                vector< array_1d<double, 3 > > TempContactForce;
                TempContactForce.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES));

                vector< double > TempContactFailureId; //M: temporarily defined as a double.. ha de ser un int.
                TempContactFailureId.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID));
                
                vector< double > TempContactDelta;
                TempContactDelta.swap((*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA));   //M: NO ELS GUARDO ELS NORMALS PER TEMP PERO ELS KE ET DESCUIDES LLUNY SI KE ESTAN A TEMP SEMPRE I PASSEM LES DELTES AIXI
 
                
                vector< double > TempRotateSpringFailType;
                //Vector & RealRotateSpringFailType = rE->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE);
                TempRotateSpringFailType.swap((*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE));
                //TempRotateSpringFailType.swap(RealRotateSpringFailType);


                vector< array_1d<double, 3 > > TempRotateSpringMoment;
                TempRotateSpringMoment.swap((*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT));
                
                
                /*
                Vector TempRotateSpringMoment;
                Vector & RealRotateSpringMoment = rE->GetValue(PARTICLE_ROTATE_SPRING_MOMENT);
                TempRotateSpringMoment.swap(RealRotateSpringMoment);
                */


                //M:in general we don't search here but we need it for the neigbours that the search calculator doesnt find.

                int n_neighbours = (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS); // the number is correct.

                // CLEARING AND INITIALITZING.

                (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).clear();

                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).clear();
                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).clear();

                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(n_neighbours);
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(n_neighbours);
                (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(n_neighbours);

                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(n_neighbours);
                (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(n_neighbours);

                //if( rCurrentProcessInfo[DUMMY_SWITCH] == 0)   (*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS).resize(n_neighbours);

                //KRATOS_WATCH("alohohohohoh1")

                // GETTING NEW NEIGHBOURS

                int neighbour_counter = 0;

                for (ResultIteratorType neighbour_it = Results.begin(); neighbour_counter != n_neighbours + 1; ++neighbour_it)
 {

                    if ((*particle_pointer_it)->Id() != (*neighbour_it)->Id()) { //the bins search finds the particle itself


                        (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);

                      
                        // LOOP TO EXTEND THE VECTORS AND SET A 0.0 VALUE EACH TIME

                        size_t Notemp = ((*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS)).size(); //Notemp its an index that will be amplified for every step of the loop as the NEIGHBOUR_ELEMENTS.size() is getting larger


                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp - 1] = ZeroVector(3);

                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp - 1] = 1;

                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[Notemp - 1] = 0.0;

                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[Notemp - 1] = 0.0;

                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp - 1] = ZeroVector(3);


                        // LOOP OVER THE OLD NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN EXISTING ONE AND COPYING THE OLD DATA

                        //R: FALTA PER REVISAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        int OldNeighbourCounter = 0;
                        for (ParticleWeakIterator old_neighbour = TempNeighbours.begin(); old_neighbour != TempNeighbours.end(); old_neighbour++)
                        {
                            // if ((*particle_pointer_it)->Id() != old_neighbour->Id() ) ///WARNING: NO FA FALTA!!!

                            {

                                if ((old_neighbour.base())->expired() == false) {
                                    if ((*neighbour_it)->Id() == old_neighbour->Id()) // MIQUEL COMPROBA SI EL TROBES O NO
                                    {
                             
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp-1][0] = TempContactForce[OldNeighbourCounter][0];
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp-1][1] = TempContactForce[OldNeighbourCounter][1];
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp-1][2] = TempContactForce[OldNeighbourCounter][2];

                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp-1] = TempContactFailureId[OldNeighbourCounter];

                                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp-1][0] = TempRotateSpringMoment[OldNeighbourCounter][0];
                                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp-1][1] = TempRotateSpringMoment[OldNeighbourCounter][1];
                                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp-1][2] = TempRotateSpringMoment[OldNeighbourCounter][2];

                                        (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[Notemp-1] = TempRotateSpringFailType[OldNeighbourCounter];


                                        break;
                                    } //end of its an old one??

                                } //end of expired?

                            } // end of its myself
                            OldNeighbourCounter++;
                        } //loop old neighbours


                        if (delta_OPTION) {

                            //r: FALTA REVISAR


                            // LOOP OVER THE INITIAL NEIGHBOURS FOR EVERY NEIGHBOUR TO CHECK IF IT'S AN INITIAL ONE AND THEN COPYING THE DELTA DATA
                            int InitialNeighboursCounter = 0;
     
                            if (((*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS)).size() != 0) {
                                for (ParticleWeakIterator ini_neighbour = ((*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS)).begin(); ini_neighbour != ((*particle_pointer_it)->GetValue(INITIAL_NEIGHBOUR_ELEMENTS)).end(); ini_neighbour++)
                                {
                                   
                                       if ((ini_neighbour.base())->expired() == false) {
                                   
                                            if ((*neighbour_it)->Id() == ini_neighbour->Id()) // IF IT'S AN INITIAL NEIGHBOUR //POOOOOYAN, ESTARA BÉ AIXO???? TINDRAN DIFERENTS IDS??? INI NEIGH TINDRA ID??
                                            {
                                                 (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[Notemp-1] = (*particle_pointer_it)->GetValue(PARTICLE_INITIAL_DELTA)[InitialNeighboursCounter];

                                                 break;
                                            }
                                       }
                               

                                    InitialNeighboursCounter++;

                                } // for initial neighbours
                            } //if u have some intial neigh

                        }//deltaOPTION
                        //KRATOS_WATCH("alohohohohoh4.32")
                    }//end of the: if((*particle_pointer_it)->Id() != (*neighbour_it)->Id()

                    ++neighbour_counter;

                } // for each neighbour, neighbour_it.

                //ADDING NOT FOUND NEIGHBOURS (the ones with negative identation still in tensile contact are not detected, but they are on the old neighbours list).

                ///WARNING : SHA DE REVISAR AKESTA PART ENCARA AMB UN CAS UNA MICA ESPECIAL.

                int TempNeighbourCounter = 0;

                for (ParticleWeakIterator temp_neighbour = TempNeighbours.begin(); temp_neighbour != TempNeighbours.end(); temp_neighbour++)
 {
                    if (TempContactFailureId[TempNeighbourCounter] == 0) // if they are not detached.
                    {
                        //KRATOS_WATCH("alohohohohoh5.2")
                        if ((temp_neighbour.base())->expired() == false)
                        {
                            //KRATOS_WATCH("alohohohohoh5.3")
                            if ((*particle_pointer_it)->Id() != temp_neighbour->Id()) {
                                //KRATOS_WATCH("alohohohohoh5.4")
                                bool AlreadyAdded = false; //identifying if they are already found ot not.

                                for (ParticleWeakIterator new_neighbour = (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).begin();
                                        new_neighbour != (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).end(); new_neighbour++) {
                                    //KRATOS_WATCH("alohohohohoh5.5")
                                    if (new_neighbour->Id() == (temp_neighbour)->Id()) {

                                        AlreadyAdded = true; //for the ones already found in the new search.
                                        break;
                                    }
                                }

                                if (AlreadyAdded == false) //for the ones not included!
                                {
                                          

                                    (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(TempNeighbours(TempNeighbourCounter)); //adding the not found neighbours.

                                    size_t Notemp = (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).size();
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES).resize(Notemp + 1); // adding one more space for every missing neighbour.
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][0] = TempContactForce[TempNeighbourCounter][0]; //copying properties.
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][1] = TempContactForce[TempNeighbourCounter][1];
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FORCES)[Notemp][2] = TempContactForce[TempNeighbourCounter][2];

                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID).resize(Notemp + 1);
                                    (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_FAILURE_ID)[Notemp] = TempContactFailureId[TempNeighbourCounter];

                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT).resize(Notemp + 1); // adding one more space for every missing neighbour.
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp][0] = TempRotateSpringMoment[TempNeighbourCounter][0]; //copying properties.
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp][1] = TempRotateSpringMoment[TempNeighbourCounter][1];
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[Notemp][2] = TempRotateSpringMoment[TempNeighbourCounter][2];

                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE).resize(Notemp + 1);
                                    (*particle_pointer_it)->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)[Notemp] = TempRotateSpringFailType[TempNeighbourCounter];


                                    if (delta_OPTION) {
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA).resize(Notemp + 1);
                                        (*particle_pointer_it)->GetValue(PARTICLE_CONTACT_DELTA)[Notemp] = TempContactDelta[TempNeighbourCounter];
                                    }

                                    (*particle_pointer_it)->GetValue(NUMBER_OF_NEIGHBOURS)++;

                                }
                            }// end its myself???

                        }//if not expired

                    } //if not detached
                    
                    TempNeighbourCounter++;
                }//loop over tempneigh

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



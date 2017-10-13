//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change
//something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_NEIGHBOURS_CALCULATOR_SPH )
#define  KRATOS_NEIGHBOURS_CALCULATOR_SPH

#define _OPENMPI 1

//M: we are using static bins for objects...

#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/weak_pointer_vector.h"
#include "containers/pointer_vector.h"
#include "containers/pointer_vector_set.h"
#include "meshless_application.h"

#include "meshless_application_variables.h"


#include "custom_utilities/SPH_particle_configure1.h"

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

class Neighbours_Calculator_SPH {
public:
    typedef SPHParticleConfigure < 3 > ConfigureType;

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

    virtual ~Neighbours_Calculator_SPH() {
    }

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
//        double cellsize=0.0006;
        Bins particle_bin(pIteratorElements.begin(), pIteratorElements.end());
        particle_bin.SearchObjectsInRadius(pIteratorElements.begin(),NumberOfElements,Radius,Results,ResultsDistances,NumberOfResults,MaximumNumberOfResults);

    }



    void Search_Neighbours(ModelPart& r_model_part)
    {
        KRATOS_TRY;
        ContainerType& pIteratorElements = Get_Elements(r_model_part);


        unsigned int MaximumNumberOfResults = 1000; // WATCH OUT!!
        unsigned int ResultIterator = 0;

        //**************************************************************************************************************************************************************

        unsigned int NumberOfElements = pIteratorElements.end() - pIteratorElements.begin(); // I made it unsigned

        Clean_Modelpart(r_model_part);


        if(NumberOfResults.size() != NumberOfElements)
        {
            NumberOfResults.resize(NumberOfElements);
            Results.resize(NumberOfElements);
            ResultsDistances.resize(NumberOfElements);
            Radius.resize(NumberOfElements);



            for(unsigned int i=0; i<NumberOfElements; i++ )
            {
                Results[i].resize(MaximumNumberOfResults);
                ResultsDistances[i].resize(MaximumNumberOfResults);
            }
        }


//        std::vector<std::size_t>               NumberOfResults(NumberOfElements);
//        std::vector<std::vector< PointerType > > Results(NumberOfElements, std::vector<PointerType>(MaximumNumberOfResults));
//        std::vector<std::vector<double> >      ResultsDistances(NumberOfElements, std::vector<double>(MaximumNumberOfResults));
//        std::vector<double>                    Radius(NumberOfElements);

        ///Radius vector fill
        for (IteratorType particle_pointer_it = pIteratorElements.begin(); particle_pointer_it != pIteratorElements.end(); ++particle_pointer_it)
        {

            Radius[particle_pointer_it - pIteratorElements.begin()] = (*particle_pointer_it)->GetGeometry()(0)->FastGetSolutionStepValue(SEARCH_RADIUS);


        }

//        KRATOS_WATCH(NumberOfElements);

        ///Aqui es fa la cerca
        SearchNeighbours(r_model_part,pIteratorElements,NumberOfElements,MaximumNumberOfResults,NumberOfResults,Results,ResultsDistances,Radius);





        // CHAPUZA :)

        for (IteratorType particle_pointer_it = pIteratorElements.begin(); particle_pointer_it != pIteratorElements.end(); ++particle_pointer_it,++ResultIterator)
        {

            unsigned int neighbour_counter = 0;
            for (ResultIteratorType neighbour_it = Results[ResultIterator].begin(); neighbour_counter < NumberOfResults[ResultIterator]; ++neighbour_it)
            {
                (*particle_pointer_it)->GetValue(NEIGHBOUR_ELEMENTS).push_back(*neighbour_it);
                ++neighbour_counter;

            }
        }



        KRATOS_CATCH("");
    }// Search_Neighbours

    virtual std::string Info() const {
        return "neighbour_calculator_for_SPH";
    }

    virtual void PrintInfo(std::ostream& rOStream) const {
    }

    virtual void PrintData(std::ostream& rOStream) const {
    }

protected:


private:
    std::vector<std::size_t>               NumberOfResults;
    std::vector<std::vector< PointerType > > Results;
    std::vector<std::vector<double> >      ResultsDistances;
    std::vector<double>                    Radius;


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

    Neighbours_Calculator_SPH & operator=(Neighbours_Calculator_SPH const& rOther);

}; // Class Neighbours_calculator_SPH


} // namespace Kratos.

#endif // KRATOS_NEIGHBOURS_CALCULATOR_SPH  defined



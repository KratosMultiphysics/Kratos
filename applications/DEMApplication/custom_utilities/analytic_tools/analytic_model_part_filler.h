//
// Author: Guillermo Casas (gcasas@cimne.upc.edu)
//
#ifndef ANALYTIC_MODEL_PART_FILLER_H
#define ANALYTIC_MODEL_PART_FILLER_H

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../custom_elements/analytic_spheric_particle.h"
#include "../create_and_destroy.h"


/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) AnalyticModelPartFiller
{

public:   

KRATOS_CLASS_POINTER_DEFINITION(AnalyticModelPartFiller);

typedef ModelPart::ElementsContainerType ElementsArrayType;
typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;

/// Default constructor
AnalyticModelPartFiller(){}
/// Destructor

virtual ~AnalyticModelPartFiller(){}

// This algorithm implements the Fisher–Yates shuffle to generate a random subset of elements of a vector class
template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random){
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

void GetRandomSample(std::vector<int>& random_positions_vector_to_fill,
                     const int n_positions,
                     const int n_random_positions);

void FillAnalyticModelPartGivenFractionOfParticlesToTransform(const double fraction_of_particles_to_convert,
                                                              ModelPart& spheres_model_part,
                                                              ParticleCreatorDestructor particle_creator_destructor,
                                                              std::string analytic_sub_model_part_name = "");


/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

/// Assignment operator
AnalyticModelPartFiller & operator=(AnalyticModelPartFiller const& rOther);

}; // Class AnalyticModelPartFiller

} // namespace Kratos.

#endif // ANALYTIC_MODEL_PART_FILLER_H

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


/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "boost/python/list.hpp"

namespace Kratos
{
class AnalyticModelPartFiller {

public:   

KRATOS_CLASS_POINTER_DEFINITION(AnalyticModelPartFiller);

/// Default constructor

AnalyticModelPartFiller(){}

/// Destructor

virtual ~AnalyticModelPartFiller(){}

public:


// This algorithm implements the Fisher–Yates shuffle to generate a random subset of elements of a vector class
template<class bidiiter>
bidiiter random_unique(bidiiter begin, bidiiter end, size_t num_random) {
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


typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;

void GetRandomSample(std::vector<int>& random_positions_vector_to_fill, const int n_positions, const int n_random_positions, const int n_iterations = 100)
{
    assert (n_random_positions <= n_positions);
    random_positions_vector_to_fill.resize(n_positions);
    std::iota (std::begin(random_positions_vector_to_fill), std::end(random_positions_vector_to_fill), 0); // Fill with 0, 1, ..., n_elements - 1.
    for (int i = 0; i < n_iterations; ++i){
      random_unique(random_positions_vector_to_fill.begin(), random_positions_vector_to_fill.end(), n_random_positions);
    }
    random_positions_vector_to_fill.resize(n_random_positions);
}


void FillAnalyticModelPartGivenFractionOfParticlesToTransform(const double fraction_of_particles_to_convert, ModelPart& spheres_model_part, ModelPart& analytic_spheres_model_part)
{
    const int n_elements = spheres_model_part.NumberOfElements();
    const int n_analytic_particles_needed =  int(fraction_of_particles_to_convert * n_elements);
    std::vector<int> random_positions;
    GetRandomSample(random_positions, n_elements, n_analytic_particles_needed);

    ElementsIteratorType i_begin = spheres_model_part.ElementsBegin();
    ElementsIteratorType i_elem = i_begin;
    AnalyticSphericParticle analytic_sample_element;

    for (int i = 0; i < n_analytic_particles_needed; ++i){
        i_elem = i_begin + random_positions[i];
        analytic_spheres_model_part.AddElement(analytic_sample_element.Create(i_elem->Id(), i_elem->GetGeometry(), i_elem->pGetProperties()));
    }
}


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

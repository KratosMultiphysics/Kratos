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
class AnalyticModelPartFiller
{

public:   

KRATOS_CLASS_POINTER_DEFINITION(AnalyticModelPartFiller);

typedef ModelPart::ElementsContainerType ElementsArrayType;

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


void FillAnalyticModelPartGivenFractionOfParticlesToTransform(const double fraction_of_particles_to_convert,
                                                              ModelPart& spheres_model_part,
                                                              ParticleCreatorDestructor particle_creator_destructor,
                                                              std::string analytic_sub_model_part_name = "")
{
    const int n_elements = spheres_model_part.NumberOfElements();
    const int n_analytic_particles_needed =  int(fraction_of_particles_to_convert * n_elements);
    std::vector<int> random_positions;
    GetRandomSample(random_positions, n_elements, n_analytic_particles_needed);
    if (analytic_sub_model_part_name == ""){
        std::string prefix("Analytic");
        analytic_sub_model_part_name = prefix + spheres_model_part.Name();
    }

    ModelPart& analytic_spheres_model_part = spheres_model_part.GetSubModelPart(analytic_sub_model_part_name);
    ElementsIteratorType i_begin = spheres_model_part.ElementsBegin();
    //ElementsIteratorType i_elem = i_begin;
    const Element& sample_element = KratosComponents<Element>::Get("AnalyticSphericParticle3D");
    ElementsArrayType elements_to_add;

    for (int i = 0; i < n_analytic_particles_needed; ++i){
        ElementsIteratorType i_elem = i_begin + random_positions[i];
        Geometry<Node<3> >::PointsArrayType nodelist;
        nodelist.push_back(i_elem->GetGeometry().pGetPoint(0));
        Element::Pointer p_elem = sample_element.Create(i_elem->Id(), nodelist, i_elem->pGetProperties());

        AnalyticSphericParticle* analytic_sample_element = dynamic_cast<AnalyticSphericParticle*> (p_elem.get());
        SphericParticle* regular_sample_element = dynamic_cast<SphericParticle*> ((*(i_elem.base())).get());

        analytic_sample_element->SetFastProperties(regular_sample_element->GetFastProperties());
        analytic_sample_element->SetRadius(nodelist[0].FastGetSolutionStepValue(RADIUS));
        analytic_sample_element->SetSearchRadius(nodelist[0].FastGetSolutionStepValue(RADIUS));
        analytic_sample_element->SetSearchRadiusWithFem(nodelist[0].FastGetSolutionStepValue(RADIUS));
        analytic_sample_element->Set(DEMFlags::HAS_ROLLING_FRICTION, false);
        analytic_sample_element->Set(DEMFlags::BELONGS_TO_A_CLUSTER, false);
        analytic_sample_element->CreateDiscontinuumConstitutiveLaws(spheres_model_part.GetProcessInfo());

        for (int i_neigh = 0; i < int(regular_sample_element->mNeighbourElements.size()); ++i_neigh){
            analytic_sample_element->mNeighbourElements.push_back(regular_sample_element->mNeighbourElements[i]);
            analytic_sample_element->mNeighbourElasticContactForces.push_back(regular_sample_element->mNeighbourElasticContactForces[i]);
            analytic_sample_element->mNeighbourElasticExtraContactForces.push_back(regular_sample_element->mNeighbourElasticExtraContactForces[i]);
        }
        analytic_sample_element->Initialize(spheres_model_part.GetProcessInfo());

        elements_to_add.push_back(p_elem);
        i_elem->Set(TO_ERASE);
        spheres_model_part.RemoveElements(TO_ERASE);
    }

    analytic_spheres_model_part.AddElements(elements_to_add.begin(), elements_to_add.end());
    KRATOS_WATCH(spheres_model_part)
    KRATOS_WATCH(analytic_spheres_model_part)

    //particle_creator_destructor.DestroyParticles(spheres_model_part);
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

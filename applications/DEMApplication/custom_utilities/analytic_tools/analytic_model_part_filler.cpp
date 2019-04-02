//   $Main author: Guillermo Casas
//

// Project includes

// System includes
#include <numeric>
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_model_part_filler.h"

namespace Kratos
{

void AnalyticModelPartFiller::GetRandomSample(std::vector<int>& random_positions_vector_to_fill,
                                              const int n_positions,
                                              const int n_random_positions)
{
    assert (n_random_positions <= n_positions);
    random_positions_vector_to_fill.resize(n_positions);
    std::iota (std::begin(random_positions_vector_to_fill), std::end(random_positions_vector_to_fill), 0); // Fill with 0, 1, ..., n_elements - 1.
    random_unique(random_positions_vector_to_fill.begin(), random_positions_vector_to_fill.end(), n_random_positions);
    random_positions_vector_to_fill.resize(n_random_positions);
}

void AnalyticModelPartFiller::FillAnalyticModelPartGivenFractionOfParticlesToTransform(const double fraction_of_particles_to_convert,
                                                              ModelPart& spheres_model_part,
                                                              ParticleCreatorDestructor particle_creator_destructor,
                                                              std::string analytic_sub_model_part_name)
{
    const int n_elements = spheres_model_part.NumberOfElements();
    const int n_analytic_particles_needed =  int(fraction_of_particles_to_convert * n_elements);
    std::vector<int> random_positions;
    GetRandomSample(random_positions, n_elements, n_analytic_particles_needed);
    std::string prefix("Analytic");
    if (analytic_sub_model_part_name == ""){
        analytic_sub_model_part_name = prefix + spheres_model_part.Name();
    }

    ModelPart& analytic_spheres_model_part = spheres_model_part.GetSubModelPart(analytic_sub_model_part_name);
    ElementsArrayType elements_to_add;
    ElementsIteratorType i_begin = spheres_model_part.ElementsBegin();
    ElementsIteratorType i_elem = i_begin;
    std::string element_type_name = i_begin->Info() + "3D";
    const Element& sample_element = KratosComponents<Element>::Get(prefix + element_type_name);

    // Producing replacements with the same Ids
    for (int i = 0; i < n_analytic_particles_needed; ++i){
        i_elem = i_begin + random_positions[i];
        Geometry<Node<3> >::PointsArrayType nodelist;
        nodelist.push_back(i_elem->GetGeometry().pGetPoint(0));
        Element::Pointer p_elem = particle_creator_destructor.GetAnalyticReplacement(sample_element,
                                                                                     nodelist,
                                                                                     *(i_elem.base()),
                                                                                     spheres_model_part);

        elements_to_add.push_back(p_elem);
        i_elem->Set(TO_ERASE);
    }

    // Eliminating replaced particles
    particle_creator_destructor.DestroyParticleElements(spheres_model_part, TO_ERASE);

    // Adding replacements to model_part without repeating Ids, now that the old ones have been eliminated
    analytic_spheres_model_part.AddElements(elements_to_add.begin(), elements_to_add.end());

}

/// Turn back information as a string.
std::string AnalyticModelPartFiller::Info() const {
        return "AnalyticModelPartFiller";
}

/// Print information about this object.
void AnalyticModelPartFiller::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticModelPartFiller::PrintData(std::ostream& rOStream) const {}

} // namespace Kratos

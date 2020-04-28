//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_PARTICLE_GENERATOR_UTILITY
#define KRATOS_MPM_PARTICLE_GENERATOR_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "particle_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "utilities/binbased_fast_point_locator.h"

// Quadrature point imports
#include "geometries/quadrature_point_geometry.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Node<3> NodeType;
    typedef Geometry<Node<3>> GeometryType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * @brief Function that return matrix of shape function value for 16 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    Matrix MP16ShapeFunctions();


    /**
     * @brief Function that return matrix of shape function value for 33 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    Matrix MP33ShapeFunctions();


    /**
     * @brief Construct material points or particles from given initial mesh
     * @details Generating particles using a designated shape functions
     */
    template<std::size_t TDimension>
    void GenerateMaterialPointElement(ModelPart& rBackgroundGridModelPart,
        ModelPart& rInitialModelPart,
        ModelPart& rMPMModelPart,
        bool IsAxisSymmetry,
        bool IsMixedFormulation)
    {
        // Initialize zero the variables needed
        array_1d<double, 3> xg = ZeroVector(3);
        array_1d<double, 3> mp_displacement = ZeroVector(3);
        array_1d<double, 3> mp_velocity = ZeroVector(3);
        array_1d<double, 3> mp_acceleration = ZeroVector(3);
        array_1d<double, 3> mp_volume_acceleration = ZeroVector(3);

        std::vector<Vector> mp_cauchy_stress_vector = { ZeroVector(6) };
        std::vector<Vector> mp_almansi_strain_vector = { ZeroVector(6) };
        double mp_pressure = 0.0;

        double mp_mass;
        double mp_volume;

        // Determine element index: This convention is done in order for the purpose of visualization in GiD
        const SizeType number_elements = rBackgroundGridModelPart.NumberOfElements() + rInitialModelPart.NumberOfElements();
        const SizeType number_nodes = rBackgroundGridModelPart.NumberOfNodes();
        IndexType last_element_id = (number_nodes > number_elements) ? (number_nodes + 1) : (number_elements + 1);


        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        SearchStructure.UpdateSearchDatabase();

        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);

        // Loop over the submodelpart of rInitialModelPart
        for (ModelPart::SubModelPartIterator submodelpart_it = rInitialModelPart.SubModelPartsBegin();
            submodelpart_it != rInitialModelPart.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart& submodelpart = *submodelpart_it;
            std::string submodelpart_name = submodelpart.Name();

            rMPMModelPart.CreateSubModelPart(submodelpart_name);

            // Loop over the element of submodelpart and generate mpm element to be appended to the rMPMModelPart
            for (ModelPart::ElementIterator i = submodelpart.ElementsBegin();
                i != submodelpart.ElementsEnd(); i++)
            {
                if (i->IsDefined(ACTIVE))
                {
                    Properties::Pointer p_properties = i->pGetProperties();

                    // Check number of particles per element to be created
                    const SizeType particles_per_element = (p_properties->Has(PARTICLES_PER_ELEMENT))
                        ? i->GetProperties()[PARTICLES_PER_ELEMENT]
                        : 1;
                    KRATOS_WARNING_IF("MPMParticleGeneratorUtility", !p_properties->Has(PARTICLES_PER_ELEMENT))
                        << "PARTICLES_PER_ELEMENT is not specified in Properties. 1 particle per element is assumed." << std::endl;

                    // Get geometry and dimension of the background grid
                    const SizeType working_space_dimension = rBackgroundGridModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
                    const SizeType num_nodes = rBackgroundGridModelPart.ElementsBegin()->GetGeometry().size();
                    const SizeType domain_size = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];

                    // Current element's geometry
                    const Geometry<Node<3>>& r_geometry = i->GetGeometry();

                    IntegrationPointsArrayType integration_points;
                    MPMParticleGeneratorUtility::GetIntegrationPoints(integration_points, r_geometry, working_space_dimension, num_nodes, particles_per_element);

                    // Number of MP per elements
                    const SizeType integration_point_per_elements = integration_points.size();

                    const double density = i->GetProperties()[DENSITY];
                    //// Evaluation of element area/volume
                    const double area = r_geometry.Area();
                    if (domain_size == 2 && i->GetProperties().Has(THICKNESS)) {
                        const double thickness = i->GetProperties()[THICKNESS];
                        mp_mass = area * thickness * density / integration_point_per_elements;
                    }
                    else {
                        mp_mass = area * density / integration_point_per_elements;
                    }
                    mp_volume = area / integration_point_per_elements;

                    // Loop over the material points that fall in each grid element
                    IndexType new_element_id = 0;
                    for (IndexType PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++)
                    {
                        array_1d<double, 3> coords;
                        r_geometry.GlobalCoordinates(coords, integration_points[PointNumber]);

                        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                        Element::Pointer pelem;
                        Vector N;

                        // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                        bool is_found = SearchStructure.FindPointOnMesh(coords, N, pelem, result_begin, 10, 1e-3);

                        if (is_found == true) {
                            pelem->Set(ACTIVE);
                            auto p_new_geometry = CreateQuadraturePointsUtility<NodeType>::CreateFromCoordinates(
                                pelem->pGetGeometry(), coords,
                                integration_points[PointNumber].Weight());

                            const Element& new_element = KratosComponents<Element>::Get("UpdatedLagrangianElement");

                            // Create new material point element
                            new_element_id = last_element_id + PointNumber;
                            Element::Pointer p_element = new_element.Create(
                                new_element_id, p_new_geometry, p_properties);

                            const ProcessInfo process_info = ProcessInfo();
                            std::vector<double> mp_mass_vector(1);
                            mp_mass_vector[0] = mp_mass;
                            std::vector<double> mp_volume_vector(1);
                            mp_volume_vector[0] = mp_volume;

                            // Setting particle element's initial condition
                            p_element->SetValuesOnIntegrationPoints(MP_MASS, mp_mass_vector, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_VOLUME, mp_volume_vector, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_COORD, { coords }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_DISPLACEMENT, { mp_displacement }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_VELOCITY, { mp_velocity }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_ACCELERATION, { mp_acceleration }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, { mp_volume_acceleration }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, { mp_cauchy_stress_vector }, process_info);
                            p_element->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, { mp_almansi_strain_vector }, process_info);

                            if (IsMixedFormulation) {
                                std::vector<double> mp_pressure = { 0.0 };

                                p_element->SetValuesOnIntegrationPoints(MP_PRESSURE, mp_pressure, process_info);
                            }

                            // Add the MP Element to the model part
                            rMPMModelPart.GetSubModelPart(submodelpart_name).AddElement(p_element);
                        }
                        else {
                            KRATOS_INFO("MPMSearchElementUtility")
                                << "WARNING: Search backgroundgrid for material point at " << coords
                                << " failed. No MPM element created." << std::endl;
                        }
                    }

                    last_element_id += integration_point_per_elements;
                }
            }
        }
    }

    /**
     * @brief Function to Initiate material point condition.
     * @details Generating particle condition using a designated shape functions
     */
    void GenerateMaterialPointCondition(ModelPart& rBackgroundGridModelPart,
                                        ModelPart& rInitialModelPart,
                                        ModelPart& rMPMModelPart);

    void GetIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const GeometryType& rGeometry,
        SizeType WorkingSpaceDimension,
        SizeType NumNodesPerGeometry,
        SizeType ParticlesPerElement);

}; // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos

#endif // KRATOS_MPM_PARTICLE_GENERATOR_UTILITY



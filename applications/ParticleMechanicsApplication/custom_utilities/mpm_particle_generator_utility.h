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

        Vector mp_cauchy_stress_vector = ZeroVector(6);
        Vector mp_almansi_strain_vector = ZeroVector(6);
        double mp_pressure = 0.0;

        double mp_mass;
        double mp_volume;

        // Determine element index: This convention is done in order for the purpose of visualization in GiD
        const unsigned int number_elements = rBackgroundGridModelPart.NumberOfElements() + rInitialModelPart.NumberOfElements();
        const unsigned int number_nodes = rBackgroundGridModelPart.NumberOfNodes();
        unsigned int last_element_id = (number_nodes > number_elements) ? (number_nodes + 1) : (number_elements + 1);


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
                    unsigned int particles_per_element;
                    if (i->GetProperties().Has(PARTICLES_PER_ELEMENT)) {
                        particles_per_element = i->GetProperties()[PARTICLES_PER_ELEMENT];
                    }
                    else {
                        std::string warning_msg = "PARTICLES_PER_ELEMENT is not specified in Properties, ";
                        warning_msg += "1 Particle per element is assumed.";
                        KRATOS_WARNING("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                        particles_per_element = 1;
                    }

                    // Get geometry and dimension of the background grid
                    const GeometryData::KratosGeometryType background_geo_type = rBackgroundGridModelPart.ElementsBegin()->GetGeometry().GetGeometryType();
                    const std::size_t domain_size = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];

                    const Geometry< Node < 3 > >& r_geometry = i->GetGeometry(); // current element's geometry

                    GeometryData::IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_1;
                    switch (particles_per_element)
                    {
                    case 3:
                        this_integration_method = GeometryData::GI_GAUSS_2;
                        break;
                    case 6:
                        this_integration_method = GeometryData::GI_GAUSS_3;
                        break;
                    case 12:
                        this_integration_method = GeometryData::GI_GAUSS_4;
                        break;
                    default:
                        std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                        warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                        warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                        warning_msg += "The default number of particle: 3 is currently assumed.";
                        KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                        break;
                    }

                    GeometryType::IntegrationPointsArrayType integration_points = r_geometry.IntegrationPoints(this_integration_method);

                    // Number of MP per elements
                    const unsigned int integration_point_per_elements = integration_points.size();

                    //CAN WE DO THIS ON THE ELEMENT?

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
                    unsigned int new_element_id = 0;
                    for (unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++)
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
                                pelem->pGetGeometry(),
                                coords,
                                integration_points[PointNumber].Weight());

                            //for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            //    r_geometry[j].Set(ACTIVE);

                            const Element& new_element = KratosComponents<Element>::Get("UpdatedLagrangianElement");

                            // Create new material point element
                            new_element_id = last_element_id + PointNumber;
                            Element::Pointer p_element = new_element.Create(
                                new_element_id,
                                p_new_geometry,
                                p_properties);

                            // Setting particle element's initial condition
                            p_element->SetValue(MP_MASS, mp_mass);
                            p_element->SetValue(MP_VOLUME, mp_volume);
                            p_element->SetValue(MP_COORD, coords);
                            p_element->SetValue(MP_DISPLACEMENT, mp_displacement);
                            p_element->SetValue(MP_VELOCITY, mp_velocity);
                            p_element->SetValue(MP_ACCELERATION, mp_acceleration);
                            p_element->SetValue(MP_VOLUME_ACCELERATION, mp_volume_acceleration);
                            p_element->SetValue(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress_vector);
                            p_element->SetValue(MP_ALMANSI_STRAIN_VECTOR, mp_almansi_strain_vector);

                            if (IsMixedFormulation)
                            {
                                p_element->SetValue(MP_PRESSURE, mp_pressure);
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
    void GenerateMaterialPointCondition(    ModelPart& rBackgroundGridModelPart,
                                            ModelPart& rInitialModelPart,
                                            ModelPart& rMPMModelPart);

}; // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos

#endif // KRATOS_MPM_PARTICLE_GENERATOR_UTILITY



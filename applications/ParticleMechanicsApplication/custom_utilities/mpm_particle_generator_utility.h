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


namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

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
    template<SizeType TDimension>
    void GenerateMaterialPointElement(  ModelPart& rBackgroundGridModelPart,
                                        ModelPart& rInitialModelPart,
                                        ModelPart& rMPMModelPart,
                                        bool IsMixedFormulation = false) {
        const bool IsAxisSymmetry = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_AXISYMMETRIC))
            ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_AXISYMMETRIC)
            : false;

        // Initialize zero the variables needed
        std::vector<array_1d<double, 3>> xg = { ZeroVector(3) };
        std::vector<array_1d<double, 3>> mp_displacement = { ZeroVector(3) };
        std::vector<array_1d<double, 3>> mp_velocity = { ZeroVector(3) };
        std::vector<array_1d<double, 3>> mp_acceleration = { ZeroVector(3) };
        std::vector<array_1d<double, 3>> mp_volume_acceleration = { ZeroVector(3) };

        std::vector<Vector> mp_cauchy_stress_vector = { ZeroVector(6) };
        std::vector<Vector> mp_almansi_strain_vector = { ZeroVector(6) };
        std::vector<double> mp_pressure = { 0.0 };

        std::vector<double> mp_mass(1);
        std::vector<double> mp_volume(1);

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
                    Properties::Pointer properties = i->pGetProperties();
                    const double density = i->GetProperties()[DENSITY];

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
                    const GeometryData::KratosGeometryType geo_type = r_geometry.GetGeometryType();
                    Matrix shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                    if (geo_type == GeometryData::Kratos_Tetrahedra3D4 || geo_type == GeometryData::Kratos_Triangle2D3)
                    {
                        switch (particles_per_element)
                        {
                        case 1:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                            break;
                        case 3:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            break;
                        case 6:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                            break;
                        case 12:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_5);
                            break;
                        case 16:
                            if (domain_size == 2) {
                                shape_functions_values = MP16ShapeFunctions();
                                break;
                            }
                        case 33:
                            if (domain_size == 2) {
                                shape_functions_values = MP33ShapeFunctions();
                                break;
                            }
                        default:
                            std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                            warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                            warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                            warning_msg += "The default number of particle: 3 is currently assumed.";
                            KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                            break;
                        }
                    }
                    else if (geo_type == GeometryData::Kratos_Hexahedra3D8 || geo_type == GeometryData::Kratos_Quadrilateral2D4)
                    {
                        switch (particles_per_element)
                        {
                        case 1:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                            break;
                        case 4:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                            break;
                        case 9:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_3);
                            break;
                        case 16:
                            shape_functions_values = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                            break;
                        default:
                            std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                            warning_msg += " is not available for Quadrilateral" + std::to_string(domain_size) + "D.\n";
                            warning_msg += "Available options are: 1, 4, 9, 16.\n";
                            warning_msg += "The default number of particle: 4 is currently assumed.";
                            KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                            break;
                        }
                    }

                    // Set element type
                    std::string element_type_name = "UpdatedLagrangian";
                    if (IsMixedFormulation) {
                        if (background_geo_type == GeometryData::Kratos_Triangle2D3) element_type_name = "UpdatedLagrangianUP";
                        else KRATOS_ERROR << "Element for mixed U-P formulation is only implemented for 2D Triangle Elements." << std::endl;
                    }
                    else if (IsAxisSymmetry && domain_size == 3) KRATOS_ERROR << "Axisymmetric elements must be used in a 2D domain. You specified a 3D domain." << std::endl;

                    // Get new element
                    const Element& new_element = KratosComponents<Element>::Get(element_type_name);

                    // Number of MP per elements
                    const unsigned int integration_point_per_elements = shape_functions_values.size1();

                    // Evaluation of element area/volume
                    // TODO this should be corrected, in the general case the mass is not evenly distributed.
                    const double area = r_geometry.Area();
                    if (domain_size == 2 && i->GetProperties().Has(THICKNESS)) {
                        const double thickness = i->GetProperties()[THICKNESS];
                        mp_mass[0] = area * thickness * density / integration_point_per_elements;
                    }
                    else {
                        mp_mass[0] = area * density / integration_point_per_elements;
                    }
                    mp_volume[0] = area / integration_point_per_elements;

                    // Loop over the material points that fall in each grid element
                    unsigned int new_element_id = 0;
                    for (unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++)
                    {
                        std::vector<double> MP_density = { density };

                        xg[0].clear();

                        // Loop over the nodes of the grid element
                        for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++)
                        {
                            for (unsigned int j = 0; j < r_geometry.size(); j++)
                            {
                                xg[0][dimension] = xg[0][dimension] + shape_functions_values(PointNumber, j) * r_geometry[j].Coordinates()[dimension];
                            }
                        }
                        if (IsAxisSymmetry)  mp_mass[0] = mp_volume[0] * 2.0 * Globals::Pi * xg[0][0] * density;

                        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                        Element::Pointer pelem;
                        Vector N;

                        // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                        bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin);
                        if (!is_found) KRATOS_WARNING("MPM particle generator utility") << "::search failed." << std::endl;

                        pelem->Set(ACTIVE);
                        auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                            pelem->pGetGeometry(), xg[0],
                            mp_volume[0]);

                        // Create new material point element
                        new_element_id = last_element_id + PointNumber;
                        Element::Pointer p_element = new_element.Create(
                            new_element_id, p_new_geometry, properties);

                        const ProcessInfo process_info = ProcessInfo();

                        // Setting particle element's initial condition
                        p_element->SetValuesOnIntegrationPoints(MP_DENSITY, MP_density, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_MASS, mp_mass, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_VOLUME, mp_volume, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_COORD, xg, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_DISPLACEMENT, mp_displacement, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_VELOCITY, mp_velocity, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_ACCELERATION, mp_acceleration, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, mp_volume_acceleration, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress_vector, process_info);
                        p_element->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, mp_almansi_strain_vector, process_info);

                        if (IsMixedFormulation)
                        {
                            p_element->SetValuesOnIntegrationPoints(MP_PRESSURE, mp_pressure, process_info);
                        }

                        // Add the MP Element to the model part
                        rMPMModelPart.GetSubModelPart(submodelpart_name).AddElement(p_element);
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



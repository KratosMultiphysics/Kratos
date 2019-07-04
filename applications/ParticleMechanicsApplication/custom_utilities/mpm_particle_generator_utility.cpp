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

// System includes

// External includes

// Project includes
#include "custom_utilities/mpm_particle_generator_utility.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "particle_mechanics_application_variables.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"


namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

    void GenerateMaterialPointElement(  ModelPart& rBackgroundGridModelPart,
                                        ModelPart& rInitialModelPart,
                                        ModelPart& rMPMModelPart,
                                        bool IsAxisSymmetry,
                                        bool IsMixedFormulation)
    {
        // Initialize zero the variables needed
        array_1d<double,3> xg = ZeroVector(3);
        array_1d<double,3> mp_displacement = ZeroVector(3);
        array_1d<double,3> mp_velocity = ZeroVector(3);
        array_1d<double,3> mp_acceleration = ZeroVector(3);
        array_1d<double,3> mp_volume_acceleration = ZeroVector(3);

        Vector mp_cauchy_stress_vector = ZeroVector(6);
        Vector mp_almansi_strain_vector = ZeroVector(6);
        double mp_pressure = 0.0;

        double mp_mass;
        double mp_volume;

        // Determine element index: This convention is done in order for the purpose of visualization in GiD
        const unsigned int number_elements = rBackgroundGridModelPart.NumberOfElements() + rInitialModelPart.NumberOfElements();
        const unsigned int number_nodes = rBackgroundGridModelPart.NumberOfNodes();
        unsigned int last_element_id = (number_nodes > number_elements) ? (number_nodes + 1) : (number_elements+1);

        // Loop over the submodelpart of rInitialModelPart
        for (ModelPart::SubModelPartIterator submodelpart_it = rInitialModelPart.SubModelPartsBegin();
                submodelpart_it != rInitialModelPart.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart&  submodelpart = *submodelpart_it;
            std::string submodelpart_name = submodelpart.Name();

            rMPMModelPart.CreateSubModelPart(submodelpart_name);

            // Loop over the element of submodelpart and generate mpm element to be appended to the rMPMModelPart
            for (ModelPart::ElementIterator i = submodelpart.ElementsBegin();
                    i != submodelpart.ElementsEnd(); i++)
            {
                if(i->IsDefined(ACTIVE))
                {
                    Properties::Pointer properties = i->pGetProperties();
                    const int material_id = i->GetProperties().Id();
                    const double density  = i->GetProperties()[DENSITY];

                    // Check number of particles per element to be created
                    unsigned int particles_per_element;
                    if (i->GetProperties().Has( PARTICLES_PER_ELEMENT )){
                        particles_per_element = i->GetProperties()[PARTICLES_PER_ELEMENT];
                    }
                    else{
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
                    Matrix shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                    if (geo_type == GeometryData::Kratos_Tetrahedra3D4  || geo_type == GeometryData::Kratos_Triangle2D3)
                    {
                        switch (particles_per_element)
                        {
                            case 1:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                break;
                            case 3:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                break;
                            case 6:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                break;
                            case 12:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                break;
                            case 16:
                                if (domain_size==2){
                                    shape_functions_values = MP16ShapeFunctions();
                                    break;
                                }
                            case 33:
                                if (domain_size==2) {
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
                    else if(geo_type == GeometryData::Kratos_Hexahedra3D8  || geo_type == GeometryData::Kratos_Quadrilateral2D4)
                    {
                        switch (particles_per_element)
                        {
                            case 1:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                break;
                            case 4:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                break;
                            case 9:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                break;
                            case 16:
                                shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
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

                    // Check element type
                    std::string element_type_name;
                    if (domain_size==2){
                        if (background_geo_type == GeometryData::Kratos_Triangle2D3){
                            if (IsMixedFormulation)
                                element_type_name = "UpdatedLagrangianUP2D3N";
                            else{
                                if (IsAxisSymmetry)
                                    element_type_name = "UpdatedLagrangianAxisymmetry2D3N";
                                else
                                    element_type_name = "UpdatedLagrangian2D3N";
                            }
                        }
                        else if (background_geo_type == GeometryData::Kratos_Quadrilateral2D4){
                            if (IsMixedFormulation)
                                KRATOS_ERROR << "Element for mixed U-P formulation in 2D for Quadrilateral Element is not yet implemented." << std::endl;
                            else{
                                if (IsAxisSymmetry)
                                    element_type_name = "UpdatedLagrangianAxisymmetry2D4N";
                                else
                                    element_type_name = "UpdatedLagrangian2D4N";
                            }
                        }
                    }
                    else if (domain_size==3){
                        if (background_geo_type == GeometryData::Kratos_Tetrahedra3D4){
                            if (IsMixedFormulation)
                                KRATOS_ERROR << "Element for mixed U-P formulation in 3D for Tetrahedral Element is not yet implemented." << std::endl;
                            else
                                element_type_name = "UpdatedLagrangian3D4N";
                        }
                        else if (background_geo_type == GeometryData::Kratos_Hexahedra3D8){
                            if (IsMixedFormulation)
                                KRATOS_ERROR << "Element for mixed U-P formulation in 3D for Hexahedral Element is not yet implemented." << std::endl;
                            else
                                element_type_name = "UpdatedLagrangian3D8N";
                        }
                    }

                    // Get new element
                    const Element& new_element = KratosComponents<Element>::Get(element_type_name);

                    // Number of MP per elements
                    const unsigned int integration_point_per_elements = shape_functions_values.size1();

                    // Evaluation of element area/volume
                    const double area = r_geometry.Area();
                    if(domain_size == 2 && i->GetProperties().Has( THICKNESS )){
                        const double thickness = i->GetProperties()[THICKNESS];
                        mp_mass = area * thickness * density / integration_point_per_elements;
                    }
                    else {
                        mp_mass = area * density / integration_point_per_elements;
                    }
                    mp_volume = area / integration_point_per_elements;

                    // Loop over the material points that fall in each grid element
                    unsigned int new_element_id = 0;
                    for ( unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++ )
                    {
                        // Create new material point element
                        new_element_id = last_element_id + PointNumber;
                        Element::Pointer p_element = new_element.Create(new_element_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                        const double MP_density  = density;
                        const int MP_material_id = material_id;

                        xg.clear();

                        // Loop over the nodes of the grid element
                        for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++)
                        {
                            for ( unsigned int j = 0; j < r_geometry.size(); j ++)
                            {
                                xg[dimension] = xg[dimension] + shape_functions_values(PointNumber, j) * r_geometry[j].Coordinates()[dimension];
                            }
                        }

                        // Setting particle element's initial condition
                        p_element->SetValue(MP_MATERIAL_ID, MP_material_id);
                        p_element->SetValue(MP_DENSITY, MP_density);
                        p_element->SetValue(MP_MASS, mp_mass);
                        p_element->SetValue(MP_VOLUME, mp_volume);
                        p_element->SetValue(MP_COORD, xg);
                        p_element->SetValue(MP_DISPLACEMENT, mp_displacement);
                        p_element->SetValue(MP_VELOCITY, mp_velocity);
                        p_element->SetValue(MP_ACCELERATION, mp_acceleration);
                        p_element->SetValue(MP_VOLUME_ACCELERATION, mp_volume_acceleration);
                        p_element->SetValue(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress_vector);
                        p_element->SetValue(MP_ALMANSI_STRAIN_VECTOR, mp_almansi_strain_vector);

                        if(IsMixedFormulation)
                        {
                            p_element->SetValue(MP_PRESSURE, mp_pressure);
                        }

                        // Add the MP Element to the model part
                        rMPMModelPart.GetSubModelPart(submodelpart_name).AddElement(p_element);
                    }

                    last_element_id += integration_point_per_elements;

                }
            }
        }
    }

/***********************************************************************************/
/***********************************************************************************/

    void GenerateMaterialPointCondition(    ModelPart& rBackgroundGridModelPart,
                                            ModelPart& rInitialModelPart,
                                            ModelPart& rMPMModelPart)
    {
        // Initialize zero the variables needed
        array_1d<double,3> mpc_xg = ZeroVector(3);
        array_1d<double,3> mpc_normal = ZeroVector(3);
        array_1d<double,3> mpc_displacement = ZeroVector(3);
        array_1d<double,3> mpc_imposed_displacement = ZeroVector(3);
        array_1d<double,3> mpc_velocity = ZeroVector(3);
        array_1d<double,3> mpc_imposed_velocity = ZeroVector(3);
        array_1d<double,3> mpc_acceleration = ZeroVector(3);
        array_1d<double,3> mpc_imposed_acceleration = ZeroVector(3);
        array_1d<double,3> mpc_contact_force = ZeroVector(3);
        array_1d<double, 3 > point_load = ZeroVector(3);

        double mpc_area = 0.0;
        double mpc_penalty_factor = 0.0;

        // Determine condition index: This convention is done in order for the purpose of visualization in GiD
        const unsigned int number_conditions = rBackgroundGridModelPart.NumberOfConditions();
        const unsigned int number_elements = rBackgroundGridModelPart.NumberOfElements() + rInitialModelPart.NumberOfElements() + rMPMModelPart.NumberOfElements();
        const unsigned int number_nodes = rBackgroundGridModelPart.NumberOfNodes();
        unsigned int last_condition_id;
        if (number_elements > number_nodes && number_elements > number_conditions)
            last_condition_id = number_elements + 1;
        else if (number_nodes > number_elements && number_nodes > number_conditions)
            last_condition_id = number_nodes + 1;
        else
            last_condition_id = number_conditions + 1;

        // Loop over the submodelpart of rBackgroundGridModelPart
        for (ModelPart::SubModelPartIterator submodelpart_it = rBackgroundGridModelPart.SubModelPartsBegin();
                submodelpart_it != rBackgroundGridModelPart.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart& submodelpart = *submodelpart_it;

            // For submodelpart without condition, exit
            if (submodelpart.NumberOfConditions() != 0){

                std::string submodelpart_name = submodelpart.Name();
                rMPMModelPart.CreateSubModelPart(submodelpart_name);

                // For regular conditions: straight copy all conditions
                if (!submodelpart.ConditionsBegin()->Is(BOUNDARY)){
                    rMPMModelPart.SetConditions(submodelpart.pConditions());
                    rMPMModelPart.GetSubModelPart(submodelpart_name).SetConditions(submodelpart.pConditions());
                }
                // For boundary conditions: create particle conditions for all the necessary conditions
                else{

                    // NOTE: To create Particle Condition, we consider both the nodal position as well as the position of integration point
                    // Loop over the conditions of submodelpart and generate mpm condition to be appended to the rMPMModelPart
                    for (ModelPart::ConditionIterator i = submodelpart.ConditionsBegin();
                            i != submodelpart.ConditionsEnd(); i++)
                    {
                        Properties::Pointer properties = i->pGetProperties();
                        const int mpc_condition_id = i->GetProperties().Id();

                        // Flag whether condition is Neumann or Dirichlet
                        const bool is_neumann_condition = i->GetValue(MPC_IS_NEUMANN);

                        // Check number of particles per condition to be created
                        unsigned int particles_per_condition = 0; // Default zero
                        if (i->Has( PARTICLES_PER_CONDITION )){
                            particles_per_condition = i->GetValue(PARTICLES_PER_CONDITION);
                        }
                        else{
                            std::string warning_msg = "PARTICLES_PER_CONDITION is not specified, ";
                            warning_msg += "Only using nodal position is assumed: 1 (Point), 2 (Line), 3 (Triangular), 4 (Quadrilateral)";
                            KRATOS_WARNING("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                        }

                        // Get condition variables:
                        // Normal vector
                        if (i->Has(NORMAL)) mpc_normal = i->GetValue(NORMAL);
                        ParticleMechanicsMathUtilities<double>::Normalize(mpc_normal);

                        // Get shape_function_values from defined particle_per_condition
                        auto& r_geometry = i->GetGeometry(); // current condition's geometry
                        const GeometryData::KratosGeometryType geo_type = r_geometry.GetGeometryType();
                        Matrix shape_functions_values;

                        // Get geometry and dimension of the background grid
                        std::string condition_type_name;
                        const GeometryData::KratosGeometryType background_geo_type = rBackgroundGridModelPart.ElementsBegin()->GetGeometry().GetGeometryType();
                        const std::size_t domain_size = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];

                        if (geo_type == GeometryData::Kratos_Point2D  || geo_type == GeometryData::Kratos_Point3D)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 1: // Only nodal
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Point" + std::to_string(domain_size) + "D.\n";
                                    warning_msg += "Available option is: 1 (default).\n";
                                    warning_msg += "The default number of particle: 1 is currently assumed.";
                                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition){
                                if (domain_size==2){
                                    if (background_geo_type == GeometryData::Kratos_Triangle2D3)
                                        condition_type_name = "MPMParticlePointLoadCondition2D3N";
                                    else if (background_geo_type == GeometryData::Kratos_Quadrilateral2D4)
                                        condition_type_name = "MPMParticlePointLoadCondition2D4N";
                                }
                                else if (domain_size==3){
                                    if (background_geo_type == GeometryData::Kratos_Tetrahedra3D4)
                                        condition_type_name = "MPMParticlePointLoadCondition3D4N";
                                    else if (background_geo_type == GeometryData::Kratos_Hexahedra3D8)
                                        condition_type_name = "MPMParticlePointLoadCondition3D8N";
                                }

                                if( i->Has( POINT_LOAD ) )
                                    point_load = i->GetValue( POINT_LOAD );
                            }

                        }
                        else if (geo_type == GeometryData::Kratos_Line2D2  || geo_type == GeometryData::Kratos_Line3D2)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 2: // Only nodal
                                    break;
                                case 3:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 4:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 5:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                    break;
                                case 6:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                case 7:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Line" + std::to_string(domain_size) + "D.\n";
                                    warning_msg += "Available options are: 2 (default), 3, 4, 5, 6, 7.\n";
                                    warning_msg += "The default number of particle: 2 is currently assumed.";
                                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle line load condition is not yet implemented." << std::endl;

                        }
                        else if(geo_type == GeometryData::Kratos_Triangle3D3)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 3: // Only nodal
                                    break;
                                case 4:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 6:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 9:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                case 15:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                    break;
                                case 19:
                                    shape_functions_values = MP16ShapeFunctions();
                                    break;
                                case 36:
                                    shape_functions_values = MP33ShapeFunctions();
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Triangular" + std::to_string(domain_size) + "D.\n";
                                    warning_msg += "Available options are: 3 (default), 4, 6, 9, 15, 19 and 36.\n";
                                    warning_msg += "The default number of particle: 3 is currently assumed.";
                                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                        }
                        else if(geo_type == GeometryData::Kratos_Quadrilateral3D4)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 4: // Only nodal
                                    break;
                                case 5:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 8:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 13:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                    break;
                                case 20:
                                    shape_functions_values = r_geometry.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Quadrilateral" + std::to_string(domain_size) + "D.\n";
                                    warning_msg += "Available options are: 4 (default), 5, 8, 13, and 20.\n";
                                    warning_msg += "The default number of particle: 4 is currently assumed.";
                                    KRATOS_INFO("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                        }
                        else{
                            std::string error_msg = "The Geometry type of the Condition given is invalid or currently not available. ";
                            error_msg += "Please remesh the problem domain to Point2D1N or Line2D2N for 2D or ";
                            error_msg += "Point3D1N, Line3D2N, Triangle3D3N or Quadrilateral3D4N for 3D.";
                            KRATOS_ERROR << error_msg << std::endl;
                        }

                        // Number of integration point per condition
                        const unsigned int integration_point_per_conditions = shape_functions_values.size1();

                        // Evaluation of geometric length/area
                        const double area = r_geometry.Area();
                        mpc_area = area / (1 + integration_point_per_conditions);
                        const double mpc_nodal_area = mpc_area / r_geometry.size();

                        // Check condition variables
                        if (i->Has(DISPLACEMENT))
                            mpc_imposed_displacement = i->GetValue(DISPLACEMENT);
                        if (i->Has(VELOCITY))
                            mpc_imposed_velocity = i->GetValue(VELOCITY);
                        if (i->Has(ACCELERATION))
                            mpc_imposed_acceleration = i->GetValue(ACCELERATION);
                        if (i->Has(PENALTY_FACTOR))
                            mpc_penalty_factor = i->GetValue(PENALTY_FACTOR);

                        const bool is_slip = i->Is(SLIP);
                        const bool is_contact = i->Is(CONTACT);
                        const bool is_interface = i->Is(INTERFACE);
                        const bool flip_normal_direction = i->Is(MODIFIED);

                        // If dirichlet boundary or coupling interface
                        if (!is_neumann_condition){
                            if(!is_interface){
                                if (domain_size==2){
                                    if (background_geo_type == GeometryData::Kratos_Triangle2D3)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
                                    else if (background_geo_type == GeometryData::Kratos_Quadrilateral2D4)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition2D4N";
                                }
                                else if (domain_size==3){
                                    if (background_geo_type == GeometryData::Kratos_Tetrahedra3D4)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition3D4N";
                                    else if (background_geo_type == GeometryData::Kratos_Hexahedra3D8)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition3D8N";
                                }
                            }
                            else{
                                if (domain_size==2){
                                    if (background_geo_type == GeometryData::Kratos_Triangle2D3)
                                        condition_type_name = "MPMParticlePenaltyCouplingInterfaceCondition2D3N";
                                    else if (background_geo_type == GeometryData::Kratos_Quadrilateral2D4)
                                        condition_type_name = "MPMParticlePenaltyCouplingInterfaceCondition2D4N";
                                }
                                else if (domain_size==3){
                                    if (background_geo_type == GeometryData::Kratos_Tetrahedra3D4)
                                        condition_type_name = "MPMParticlePenaltyCouplingInterfaceCondition3D4N";
                                    else if (background_geo_type == GeometryData::Kratos_Hexahedra3D8)
                                        condition_type_name = "MPMParticlePenaltyCouplingInterfaceCondition3D8N";
                                }
                            }
                        }

                        // Get new condition
                        const Condition& new_condition = KratosComponents<Condition>::Get(condition_type_name);

                        // 1. Loop over the conditions to create inner particle condition
                        unsigned int new_condition_id = 0;
                        for ( unsigned int point_number = 0; point_number < integration_point_per_conditions; point_number++ )
                        {
                            // Create new material point condition
                            new_condition_id = last_condition_id + point_number;
                            Condition::Pointer p_condition = new_condition.Create(new_condition_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                            mpc_xg.clear();

                            // Loop over the nodes of the grid condition
                            for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++){
                                for ( unsigned int j = 0; j < r_geometry.size(); j ++){
                                    mpc_xg[dimension] = mpc_xg[dimension] + shape_functions_values(point_number, j) * r_geometry[j].Coordinates()[dimension];
                                }
                            }

                            // Check Normal direction
                            if (flip_normal_direction) mpc_normal *= -1.0;

                            // Setting particle condition's initial condition
                            // TODO: If any variable is added or remove here, please add and remove also at the second loop below
                            p_condition->SetValue(MPC_CONDITION_ID, mpc_condition_id);
                            p_condition->SetValue(MPC_COORD, mpc_xg);
                            p_condition->SetValue(MPC_AREA, mpc_area);
                            p_condition->SetValue(MPC_NORMAL, mpc_normal);
                            p_condition->SetValue(MPC_DISPLACEMENT, mpc_displacement);
                            p_condition->SetValue(MPC_IMPOSED_DISPLACEMENT, mpc_imposed_displacement);
                            p_condition->SetValue(MPC_VELOCITY, mpc_velocity);
                            p_condition->SetValue(MPC_IMPOSED_VELOCITY, mpc_imposed_velocity);
                            p_condition->SetValue(MPC_ACCELERATION, mpc_acceleration);
                            p_condition->SetValue(MPC_IMPOSED_ACCELERATION, mpc_imposed_acceleration);

                            if (is_neumann_condition)
                                p_condition->SetValue(POINT_LOAD, point_load);
                            else{
                                p_condition->SetValue(PENALTY_FACTOR, mpc_penalty_factor);
                                if (is_slip)
                                    p_condition->Set(SLIP);
                                if (is_contact)
                                    p_condition->Set(CONTACT);
                                if (is_interface)
                                {
                                    p_condition->Set(INTERFACE);
                                    p_condition->SetValue(MPC_CONTACT_FORCE, mpc_contact_force);
                                }
                            }

                            // Add the MP Condition to the model part
                            rMPMModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                        }

                        last_condition_id += integration_point_per_conditions;

                        // 2. Loop over the nodes associated to each condition to create nodal particle condition
                        for ( unsigned int j = 0; j < r_geometry.size(); j ++)
                        {
                            // Nodal normal vector is used
                            if (r_geometry[j].Has(NORMAL)) mpc_normal = r_geometry[j].FastGetSolutionStepValue(NORMAL);
                            const double denominator = std::sqrt(mpc_normal[0]*mpc_normal[0] + mpc_normal[1]*mpc_normal[1] + mpc_normal[2]*mpc_normal[2]);
                            if (std::abs(denominator) > std::numeric_limits<double>::epsilon() ) mpc_normal *= 1.0 / denominator;

                            // Check Normal direction
                            if (flip_normal_direction) mpc_normal *= -1.0;

                            // Create new material point condition
                            new_condition_id = last_condition_id + j;
                            Condition::Pointer p_condition = new_condition.Create(new_condition_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                            mpc_xg.clear();
                            for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++){
                                mpc_xg[dimension] = r_geometry[j].Coordinates()[dimension];
                            }

                            // Setting particle condition's initial condition
                            // TODO: If any variable is added or remove here, please add and remove also at the first loop above
                            p_condition->SetValue(MPC_CONDITION_ID, mpc_condition_id);
                            p_condition->SetValue(MPC_COORD, mpc_xg);
                            p_condition->SetValue(MPC_AREA, mpc_nodal_area);
                            p_condition->SetValue(MPC_NORMAL, mpc_normal);
                            p_condition->SetValue(MPC_DISPLACEMENT, mpc_displacement);
                            p_condition->SetValue(MPC_IMPOSED_DISPLACEMENT, mpc_imposed_displacement);
                            p_condition->SetValue(MPC_VELOCITY, mpc_velocity);
                            p_condition->SetValue(MPC_IMPOSED_VELOCITY, mpc_imposed_velocity);
                            p_condition->SetValue(MPC_ACCELERATION, mpc_acceleration);
                            p_condition->SetValue(MPC_IMPOSED_ACCELERATION, mpc_imposed_acceleration);

                            if (is_neumann_condition)
                                p_condition->SetValue(POINT_LOAD, point_load);
                            else{
                                p_condition->SetValue(PENALTY_FACTOR, mpc_penalty_factor);
                                if (is_slip)
                                    p_condition->Set(SLIP);
                                if (is_contact)
                                    p_condition->Set(CONTACT);
                                if (is_interface)
                                {
                                    p_condition->Set(INTERFACE);
                                    p_condition->SetValue(MPC_CONTACT_FORCE, mpc_contact_force);
                                }
                            }

                            // Add the MP Condition to the model part
                            rMPMModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                        }

                        last_condition_id += r_geometry.size();
                    }
                }
            }
        }
    }

/***********************************************************************************/
/***********************************************************************************/

    Matrix MP16ShapeFunctions()
    {
        const double Na1 = 0.33333333333333;
        const double Nb1 = 0.45929258829272;
        const double Nb2 = 0.08141482341455;
        const double Nc1 = 0.17056930775176;
        const double Nc2 = 0.65886138449648;

        const double Nd1 = 0.05054722831703;
        const double Nd2 = 0.89890554336594;

        const double Ne1 = 0.26311282963464;
        const double Ne2 = 0.72849239295540;
        const double Ne3 = 0.00839477740996;

        BoundedMatrix<double,16,3> MP_shape_functions;

        MP_shape_functions(0,0) = Na1;
        MP_shape_functions(0,1) = Na1;
        MP_shape_functions(0,2) = Na1;

        MP_shape_functions(1,0) = Nb1;
        MP_shape_functions(1,1) = Nb1;
        MP_shape_functions(1,2) = Nb2;

        MP_shape_functions(2,0) = Nb1;
        MP_shape_functions(2,1) = Nb2;
        MP_shape_functions(2,2) = Nb1;

        MP_shape_functions(3,0) = Nb2;
        MP_shape_functions(3,1) = Nb1;
        MP_shape_functions(3,2) = Nb1;

        MP_shape_functions(4,0) = Nc1;
        MP_shape_functions(4,1) = Nc1;
        MP_shape_functions(4,2) = Nc2;

        MP_shape_functions(5,0) = Nc1;
        MP_shape_functions(5,1) = Nc2;
        MP_shape_functions(5,2) = Nc1;

        MP_shape_functions(6,0) = Nc2;
        MP_shape_functions(6,1) = Nc1;
        MP_shape_functions(6,2) = Nc1;

        MP_shape_functions(7,0) = Nd1;
        MP_shape_functions(7,1) = Nd1;
        MP_shape_functions(7,2) = Nd2;

        MP_shape_functions(8,0) = Nd1;
        MP_shape_functions(8,1) = Nd2;
        MP_shape_functions(8,2) = Nd1;

        MP_shape_functions(9,0) = Nd2;
        MP_shape_functions(9,1) = Nd1;
        MP_shape_functions(9,2) = Nd1;

        MP_shape_functions(10,0) = Ne1;
        MP_shape_functions(10,1) = Ne2;
        MP_shape_functions(10,2) = Ne3;

        MP_shape_functions(11,0) = Ne2;
        MP_shape_functions(11,1) = Ne3;
        MP_shape_functions(11,2) = Ne1;

        MP_shape_functions(12,0) = Ne3;
        MP_shape_functions(12,1) = Ne1;
        MP_shape_functions(12,2) = Ne2;

        MP_shape_functions(13,0) = Ne2;
        MP_shape_functions(13,1) = Ne1;
        MP_shape_functions(13,2) = Ne3;

        MP_shape_functions(14,0) = Ne1;
        MP_shape_functions(14,1) = Ne3;
        MP_shape_functions(14,2) = Ne2;

        MP_shape_functions(15,0) = Ne3;
        MP_shape_functions(15,1) = Ne2;
        MP_shape_functions(15,2) = Ne1;

        //MP_shape_functions = [(Na1, Na1, Na1),(Nb1, Nb1, Nb2),(Nb1, Nb2, Nb1),(Nb2, Nb1, Nb1),
        //                    (Nc1, Nc1, Nc2),(Nc1, Nc2, Nc1),(Nc2, Nc1, Nc1),(Nd1, Nd1, Nd2),
        //                    (Nd1, Nd2, Nd1),(Nd2, Nd1, Nd1),(Ne1, Ne2, Ne3),(Ne2, Ne3, Ne1),
        //                    (Ne3, Ne1, Ne2),(Ne2, Ne1, Ne3),(Ne1, Ne3, Ne2),(Ne3, Ne2, Ne1)];

        return MP_shape_functions;
    }

/***********************************************************************************/
/***********************************************************************************/

    Matrix MP33ShapeFunctions()
    {
        const double Na2 = 0.02356522045239;
        const double Na1 = 0.488217389773805;

        const double Nb2 = 0.120551215411079;
        const double Nb1 = 0.43972439229446;

        const double Nc2 = 0.457579229975768;
        const double Nc1 = 0.271210385012116;

        const double Nd2 = 0.744847708916828;
        const double Nd1 = 0.127576145541586;

        const double Ne2 = 0.957365299093579;
        const double Ne1 = 0.021317350453210;

        const double Nf1 = 0.115343494534698;
        const double Nf2 = 0.275713269685514;
        const double Nf3 = 0.608943235779788;

        const double Ng1 = 0.022838332222257;
        const double Ng2 = 0.281325580989940;
        const double Ng3 = 0.695836086787803;

        const double Nh1 = 0.025734050548330;
        const double Nh2 = 0.116251915907597;
        const double Nh3 = 0.858014033544073;

        BoundedMatrix<double,33,3> MP_shape_functions;

        MP_shape_functions(0,0) = Na1;
        MP_shape_functions(0,1) = Na1;
        MP_shape_functions(0,2) = Na2;

        MP_shape_functions(1,0) = Na1;
        MP_shape_functions(1,1) = Na2;
        MP_shape_functions(1,2) = Na1;

        MP_shape_functions(2,0) = Na2;
        MP_shape_functions(2,1) = Na1;
        MP_shape_functions(2,2) = Na1;


        MP_shape_functions(3,0) = Nb1;
        MP_shape_functions(3,1) = Nb1;
        MP_shape_functions(3,2) = Nb2;

        MP_shape_functions(4,0) = Nb1;
        MP_shape_functions(4,1) = Nb2;
        MP_shape_functions(4,2) = Nb1;

        MP_shape_functions(5,0) = Nb2;
        MP_shape_functions(5,1) = Nb1;
        MP_shape_functions(5,2) = Nb1;

        MP_shape_functions(6,0) = Nc1;
        MP_shape_functions(6,1) = Nc1;
        MP_shape_functions(6,2) = Nc2;

        MP_shape_functions(7,0) = Nc1;
        MP_shape_functions(7,1) = Nc2;
        MP_shape_functions(7,2) = Nc1;

        MP_shape_functions(8,0) = Nc2;
        MP_shape_functions(8,1) = Nc1;
        MP_shape_functions(8,2) = Nc1;

        MP_shape_functions(9,0) = Nd1;
        MP_shape_functions(9,1) = Nd1;
        MP_shape_functions(9,2) = Nd2;

        MP_shape_functions(10,0) = Nd1;
        MP_shape_functions(10,1) = Nd2;
        MP_shape_functions(10,2) = Nd1;

        MP_shape_functions(11,0) = Nd2;
        MP_shape_functions(11,1) = Nd1;
        MP_shape_functions(11,2) = Nd1;

        MP_shape_functions(12,0) = Ne1;
        MP_shape_functions(12,1) = Ne1;
        MP_shape_functions(12,2) = Ne2;

        MP_shape_functions(13,0) = Ne1;
        MP_shape_functions(13,1) = Ne2;
        MP_shape_functions(13,2) = Ne1;

        MP_shape_functions(14,0) = Ne2;
        MP_shape_functions(14,1) = Ne1;
        MP_shape_functions(14,2) = Ne1;

        MP_shape_functions(15,0) = Nf1;
        MP_shape_functions(15,1) = Nf2;
        MP_shape_functions(15,2) = Nf3;

        MP_shape_functions(16,0) = Nf2;
        MP_shape_functions(16,1) = Nf3;
        MP_shape_functions(16,2) = Nf1;

        MP_shape_functions(17,0) = Nf3;
        MP_shape_functions(17,1) = Nf1;
        MP_shape_functions(17,2) = Nf2;

        MP_shape_functions(18,0) = Nf2;
        MP_shape_functions(18,1) = Nf1;
        MP_shape_functions(18,2) = Nf3;

        MP_shape_functions(19,0) = Nf1;
        MP_shape_functions(19,1) = Nf3;
        MP_shape_functions(19,2) = Nf2;

        MP_shape_functions(20,0) = Nf3;
        MP_shape_functions(20,1) = Nf2;
        MP_shape_functions(20,2) = Nf1;

        MP_shape_functions(21,0) = Ng1;
        MP_shape_functions(21,1) = Ng2;
        MP_shape_functions(21,2) = Ng3;

        MP_shape_functions(22,0) = Ng2;
        MP_shape_functions(22,1) = Ng3;
        MP_shape_functions(22,2) = Ng1;

        MP_shape_functions(23,0) = Ng3;
        MP_shape_functions(23,1) = Ng1;
        MP_shape_functions(23,2) = Ng2;

        MP_shape_functions(24,0) = Ng2;
        MP_shape_functions(24,1) = Ng1;
        MP_shape_functions(24,2) = Ng3;

        MP_shape_functions(25,0) = Ng1;
        MP_shape_functions(25,1) = Ng3;
        MP_shape_functions(25,2) = Ng2;

        MP_shape_functions(26,0) = Ng3;
        MP_shape_functions(26,1) = Ng2;
        MP_shape_functions(26,2) = Ng1;

        MP_shape_functions(27,0) = Nh1;
        MP_shape_functions(27,1) = Nh2;
        MP_shape_functions(27,2) = Nh3;

        MP_shape_functions(28,0) = Nh2;
        MP_shape_functions(28,1) = Nh3;
        MP_shape_functions(28,2) = Nh1;

        MP_shape_functions(29,0) = Nh3;
        MP_shape_functions(29,1) = Nh1;
        MP_shape_functions(29,2) = Nh2;

        MP_shape_functions(30,0) = Nh2;
        MP_shape_functions(30,1) = Nh1;
        MP_shape_functions(30,2) = Nh3;

        MP_shape_functions(31,0) = Nh1;
        MP_shape_functions(31,1) = Nh3;
        MP_shape_functions(31,2) = Nh2;

        MP_shape_functions(32,0) = Nh3;
        MP_shape_functions(32,1) = Nh2;
        MP_shape_functions(32,2) = Nh1;

        return MP_shape_functions;
    }

} // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos




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


#ifndef KRATOS_MPM_CREATE_ELEMENT_CONDITION_UTILITY
#define KRATOS_MPM_CREATE_ELEMENT_CONDITION_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
    namespace MpmCreateElementConditionUtility
    {

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        /**
        * @brief Function to Initiate material point element.
        * @details It is designed to be called ONCE by the class constructor.
        */
        static void CreateMaterialPointElement(
            ModelPart& rBackgroundGridModelPart,
            ModelPart& rMpmModelPart,
            ModelPart& rMpmInitialModelPart,
            Element const& rNewElement,
            bool IsMixedFormulation = false)
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
            const unsigned int number_elements = rBackgroundGridModelPart.NumberOfElements() + rMpmInitialModelPart.NumberOfElements();
            const unsigned int number_nodes = rBackgroundGridModelPart.NumberOfNodes();
            unsigned int last_element_id = (number_nodes > number_elements) ? (number_nodes + 1) : (number_elements + 1);

            // Loop over the submodelpart of rMpmInitialModelPart
            for (ModelPart::SubModelPartIterator submodelpart_it = rMpmInitialModelPart.SubModelPartsBegin();
                submodelpart_it != rMpmInitialModelPart.SubModelPartsEnd(); submodelpart_it++)
            {
                ModelPart&  submodelpart = *submodelpart_it;
                std::string submodelpart_name = submodelpart.Name();

                KRATOS_WATCH(submodelpart_name)

                rMpmModelPart.CreateSubModelPart(submodelpart_name);

                KRATOS_WATCH(submodelpart.NumberOfElements())
                // Loop over the element of submodelpart and generate mpm element to be appended to the rMpmModelPart
                for (ModelPart::ElementIterator i = submodelpart.ElementsBegin();
                    i != submodelpart.ElementsEnd(); i++)
                {
                    if (true)//(i->IsDefined(ACTIVE))
                    {
                        Properties::Pointer properties = i->pGetProperties();
                        const int material_id = i->GetProperties().Id();
                        const double density = i->GetProperties()[DENSITY];

                        // Check number of particles per element to be created
                        unsigned int particles_per_element;
                        if (i->GetProperties().Has(PARTICLES_PER_ELEMENT)) {
                            particles_per_element = i->GetProperties()[PARTICLES_PER_ELEMENT];
                        }
                        else {
                            std::string warning_msg = "PARTICLES_PER_ELEMENT is not specified in Properties, ";
                            warning_msg += "1 Particle per element is assumed.";
                            KRATOS_WARNING("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                            particles_per_element = 1;
                        }

                        KRATOS_WATCH(particles_per_element)

                        const Geometry< Node < 3 > >& rGeom = i->GetGeometry(); // current element's geometry
                        const GeometryData::KratosGeometryType rGeoType = rGeom.GetGeometryType();
                        Matrix shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                        if (rGeoType == GeometryData::Kratos_Tetrahedra3D4 || rGeoType == GeometryData::Kratos_Triangle2D3)
                        {
                            switch (particles_per_element)
                            {
                            case 1:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                                break;
                            case 3:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                break;
                            case 6:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                                break;
                            case 12:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_5);
                                break;
                            //case 16:
                            //    if (TDim == 2) {
                            //        shape_functions_values = this->MP16ShapeFunctions();
                            //        break;
                            //    }
                            //case 33:
                            //    if (TDim == 2) {
                            //        shape_functions_values = this->MP33ShapeFunctions();
                            //        break;
                            //    }
                            default:
                                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                                warning_msg += " is not available for Triangular 3D.\n";
                                warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                                warning_msg += "The default number of particle: 3 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                            }
                        }
                        else if (rGeoType == GeometryData::Kratos_Hexahedra3D8 || rGeoType == GeometryData::Kratos_Quadrilateral2D4)
                        {
                            switch (particles_per_element)
                            {
                            case 1:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                                break;
                            case 4:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                break;
                            case 9:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_3);
                                break;
                            case 16:
                                shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                                break;
                            default:
                                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                                warning_msg += " is not available for Quadrilateral3D.\n";
                                warning_msg += "Available options are: 1, 4, 9, 16.\n";
                                warning_msg += "The default number of particle: 4 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                            }
                        }
                        else {
                            std::string error_msg = "The Geometry type of the Element given is invalid or currently not available. ";
                            error_msg += "Please remesh the problem domain to Triangle2D3N or Quadrilateral2D4N for 2D or ";
                            error_msg += "Tetrahedral3D4N or Hexahedral3D8N for 3D.";
                            KRATOS_ERROR << error_msg << std::endl;
                        }

                        // Number of MP per elements
                        const unsigned int integration_point_per_elements = shape_functions_values.size1();

                        // Evaluation of element area/volume
                        const double area = rGeom.Area();
                        if (rGeom.WorkingSpaceDimension() == 2 && i->GetProperties().Has(THICKNESS)) {
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
                            // Create new material point element
                            new_element_id = last_element_id + PointNumber;
                            Element::Pointer p_element = rNewElement.Create(new_element_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                            const double MP_Density = density;
                            const int MP_Material_Id = material_id;

                            xg.clear();

                            // Loop over the nodes of the grid element
                            for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++)
                            {
                                for (unsigned int j = 0; j < rGeom.size(); j++)
                                {
                                    xg[dim] = xg[dim] + shape_functions_values(PointNumber, j) * rGeom[j].Coordinates()[dim];
                                }
                            }

                            // Setting particle element's initial condition
                            p_element->SetValue(MP_MATERIAL_ID, MP_Material_Id);
                            p_element->SetValue(MP_DENSITY, MP_Density);
                            p_element->SetValue(MP_MASS, mp_mass);
                            p_element->SetValue(MP_VOLUME, mp_volume);
                            p_element->SetValue(MP_COORD, xg);
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
                            rMpmModelPart.GetSubModelPart(submodelpart_name).AddElement(p_element);
                        }

                        last_element_id += integration_point_per_elements;

                    }
                }
            }
        }

        /**
        * @brief Function to Initiate material point condition.
        * @details It is designed to be called ONCE by the class constructor.
        */
        static void CreateMaterialPointCondition(
            ModelPart& rBackgroundGridModelPart,
            ModelPart& rMpmModelPart,
            ModelPart& rMpmInitialModelPart)
        {
            // Initialize zero the variables needed
            array_1d<double, 3> mpc_xg = ZeroVector(3);
            array_1d<double, 3> mpc_normal = ZeroVector(3);
            array_1d<double, 3> mpc_displacement = ZeroVector(3);
            array_1d<double, 3> mpc_velocity = ZeroVector(3);
            array_1d<double, 3> mpc_acceleration = ZeroVector(3);
            array_1d<double, 3 > point_load = ZeroVector(3);

            double TDim = 2;

            double mpc_area = 0.0;
            double mpc_penalty_factor = 0.0;

            // Determine condition index: This convention is done in order for the purpose of visualization in GiD
            const unsigned int number_conditions = rBackgroundGridModelPart.NumberOfConditions();
            const unsigned int number_elements = rBackgroundGridModelPart.NumberOfElements() + rMpmInitialModelPart.NumberOfElements() + rMpmModelPart.NumberOfElements();
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
                if (submodelpart.NumberOfConditions() != 0) {

                    std::string submodelpart_name = submodelpart.Name();
                    rMpmModelPart.CreateSubModelPart(submodelpart_name);

                    // For regular conditions: straight copy all conditions
                    if (!submodelpart.ConditionsBegin()->Is(BOUNDARY)) {
                        rMpmModelPart.SetConditions(submodelpart.pConditions());
                        rMpmModelPart.GetSubModelPart(submodelpart_name).SetConditions(submodelpart.pConditions());
                    }
                    // For boundary conditions: create particle conditions for all the necessary conditions
                    else {

                        // NOTE: To create Particle Condition, we consider both the nodal position as well as the position of integration point
                        // Loop over the conditions of submodelpart and generate mpm condition to be appended to the rMpmModelPart
                        for (ModelPart::ConditionIterator i = submodelpart.ConditionsBegin();
                            i != submodelpart.ConditionsEnd(); i++)
                        {
                            Properties::Pointer properties = i->pGetProperties();
                            const int mpc_condition_id = i->GetProperties().Id();

                            // Flag whether condition is Neumann or Dirichlet
                            const bool is_neumann_condition = i->GetValue(MPC_IS_NEUMANN);

                            // Check number of particles per condition to be created
                            unsigned int particles_per_condition = 0; // Default zero
                            if (i->Has(PARTICLES_PER_CONDITION)) {
                                particles_per_condition = i->GetValue(PARTICLES_PER_CONDITION);
                            }
                            else {
                                std::string warning_msg = "PARTICLES_PER_CONDITION is not specified, ";
                                warning_msg += "Only using nodal position is assumed: 1 (Point), 2 (Line), 3 (Triangular), 4 (Quadrilateral)";
                                KRATOS_WARNING("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                            }

                            // Get condition variables:
                            // Normal vector (normalized)
                            if (i->Has(NORMAL)) mpc_normal = i->GetValue(NORMAL);
                            const double denominator = std::sqrt(mpc_normal[0] * mpc_normal[0] + mpc_normal[1] * mpc_normal[1] + mpc_normal[2] * mpc_normal[2]);
                            if (std::abs(denominator) > std::numeric_limits<double>::epsilon()) mpc_normal *= 1.0 / denominator;

                            // Get shape_function_values from defined particle_per_condition
                            auto& rGeom = i->GetGeometry(); // current condition's geometry
                            const GeometryData::KratosGeometryType rGeoType = rGeom.GetGeometryType();
                            Matrix shape_functions_values;

                            //Get geometry of the background grid
                            std::string condition_type_name;
                            const GeometryData::KratosGeometryType rBackgroundGeoType = rBackgroundGridModelPart.ElementsBegin()->GetGeometry().GetGeometryType();

                            if (rGeoType == GeometryData::Kratos_Point2D || rGeoType == GeometryData::Kratos_Point3D)
                            {
                                switch (particles_per_condition)
                                {
                                case 0: // Default case
                                    break;
                                case 1: // Only nodal
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Point3D.\n";
                                    warning_msg += "Available option is: 1 (default).\n";
                                    warning_msg += "The default number of particle: 1 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                                }

                                if (is_neumann_condition) {

                                    if (TDim == 2) {
                                        if (rBackgroundGeoType == GeometryData::Kratos_Triangle2D3)
                                            condition_type_name = "MPMParticlePointLoadCondition2D3N";
                                        else if (rBackgroundGeoType == GeometryData::Kratos_Quadrilateral2D4)
                                            condition_type_name = "MPMParticlePointLoadCondition2D4N";
                                    }
                                    else if (TDim == 3) {
                                        if (rBackgroundGeoType == GeometryData::Kratos_Tetrahedra3D4)
                                            condition_type_name = "MPMParticlePointLoadCondition3D4N";
                                        else if (rBackgroundGeoType == GeometryData::Kratos_Hexahedra3D8)
                                            condition_type_name = "MPMParticlePointLoadCondition3D8N";
                                    }

                                    if (i->Has(POINT_LOAD))
                                        point_load = i->GetValue(POINT_LOAD);
                                }

                            }
                            else if (rGeoType == GeometryData::Kratos_Line2D2 || rGeoType == GeometryData::Kratos_Line3D2)
                            {
                                switch (particles_per_condition)
                                {
                                case 0: // Default case
                                    break;
                                case 2: // Only nodal
                                    break;
                                case 3:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                                    break;
                                case 4:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                    break;
                                case 5:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_3);
                                    break;
                                case 6:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                                    break;
                                case 7:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_5);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Line" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 2 (default), 3, 4, 5, 6, 7.\n";
                                    warning_msg += "The default number of particle: 2 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                                }

                                if (is_neumann_condition)
                                    KRATOS_ERROR << "Particle line load condition is not yet implemented." << std::endl;

                            }
                            else if (rGeoType == GeometryData::Kratos_Triangle3D3)
                            {
                                switch (particles_per_condition)
                                {
                                case 0: // Default case
                                    break;
                                case 3: // Only nodal
                                    break;
                                case 4:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                                    break;
                                case 6:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                    break;
                                case 9:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                                    break;
                                case 15:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_5);
                                    break;
                                //case 19:
                                //    shape_functions_values = this->MP16ShapeFunctions();
                                //    break;
                                //case 36:
                                //    shape_functions_values = this->MP33ShapeFunctions();
                                //    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Triangular" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 3 (default), 4, 6, 9, 15, 19 and 36.\n";
                                    warning_msg += "The default number of particle: 3 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                                }

                                if (is_neumann_condition)
                                    KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                            }
                            else if (rGeoType == GeometryData::Kratos_Quadrilateral3D4)
                            {
                                switch (particles_per_condition)
                                {
                                case 0: // Default case
                                    break;
                                case 4: // Only nodal
                                    break;
                                case 5:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                                    break;
                                case 8:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                                    break;
                                case 13:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_3);
                                    break;
                                case 20:
                                    shape_functions_values = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Quadrilateral" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 4 (default), 5, 8, 13, and 20.\n";
                                    warning_msg += "The default number of particle: 4 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                                }

                                if (is_neumann_condition)
                                    KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                            }
                            else {
                                std::string error_msg = "The Geometry type of the Condition given is invalid or currently not available. ";
                                error_msg += "Please remesh the problem domain to Point2D1N or Line2D2N for 2D or ";
                                error_msg += "Point3D1N, Line3D2N, Triangle3D3N or Quadrilateral3D4N for 3D.";
                                KRATOS_ERROR << error_msg << std::endl;
                            }

                            // Number of integration point per condition
                            const unsigned int integration_point_per_conditions = shape_functions_values.size1();

                            // Evaluation of geometric length/area
                            const double area = rGeom.Area();
                            mpc_area = area / (rGeom.size() + integration_point_per_conditions);

                            // Check condition variables
                            if (i->Has(DISPLACEMENT))
                                mpc_displacement = i->GetValue(DISPLACEMENT);
                            if (i->Has(VELOCITY))
                                mpc_velocity = i->GetValue(VELOCITY);
                            if (i->Has(ACCELERATION))
                                mpc_acceleration = i->GetValue(ACCELERATION);
                            if (i->Has(PENALTY_FACTOR))
                                mpc_penalty_factor = i->GetValue(PENALTY_FACTOR);
                            const bool is_slip = i->Is(SLIP);
                            const bool is_contact = i->Is(CONTACT);
                            const bool flip_normal_direction = i->Is(MODIFIED);

                            // If dirichlet boundary
                            if (!is_neumann_condition) {
                                if (TDim == 2) {
                                    if (rBackgroundGeoType == GeometryData::Kratos_Triangle2D3)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
                                    else if (rBackgroundGeoType == GeometryData::Kratos_Quadrilateral2D4)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition2D4N";
                                }
                                else if (TDim == 3) {
                                    if (rBackgroundGeoType == GeometryData::Kratos_Tetrahedra3D4)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition3D4N";
                                    else if (rBackgroundGeoType == GeometryData::Kratos_Hexahedra3D8)
                                        condition_type_name = "MPMParticlePenaltyDirichletCondition3D8N";
                                }
                            }

                            // Get new condition
                            const Condition& new_condition = KratosComponents<Condition>::Get(condition_type_name);

                            // 1. Loop over the conditions to create inner particle condition
                            unsigned int new_condition_id = 0;
                            for (unsigned int point_number = 0; point_number < integration_point_per_conditions; point_number++)
                            {
                                // Create new material point condition
                                new_condition_id = last_condition_id + point_number;
                                Condition::Pointer p_condition = new_condition.Create(new_condition_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                                mpc_xg.clear();

                                // Loop over the nodes of the grid condition
                                for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++) {
                                    for (unsigned int j = 0; j < rGeom.size(); j++) {
                                        mpc_xg[dim] = mpc_xg[dim] + shape_functions_values(point_number, j) * rGeom[j].Coordinates()[dim];
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
                                p_condition->SetValue(MPC_VELOCITY, mpc_velocity);
                                p_condition->SetValue(MPC_ACCELERATION, mpc_acceleration);

                                if (is_neumann_condition)
                                    p_condition->SetValue(POINT_LOAD, point_load);
                                else {
                                    p_condition->SetValue(PENALTY_FACTOR, mpc_penalty_factor);
                                    if (is_slip)
                                        p_condition->Set(SLIP);
                                    if (is_contact)
                                        p_condition->Set(CONTACT);
                                }

                                // Add the MP Condition to the model part
                                rMpmModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                            }

                            last_condition_id += integration_point_per_conditions;

                            // 2. Loop over the nodes associated to each condition to create nodal particle condition
                            for (unsigned int j = 0; j < rGeom.size(); j++)
                            {
                                // Nodal normal vector is used
                                if (rGeom[j].Has(NORMAL)) mpc_normal = rGeom[j].FastGetSolutionStepValue(NORMAL);
                                const double denominator = std::sqrt(mpc_normal[0] * mpc_normal[0] + mpc_normal[1] * mpc_normal[1] + mpc_normal[2] * mpc_normal[2]);
                                if (std::abs(denominator) > std::numeric_limits<double>::epsilon()) mpc_normal *= 1.0 / denominator;

                                // Check Normal direction
                                if (flip_normal_direction) mpc_normal *= -1.0;

                                // Create new material point condition
                                new_condition_id = last_condition_id + j;
                                Condition::Pointer p_condition = new_condition.Create(new_condition_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);

                                mpc_xg.clear();
                                for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++) {
                                    mpc_xg[dim] = rGeom[j].Coordinates()[dim];
                                }

                                // Setting particle condition's initial condition
                                // TODO: If any variable is added or remove here, please add and remove also at the first loop above
                                p_condition->SetValue(MPC_CONDITION_ID, mpc_condition_id);
                                p_condition->SetValue(MPC_COORD, mpc_xg);
                                p_condition->SetValue(MPC_AREA, mpc_area);
                                p_condition->SetValue(MPC_NORMAL, mpc_normal);
                                p_condition->SetValue(MPC_DISPLACEMENT, mpc_displacement);
                                p_condition->SetValue(MPC_VELOCITY, mpc_velocity);
                                p_condition->SetValue(MPC_ACCELERATION, mpc_acceleration);

                                if (is_neumann_condition)
                                    p_condition->SetValue(POINT_LOAD, point_load);
                                else {
                                    p_condition->SetValue(PENALTY_FACTOR, mpc_penalty_factor);
                                    if (is_slip)
                                        p_condition->Set(SLIP);
                                    if (is_contact)
                                        p_condition->Set(CONTACT);
                                }

                                // Add the MP Condition to the model part
                                rMpmModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                            }

                            last_condition_id += rGeom.size();
                        }
                    }
                }
            }
        }

    }
} // end namespace Kratos

#endif // KRATOS_MPM_CREATE_ELEMENT_CONDITION_UTILITY



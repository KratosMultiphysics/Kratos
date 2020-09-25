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
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"
#include "particle_mechanics_application_variables.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"


namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Geometry< Node<3> > GeometryType;

    typedef GeometryData::IntegrationMethod IntegrationMethod;

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

    /// Get integration weights of the geometry for the given integration method
    void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) GetIntegrationPointVolumes(const GeometryType& rGeom, const IntegrationMethod IntegrationMethod, Vector& rIntVolumes);

    /// Get integration method and shape function values for the given element
    void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) DetermineIntegrationMethodAndShapeFunctionValues(const GeometryType& rGeom, const SizeType ParticlesPerElement,
        IntegrationMethod& rIntegrationMethod, Matrix& rN, bool& IsEqualVolumes);

    /// Get integration weights/area of the geometry for the given integration method
    void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) GetIntegrationPointArea(const GeometryType& rGeom, const IntegrationMethod IntegrationMethod, Vector& rIntVolumes);

    /// Get integration method and shape function values for the given condition
    void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) DetermineConditionIntegrationMethodAndShapeFunctionValues(const GeometryType& rGeom, const SizeType ParticlesPerCondition,
        IntegrationMethod& rIntegrationMethod, Matrix& rN, bool& IsEqualVolumes);

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

                    // Get integration method and shape function values
                    IntegrationMethod int_method = GeometryData::GI_GAUSS_1;
                    Matrix shape_functions_values;
                    bool is_equal_int_volumes = false;
                    DetermineIntegrationMethodAndShapeFunctionValues(r_geometry, particles_per_element,
                        int_method, shape_functions_values, is_equal_int_volumes);

                    // Get volumes of the material points
                    const unsigned int integration_point_per_elements = shape_functions_values.size1();
                    Vector int_volumes (integration_point_per_elements);
                    if (is_equal_int_volumes) {
                        for (size_t j = 0; j < integration_point_per_elements; ++j)  int_volumes[j] = r_geometry.DomainSize() / integration_point_per_elements;
                    }
                    else  GetIntegrationPointVolumes(r_geometry, int_method, int_volumes);
                    if (domain_size == 2 && i->GetProperties().Has(THICKNESS)) {
                        for (size_t j = 0; j < integration_point_per_elements; ++j) int_volumes[j] *= i->GetProperties()[THICKNESS];
                    }

                    // Set element type
                    std::string element_type_name = "UpdatedLagrangian";
                    if (IsMixedFormulation) {
                        if (background_geo_type == GeometryData::Kratos_Triangle2D3) element_type_name = "UpdatedLagrangianUP";
                        else KRATOS_ERROR << "Element for mixed U-P formulation is only implemented for 2D Triangle Elements." << std::endl;
                    }
                    else if (IsAxisSymmetry && domain_size == 3) KRATOS_ERROR << "Axisymmetric elements must be used in a 2D domain. You specified a 3D domain." << std::endl;
                    else if (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM)) {
                        if (rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM)) {
                            element_type_name = "UpdatedLagrangianPQ";
                            KRATOS_ERROR_IF(IsAxisSymmetry) << "PQMPM is not implemented for axisymmetric elements yet." << std::endl;
                        }
                    }

                    // Get new element
                    const Element& new_element = KratosComponents<Element>::Get(element_type_name);

                    // Loop over the material points that fall in each grid element
                    unsigned int new_element_id = 0;
                    for (unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++)
                    {
                        mp_volume[0] = int_volumes[PointNumber];
                        mp_mass[0] = int_volumes[PointNumber]*density;

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
                        if (IsAxisSymmetry) {
                            mp_volume[0] *= 2.0 * Globals::Pi * xg[0][0];
                            mp_mass[0] = mp_volume[0] * density;
                        }

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
    template<SizeType TDimension>
    void GenerateMaterialPointCondition(ModelPart& rBackgroundGridModelPart,
                                            ModelPart& rInitialModelPart,
                                            ModelPart& rMPMModelPart){
        // Initialize zero the variables needed
        std::vector<array_1d<double, 3>> mpc_xg = { ZeroVector(3) };
        array_1d<double,3> mpc_normal = ZeroVector(3);
        array_1d<double,3> mpc_displacement = ZeroVector(3);
        array_1d<double,3> mpc_imposed_displacement = ZeroVector(3);
        array_1d<double,3> mpc_velocity = ZeroVector(3);
        array_1d<double,3> mpc_imposed_velocity = ZeroVector(3);
        array_1d<double,3> mpc_acceleration = ZeroVector(3);
        array_1d<double,3> mpc_imposed_acceleration = ZeroVector(3);
        std::vector<array_1d<double, 3>> mpc_contact_force = { ZeroVector(3) };
        array_1d<double, 3 > point_load = ZeroVector(3);;

        std::vector<double> mpc_area(1);
        std::vector<double> mpc_penalty_factor(1);


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

        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        // SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);

        // Loop over the submodelpart of rBackgroundGridModelPart
        for (ModelPart::SubModelPartIterator submodelpart_it = rBackgroundGridModelPart.SubModelPartsBegin();
                submodelpart_it != rBackgroundGridModelPart.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart& submodelpart = *submodelpart_it;

            // For submodelpart without condition, exit
            if (submodelpart.NumberOfConditions() != 0){

                std::string submodelpart_name = submodelpart.Name();

                // For regular conditions: straight copy all conditions
                if (!submodelpart.ConditionsBegin()->Is(BOUNDARY)){
                    if (submodelpart.NodesBegin()->Is(SLIP)){
                        // Do nothing, this is a slip condition applied directly
                        // to the background grid nodes.
                        // Check 'apply_mpm_slip_boundary_process.py'
                    }
                    else {
                        rMPMModelPart.CreateSubModelPart(submodelpart_name);
                        rMPMModelPart.SetConditions(submodelpart.pConditions());
                    }
                }
                // For boundary conditions: create particle conditions for all the necessary conditions
                else{
                    // NOTE: To create Particle Condition, we consider the position of integration point
                    // Loop over the conditions of submodelpart and generate mpm condition to be appended to the rMPMModelPart
                    rMPMModelPart.CreateSubModelPart(submodelpart_name);
                    for (ModelPart::ConditionIterator i = submodelpart.ConditionsBegin();
                            i != submodelpart.ConditionsEnd(); i++)
                    {
                        Properties::Pointer properties = i->pGetProperties();

                        // Flag whether condition is Neumann or Dirichlet
                        const bool is_neumann_condition = i->GetValue(MPC_IS_NEUMANN);
                        const int boundary_condition_type = i->GetValue(MPC_BOUNDARY_CONDITION_TYPE);

                        // Check number of particles per condition to be created
                        unsigned int particles_per_condition = 0; // Default zero
                        if (i->Has( PARTICLES_PER_CONDITION )){
                            particles_per_condition = i->GetValue(PARTICLES_PER_CONDITION);
                        }
                        else{
                            std::string warning_msg = "PARTICLES_PER_CONDITION is not specified, ";
                            warning_msg += "Only one particle is assumed";
                            KRATOS_WARNING("MPMParticleGeneratorUtility") << "WARNING: " << warning_msg << std::endl;
                        }

                        // Get condition variables:
                        // Normal vector
                        if (i->Has(NORMAL)) mpc_normal = i->GetValue(NORMAL);
                        ParticleMechanicsMathUtilities<double>::Normalize(mpc_normal);

                        // Get shape_function_values from defined particle_per_condition
                        const Geometry< Node < 3 > >& r_geometry = i->GetGeometry(); // current condition's geometry
                        const GeometryData::KratosGeometryType geo_type = r_geometry.GetGeometryType();
                        Matrix shape_functions_values;

                        // Get geometry and dimension of the background grid
                        const std::size_t domain_size = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];
                        

                        // Get integration method and shape function values
                        IntegrationMethod int_method = GeometryData::GI_GAUSS_1;
                        bool is_equal_int_volumes = false;

                        std::string condition_type_name;


                        DetermineConditionIntegrationMethodAndShapeFunctionValues(r_geometry, particles_per_condition,
                            int_method, shape_functions_values, is_equal_int_volumes);

                        

                        // Number of integration point per condition
                        const unsigned int integration_point_per_conditions = shape_functions_values.size1();
                        Vector int_volumes (integration_point_per_conditions);

                        if (is_equal_int_volumes) {
                            for (size_t j = 0; j < integration_point_per_conditions; ++j)  int_volumes[j] = r_geometry.Area() / integration_point_per_conditions;
                        }
                        else  GetIntegrationPointArea(r_geometry, int_method, int_volumes);

                        // Check condition variables
                        if (i->Has(DISPLACEMENT))
                            mpc_imposed_displacement = i->GetValue(DISPLACEMENT);
                        if (i->Has(VELOCITY))
                            mpc_imposed_velocity = i->GetValue(VELOCITY);
                        if (i->Has(ACCELERATION))
                            mpc_imposed_acceleration = i->GetValue(ACCELERATION);
                        if (i->Has(PENALTY_FACTOR))
                            mpc_penalty_factor[0] = i->GetValue(PENALTY_FACTOR);

                            const bool is_slip = i->Is(SLIP);
                            const bool is_contact = i->Is(CONTACT);
                            const bool is_interface = i->Is(INTERFACE);
                            const bool flip_normal_direction = i->Is(MODIFIED);

                        // If dirichlet boundary or coupling interface
                        if (!is_neumann_condition){
                            if(!is_interface){
                                condition_type_name = "MPMParticlePenaltyDirichletCondition";
                            }
                            else{
                                    condition_type_name = "MPMParticlePenaltyCouplingInterfaceCondition";
                                }
                        }
                        else{
                            if( i->Has( POINT_LOAD ) ){
                                point_load = i->GetValue( POINT_LOAD );
                                condition_type_name = "MPMParticlePointLoadCondition";
                            }
                            else{
                                KRATOS_ERROR << "Particle line load / surface load condition is not yet implemented." << std::endl;
                            }
                                
                        }

                        // Get new condition
                        const Condition& new_condition = KratosComponents<Condition>::Get(condition_type_name);

                        // Check Normal direction
                        if (flip_normal_direction) mpc_normal *= -1.0;

                        unsigned int new_condition_id = 0;

                       

                        // Create Particle Point Load Condition
                        if (condition_type_name == "MPMParticlePointLoadCondition" ){
                            // create point load condition
                            mpc_area[0] = 1;

                            // Create new material point condition
                            new_condition_id = last_condition_id + 1;
                        
                            mpc_xg[0].clear();
                            
                            for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++){
                                mpc_xg[0][dimension] = r_geometry[0].Coordinates()[dimension];
                            }
                            
                            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                            Element::Pointer pelem;
                            Vector N;
                            

                            // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                            bool is_found = SearchStructure.FindPointOnMesh(mpc_xg[0], N, pelem, result_begin);
                            if (!is_found) KRATOS_WARNING("MPM particle generator utility") << "::search failed." << std::endl;

                            auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                                pelem->pGetGeometry(), mpc_xg[0],
                                mpc_area[0]);

                            Condition::Pointer p_condition = new_condition.Create(
                                new_condition_id, p_new_geometry, properties);


                            ProcessInfo process_info = ProcessInfo();

                            // Setting particle condition's initial condition
                            p_condition->SetValuesOnIntegrationPoints(MPC_COORD, mpc_xg , process_info);
                            p_condition->SetValuesOnIntegrationPoints(MPC_AREA, mpc_area, process_info);
                            p_condition->SetValuesOnIntegrationPoints(POINT_LOAD, { point_load }, process_info);
                            

                            // Add the MP Condition to the model part
                            rMPMModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);

                        }
                        // Loop over the conditions to create inner particle condition (except point load condition)
                        else{
                            for ( unsigned int point_number = 0; point_number < integration_point_per_conditions; point_number++ )
                            {
                                mpc_area[0] = int_volumes[point_number];

                                mpc_xg[0].clear();

                                // Loop over the nodes of the grid condition
                                for (unsigned int dimension = 0; dimension < r_geometry.WorkingSpaceDimension(); dimension++){
                                    for ( unsigned int j = 0; j < r_geometry.size(); j ++){
                                        mpc_xg[0][dimension] = mpc_xg[0][dimension] + shape_functions_values(point_number, j) * r_geometry[j].Coordinates()[dimension];
                                    }
                                }

                                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                                Element::Pointer pelem;
                                Vector N;

                                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                                bool is_found = SearchStructure.FindPointOnMesh(mpc_xg[0], N, pelem, result_begin);
                                if (!is_found) KRATOS_WARNING("MPM particle generator utility") << "::MPC search failed." << std::endl;

                                pelem->Set(ACTIVE);
                                auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                                    pelem->pGetGeometry(), mpc_xg[0],
                                    mpc_area[0]);

                                // Create new material point condition
                                new_condition_id = last_condition_id + point_number;
                                Condition::Pointer p_condition = new_condition.Create(new_condition_id, p_new_geometry, properties);



                                ProcessInfo process_info = ProcessInfo();

                                // Setting particle condition's initial condition
                                // TODO: If any variable is added or remove here, please add and remove also at the second loop below
                                //p_condition->SetValuesOnIntegrationPoints(MPC_CONDITION_ID, mpc_condition_id, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_COORD, mpc_xg , process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_AREA,  mpc_area  , process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_NORMAL, { mpc_normal }, process_info);

                                
                                p_condition->SetValuesOnIntegrationPoints(MPC_DISPLACEMENT, { mpc_displacement }, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, { mpc_imposed_displacement }, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_VELOCITY, { mpc_velocity }, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_VELOCITY, { mpc_imposed_velocity }, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_ACCELERATION, { mpc_acceleration }, process_info);
                                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, { mpc_imposed_acceleration }, process_info);

                                if (boundary_condition_type == 1)
                                {
                                    p_condition->SetValuesOnIntegrationPoints(PENALTY_FACTOR, mpc_penalty_factor , process_info);
                                }
                                    

                                if (is_slip)
                                    p_condition->Set(SLIP);
                                if (is_contact)
                                    p_condition->Set(CONTACT);
                                if (is_interface)
                                {
                                    p_condition->Set(INTERFACE);
                                    p_condition->SetValuesOnIntegrationPoints(MPC_CONTACT_FORCE,  mpc_contact_force , process_info);
                                }

                                // Add the MP Condition to the model part
                                rMPMModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                            }

                            last_condition_id += integration_point_per_conditions;

                        }
                            
                        
                    }
                }
            }
        }
    }

}; // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos

#endif // KRATOS_MPM_PARTICLE_GENERATOR_UTILITY



// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//


#include "cr_beam_adjoint_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"


namespace Kratos
{

    CrBeamAdjointElement3D2N::CrBeamAdjointElement3D2N(IndexType NewId,
        GeometryType::Pointer pGeometry)
        : CrBeamElementLinear3D2N(NewId, pGeometry)
    {
    }

    CrBeamAdjointElement3D2N::CrBeamAdjointElement3D2N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : CrBeamElementLinear3D2N(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer CrBeamAdjointElement3D2N::Create(IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        const GeometryType& rGeom = this->GetGeometry();
        return BaseType::Pointer(new CrBeamAdjointElement3D2N(
            NewId, rGeom.Create(rThisNodes), pProperties));
    }

    Element::Pointer CrBeamAdjointElement3D2N::Create(IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(
                new CrBeamAdjointElement3D2N(NewId, pGeom, pProperties));
        KRATOS_CATCH("")
    }

    CrBeamAdjointElement3D2N::~CrBeamAdjointElement3D2N() {}

    void CrBeamAdjointElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) {

        const int number_of_nodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const unsigned int local_size = number_of_nodes * dimension * 2;

        if (rResult.size() != local_size) rResult.resize(local_size);

        for (int i = 0; i < number_of_nodes; ++i)
        {
            int index = i * number_of_nodes * dimension;
            rResult[index] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X)
                .EquationId();
            rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y)
                .EquationId();
            rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z)
                .EquationId();

            rResult[index + 3] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_X)
                .EquationId();
            rResult[index + 4] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_Y)
                .EquationId();
            rResult[index + 5] = this->GetGeometry()[i].GetDof(ADJOINT_ROTATION_Z)
                .EquationId();
        }

    }

    void CrBeamAdjointElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) {

        const int number_of_nodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const unsigned int local_size = number_of_nodes * dimension * 2;

        if (rElementalDofList.size() != local_size) {
            rElementalDofList.resize(local_size);
        }

        for (int i = 0; i < number_of_nodes; ++i)
        {
            int index = i * number_of_nodes * dimension;
            rElementalDofList[index] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_DISPLACEMENT_X);
            rElementalDofList[index + 1] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_DISPLACEMENT_Y);
            rElementalDofList[index + 2] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_DISPLACEMENT_Z);

            rElementalDofList[index + 3] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_ROTATION_X);
            rElementalDofList[index + 4] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_ROTATION_Y);
            rElementalDofList[index + 5] = this->GetGeometry()[i]
                .pGetDof(ADJOINT_ROTATION_Z);
        }
    }

    double CrBeamAdjointElement3D2N::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
    {
        KRATOS_TRY;

        if ( this->GetProperties().Has(rDesignVariable) )
        {
            const double variable_value = this->GetProperties()[rDesignVariable];
            return variable_value;
        }
        else
            return 1.0;

        KRATOS_CATCH("")
    }

    double CrBeamAdjointElement3D2N::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
    {
        KRATOS_TRY;

        if(rDesignVariable == SHAPE_SENSITIVITY)
        {
            double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
            double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
            double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
            double L = sqrt(dx*dx + dy*dy + dz*dz);
            return L;
        }
        else
            return 1.0;

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        // define working variables
        Vector RHS_undist;
        Vector RHS_dist;
        ProcessInfo copy_of_process_info = rCurrentProcessInfo;

        // Compute RHS before disturbing
        this->CalculateRightHandSide(RHS_undist, copy_of_process_info);

        rOutput.resize(1,RHS_undist.size());

        // Get disturbance measure
        double delta = this->GetValue(DISTURBANCE_MEASURE);
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        if ( this->GetProperties().Has(rDesignVariable) )
        {
            // Save properties and its pointer
            Properties& r_global_property = this->GetProperties();
            Properties::Pointer p_global_properties = this->pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Disturb the design variable
            const double current_property_value = this->GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

            // Compute RHS after disturbance
            this->CalculateRightHandSide(RHS_dist, copy_of_process_info);

            rOutput.resize(1,RHS_dist.size());

            // Compute derivative of RHS w.r.t. design variable with finite differences
            noalias(RHS_dist) -= RHS_undist;
            RHS_dist /= delta;
            for(unsigned int i = 0; i < RHS_dist.size(); i++)
                rOutput(0, i) = RHS_dist[i];

            // Give element original properties back
            this->SetProperties(p_global_properties);
        }
        else
        {
            rOutput.clear();
        }

        KRATOS_CATCH("")

    }

    void CrBeamAdjointElement3D2N::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // define working variables
        Vector RHS_undist;
        Vector RHS_dist;
        ProcessInfo copy_of_process_info = rCurrentProcessInfo;

        // Get disturbance measure
        double delta = this->GetValue(DISTURBANCE_MEASURE);
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        if(rDesignVariable == SHAPE_SENSITIVITY)
        {
            const int number_of_nodes = GetGeometry().PointsNumber();
            const int dimension = this->GetGeometry().WorkingSpaceDimension();
            const int local_size = number_of_nodes * dimension * 2;
            unsigned int num_coord_dir = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);

            rOutput.resize(dimension * number_of_nodes, local_size);

            // compute RHS before disturbing
            this->CalculateRightHandSide(RHS_undist, copy_of_process_info);

            int index = 0;
            //TODO: look that this works also for parallel computing
            for(auto& node_i : this->GetGeometry())
            {
                for(std::size_t coord_dir_i = 0; coord_dir_i < num_coord_dir; coord_dir_i++)
                {
                    // disturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] += delta;

                    // compute RHS after disturbance
                    this->CalculateRightHandSide(RHS_dist, copy_of_process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    noalias(RHS_dist) -= RHS_undist;
                    RHS_dist /= delta;
                    for(unsigned int i = 0; i < RHS_dist.size(); i++)
                        rOutput( (coord_dir_i + index*dimension), i) = RHS_dist[i]; 

                    // Reset pertubed vector
                    RHS_dist = Vector(0);

                    // undisturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] -= delta;
                }
                index++;
            }// end loop over element nodes
        }
        else
            KRATOS_ERROR << "Unsupported design variable!" << std::endl;

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if(rVariable == STRESS_ON_GP || rVariable == STRESS_ON_NODE)
        {
            TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

            std::vector< array_1d<double, 3 > > stress_vector;
            int direction_1 = 0;
            bool stress_is_moment = true;
  
            switch (traced_stress_type)  
            { 
                case MX:
                {
                    direction_1 = 0; 
                    break;
                }
                case MY:
                {
                    direction_1 = 1; 
                    break; 
                }
                case MZ:
                {
                    direction_1 = 2; 
                    break; 
                }
                case FX:
                {
                    direction_1 = 0; 
                    stress_is_moment = false;
                    break;
                }	
                case FY:
                {
                    direction_1 = 1; 
                    stress_is_moment = false;
                    break;
                }
                case FZ:
                {
                    direction_1 = 2; 
                    stress_is_moment = false;
                    break;
                }
                default:
                    KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;  
            }

            if(stress_is_moment)
                CrBeamElementLinear3D2N::GetValueOnIntegrationPoints(MOMENT, stress_vector, rCurrentProcessInfo);
            else
                CrBeamElementLinear3D2N::GetValueOnIntegrationPoints(FORCE, stress_vector, rCurrentProcessInfo);
    
            if(rVariable == STRESS_ON_GP)
            {
                const unsigned int&  GP_num = GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);

                rOutput.resize(GP_num);
                for(unsigned int i = 0; i < GP_num ; i++)
                {
                    rOutput(i) = stress_vector[i][direction_1];
                }
            }
            else if(rVariable == STRESS_ON_NODE)
            {
                rOutput.resize(2);
                rOutput(0) = 2 * stress_vector[0][direction_1] - stress_vector[1][direction_1];
                rOutput(1) = 2 * stress_vector[2][direction_1] - stress_vector[1][direction_1];
            }
        }
        else
        {
            rOutput.resize(3);
            rOutput.clear();
        }

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::Calculate(const Variable<Matrix >& rVariable,
                           Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
    {
           KRATOS_TRY;


        if(rVariable == STRESS_DISP_DERIV_ON_GP)
        {
               this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if(rVariable == STRESS_DISP_DERIV_ON_NODE)
        {
            this->CalculateStressDisplacementDerivative(STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
        else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
        {
            std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

            if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(design_varible_name);
                this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
                this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
            }
        }
        else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_NODE)
        {
            std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

            if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(design_varible_name);
                this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
                this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
            }
        }
           else
        {
            rOutput.clear();
        }

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int num_nodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const int num_dofs = num_nodes * dimension * 2;
        ProcessInfo copy_process_info = rCurrentProcessInfo;

        Vector initial_state_variables;
        initial_state_variables.resize(num_dofs);
        Vector stress_derivatives_vector;

        Vector dummy;
        this->Calculate(rStressVariable, dummy, rCurrentProcessInfo);
        rOutput.resize(num_dofs, dummy.size() );
        rOutput.clear();
        
        // Built vector of variables containing the DOF-variables of the primal problem 
        std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>> primal_solution_variable_list;
        primal_solution_variable_list.push_back(DISPLACEMENT_X);       
        primal_solution_variable_list.push_back(DISPLACEMENT_Y);       
        primal_solution_variable_list.push_back(DISPLACEMENT_Z);       
        primal_solution_variable_list.push_back(ROTATION_X);       
        primal_solution_variable_list.push_back(ROTATION_Y);       
        primal_solution_variable_list.push_back(ROTATION_Z);       
        
        // Concept A: Analytic apporoch ###################################################

        for (int i = 0; i < num_nodes; i++) 
        {	
            int index = i * dimension * 2;
            for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
            {
                initial_state_variables[index + j] = this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
            }
        }

        for (int i = 0; i < num_nodes; i++) 
        {	
            int index = i * dimension * 2;
            for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
            {
                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 1.0;
                
                this->Calculate(rStressVariable, stress_derivatives_vector, rCurrentProcessInfo);
                
                for(unsigned int k = 0; k < stress_derivatives_vector.size(); k++)
                    rOutput(index+j, k) = stress_derivatives_vector[k];
                
                stress_derivatives_vector.clear();
                
                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
            }
        }
       
        for (int i = 0; i < num_nodes; i++) 
        {	
            int index = i * dimension * 2;
            for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_state_variables[index + j];
        }

        // Concept B: Derive by finite differences ###################################################

        /*this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);
        Vector stress_vector_undist;
        Vector stress_vector_dist;
        double initial_value_of_state_variable = 0.0;

        // Get disturbance measure
        //double dist_measure = this->GetValue(DISTURBANCE_MEASURE);

        unsigned int size_stress_vec = stress_vector_undist.size();	
            
        rOutput.resize(num_dofs, size_stress_vec);
        rOutput.clear();
        int index = 0;
        for (int i = 0; i < num_nodes; i++) 
        {	
            for(unsigned int j = 0; j < primal_solution_variable_list.size(); j++)
            {
                initial_value_of_state_variable = this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
                
                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_value_of_state_variable + dist_measure;
                
                this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);
            
                for(unsigned int k = 0; k < size_stress_vec; k++)
                {
                    stress_vector_dist[k] -= stress_vector_undist[k];
                    stress_vector_dist[k] /= dist_measure;
                    rOutput(index,k) = stress_vector_dist[k];
                }

                this->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_value_of_state_variable;

                stress_vector_dist.clear();
                index++;
            }
        }*/

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                        const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
         KRATOS_TRY;

        // Define working variables
        Vector stress_vector_undist;
        Vector stress_vector_dist;
        ProcessInfo copy_process_info = rCurrentProcessInfo;

        // Get disturbance measure
        double delta= this->GetValue(DISTURBANCE_MEASURE);
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        // Compute stress before disturbance
        this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        const int stress_vector_size = stress_vector_undist.size();
        rOutput.resize(1, stress_vector_size);

        if( this->GetProperties().Has(rDesignVariable) )
        {
            // Save properties and its pointer
            Properties& r_global_property = this->GetProperties();
            Properties::Pointer p_global_properties = this->pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Disturb the design variable
            const double current_property_value = this->GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

            // Compute stress on GP after disturbance
            this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

            // Compute derivative of stress w.r.t. design variable with finite differences
            noalias(stress_vector_dist)  -= stress_vector_undist;
            stress_vector_dist  /= delta;

            for(int j = 0; j < stress_vector_size; j++)
                rOutput(0, j) = stress_vector_dist[j];

            // Give element original properties back
            this->SetProperties(p_global_properties);
        }
        else
            rOutput.clear();

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // define working variables
        Vector stress_vector_undist;
        Vector stress_vector_dist;
        ProcessInfo copy_process_info = rCurrentProcessInfo;

        // Get disturbance measure
        double delta= this->GetValue(DISTURBANCE_MEASURE);
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        if(rDesignVariable == SHAPE_SENSITIVITY)
        {
            const int number_of_nodes = GetGeometry().PointsNumber();
            const int dimension = this->GetGeometry().WorkingSpaceDimension();
            unsigned int num_coord_dir = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);

            // Compute stress on GP before disturbance
            this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

            const int stress_vector_size = stress_vector_undist.size();

            rOutput.resize(dimension * number_of_nodes, stress_vector_size);

            int index = 0;
            //TODO: look that this works also for parallel computing
            for(auto& node_i : this->GetGeometry())
            {
                for(std::size_t coord_dir_i = 0; coord_dir_i < num_coord_dir; coord_dir_i++)
                {
                    // Disturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] += delta;

                    // Compute stress on GP after disturbance
                    this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

                    // Compute derivative of stress w.r.t. design variable with finite differences
                    noalias(stress_vector_dist)  -= stress_vector_undist;
                    stress_vector_dist  /= delta;

                    for(int i = 0; i < stress_vector_size; i++)
                        rOutput( (coord_dir_i + index*dimension), i) = stress_vector_dist[i];

                    // Reset pertubed vector
                    stress_vector_dist = Vector(0);

                    // Undisturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] -= delta;
                }
                index++;
            }// end loop over element nodes
        }
        else
            KRATOS_ERROR << "Unsupported design variable!" << std::endl;

        KRATOS_CATCH("")
    }


    void CrBeamAdjointElement3D2N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                          std::vector<double>& rOutput,
                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if(this->Has(rVariable))
        {
            // Get result value for output
            double output_value = this->GetValue(rVariable);

            // Resize Output
            const unsigned int&  write_points_number = GetGeometry()
                .IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
            if (rOutput.size() != write_points_number)
            {
                rOutput.resize(write_points_number);
            }

            // Write scalar result value on all Gauss-Points
            for(unsigned int i = 0; i < write_points_number; ++i)
            {
                rOutput[i] = output_value;
            }
        }
        else
            KRATOS_ERROR << "Unsupported output variable." << std::endl;


        KRATOS_CATCH("")

    }

    void CrBeamAdjointElement3D2N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                         std::vector<double>& rValues,
                         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::GetValuesVector(Vector& rValues, int Step) {

        KRATOS_TRY
            const int number_of_nodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const unsigned int element_size = number_of_nodes * dimension * 2;

        if (rValues.size() != element_size) rValues.resize(element_size, false);

        for (int i = 0; i < number_of_nodes; ++i)
        {
            int index = i * dimension * 2;
            rValues[index] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X, Step);
            rValues[index + 1] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Y, Step);
            rValues[index + 2] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Z, Step);

            rValues[index + 3] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_ROTATION_X, Step);
            rValues[index + 4] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_ROTATION_Y, Step);
            rValues[index + 5] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ADJOINT_ROTATION_Z, Step);
        }
        KRATOS_CATCH("")
    }

    int CrBeamAdjointElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
        << "The beam element works only in 3D and with 2 noded elements" << "" << std::endl;

        // verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
        KRATOS_CHECK_VARIABLE_KEY(DENSITY);
        KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
        KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

        // check properties
        KRATOS_ERROR_IF(this->GetProperties().Has(CROSS_AREA) == false || this->GetProperties()[CROSS_AREA] == 0)
        << "CROSS_AREA not provided for this element" << this->Id() << std::endl;

        KRATOS_ERROR_IF(this->GetProperties().Has(YOUNG_MODULUS) == false || this->GetProperties()[YOUNG_MODULUS] == 0)
        << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT( this->GetProperties().Has(DENSITY) )
        << "DENSITY not provided for this element" << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT( this->GetProperties().Has(POISSON_RATIO) )
        << "POISSON_RATIO not provided for this element" << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT( this->GetProperties().Has(TORSIONAL_INERTIA) )
        << "TORSIONAL_INERTIA not provided for this element" << this->Id() << std::endl;
    
        KRATOS_ERROR_IF_NOT( this->GetProperties().Has(I22) )
        << "I22 not provided for this element" << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT( this->GetProperties().Has(I33) )
        << "I33 not provided for this element" << this->Id() << std::endl;

        // Check dofs
        GeometryType& r_geom = GetGeometry();
        for (unsigned int i = 0; i < r_geom.size(); i++)
        {
            auto& r_node = r_geom[i];

            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
        }

        return 0;

        KRATOS_CATCH("")
    }

    void CrBeamAdjointElement3D2N::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElementLinear3D2N);
    }

    void CrBeamAdjointElement3D2N::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElementLinear3D2N);
    }

} // namespace Kratos.



// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Martin Fusseder
//                   
//                   
//
#include "custom_elements/adjoint_elements/cr_beam_adjoint_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"



namespace Kratos
{

	CrBeamAdjointElement3D2N::CrBeamAdjointElement3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry, bool rLinear)
		: CrBeamElement3D2N(NewId, pGeometry, rLinear)
	{
		this->mIsLinearElement = rLinear;
	}

	CrBeamAdjointElement3D2N::CrBeamAdjointElement3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties, bool rLinear)
		: CrBeamElement3D2N(NewId, pGeometry, pProperties, rLinear)
	{
		this->mIsLinearElement = rLinear;
	}

	Element::Pointer CrBeamAdjointElement3D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new CrBeamAdjointElement3D2N(
			NewId, rGeom.Create(rThisNodes), pProperties, true));//this->mIsLinearElement));
	}

	Element::Pointer CrBeamAdjointElement3D2N::Create(IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
        return Element::Pointer(
                new CrBeamAdjointElement3D2N(NewId, pGeom, pProperties, true)); //this->mIsLinearElement));
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
		ProcessInfo testProcessInfo = rCurrentProcessInfo;

		// Compute RHS before disturbing
		this->CalculateRightHandSide(RHS_undist, testProcessInfo); 

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
			this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

			rOutput.resize(1,RHS_dist.size());

			// Compute derivative of RHS w.r.t. design variable with finite differences
			RHS_dist -= RHS_undist;
			RHS_dist /= delta;
			for(unsigned int i = 0; i < RHS_dist.size(); i++)
				rOutput(0, i) = RHS_dist[i];

            // Give element original properties back
            this->SetProperties(p_global_properties);	

			// Compute RHS again in order to ensure that changed member variables like mLHS get back their origin values
			this->CalculateRightHandSide(RHS_dist, testProcessInfo); 
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
		ProcessInfo testProcessInfo = rCurrentProcessInfo;
		
		// Get disturbance measure
        double delta = this->GetValue(DISTURBANCE_MEASURE); 
		double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
		delta *= correction_factor;

		if(rDesignVariable == SHAPE_SENSITIVITY) 
		{
			const int number_of_nodes = GetGeometry().PointsNumber();
			const int dimension = this->GetGeometry().WorkingSpaceDimension();
			const int local_size = number_of_nodes * dimension * 2;
 
			rOutput.resize(dimension * number_of_nodes, local_size);

			// compute RHS before disturbing
			this->CalculateRightHandSide(RHS_undist, testProcessInfo); 

            //TODO: look that this works also for parallel computing
			for(int j = 0; j < number_of_nodes; j++)
			{
				//begin: derive w.r.t. x-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].X0() += delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();

				// compute RHS after disturbance
				this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

				//compute derivative of RHS w.r.t. design variable with finite differences
				RHS_dist -= RHS_undist;
				RHS_dist /= delta;
				for(unsigned int i = 0; i < RHS_dist.size(); i++)  
					rOutput( (0 + j*dimension), i) = RHS_dist[i]; 

				// Reset pertubed vector
				RHS_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].X0() -= delta;

				//end: derive w.r.t. x-coordinate-----------------------------------------------------

				//begin: derive w.r.t. y-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].Y0() += delta;
                // Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();

				// compute RHS after disturbance
				this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

				//compute derivative of RHS w.r.t. design variable with finite differences
				RHS_dist -= RHS_undist;
				RHS_dist /= delta;
				for(unsigned int i = 0; i < RHS_dist.size(); i++) 
					rOutput((1 + j*dimension),i) = RHS_dist[i]; 

				// Reset pertubed vector
				RHS_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].Y0() -= delta;

				//end: derive w.r.t. y-coordinate-----------------------------------------------------

				//begin: derive w.r.t. z-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].Z0() += delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();

				// compute RHS after disturbance
				this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

				//compute derivative of RHS w.r.t. design variable with finite differences
				RHS_dist -= RHS_undist;
				RHS_dist /= delta;
				for(unsigned int i = 0; i < RHS_dist.size(); i++) 
					rOutput((2 + j*dimension),i) = RHS_dist[i];

				// Reset pertubed vector
				RHS_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].Z0() -= delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();

				// Compute RHS again in order to ensure that changed member variables like mLHS get back their origin values
				this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

				//end: derive w.r.t. z-coordinate-----------------------------------------------------

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

		const Variable<Vector> & rSTRESS_ON_GP =
           	  KratosComponents<Variable<Vector> >::Get("STRESS_ON_GP");  
		const Variable<Vector> & rSTRESS_ON_NODE =
           	  KratosComponents<Variable<Vector> >::Get("STRESS_ON_NODE");  	
		const Variable<std::string>& rTRACED_STRESS_TYPE =
           	  KratosComponents<Variable<std::string>>::Get("TRACED_STRESS_TYPE");	 

		if(rVariable == rSTRESS_ON_GP || rVariable == rSTRESS_ON_NODE)
		{
			std::string traced_stress_type = this->GetValue(rTRACED_STRESS_TYPE);

    		const char item_1 = traced_stress_type.at(0);
    		const char item_2 = traced_stress_type.at(1);

    		int direction_1 = 0;
    		std::vector< array_1d<double, 3 > > stress_vector;
  
    		if(item_1 == 'M') 
        		CrBeamElement3D2N::GetValueOnIntegrationPoints(MOMENT, stress_vector, rCurrentProcessInfo);
    		else if(item_1 == 'F') 
        		CrBeamElement3D2N::GetValueOnIntegrationPoints(FORCE, stress_vector, rCurrentProcessInfo);
    		else 
        		KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;  

    		if(item_2 == 'X')  
        		direction_1 = 0; 
    		else if(item_2 == 'Y')  
        		direction_1 = 1; 
    		else if(item_2 == 'Z')  
        		direction_1 = 2;   
    		else 
        		KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;              

			if(rVariable == rSTRESS_ON_GP)
			{
				const unsigned int&  GP_num = GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);	

    			rOutput.resize(GP_num);   
    			for(unsigned int i = 0; i < GP_num ; i++)
    			{
        			rOutput(i) = stress_vector[i][direction_1];
    			}
			}
			else if(rVariable == rSTRESS_ON_NODE)
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

		const Variable<Matrix> & rSTRESS_DISP_DERIV_ON_GP =
           	  KratosComponents<Variable<Matrix> >::Get("STRESS_DISP_DERIV_ON_GP");  
		const Variable<Matrix> & rSTRESS_DISP_DERIV_ON_NODE =
           	  KratosComponents<Variable<Matrix> >::Get("STRESS_DISP_DERIV_ON_NODE"); 
		const Variable<Matrix> & rSTRESS_DV_DERIV_ON_GP =
           	  KratosComponents<Variable<Matrix> >::Get("STRESS_DV_DERIV_ON_GP");  
		const Variable<Matrix> & rSTRESS_DV_DERIV_ON_NODE =
           	  KratosComponents<Variable<Matrix> >::Get("STRESS_DV_DERIV_ON_NODE");   
		const Variable<Vector>& rSTRESS_ON_GP =
        	  KratosComponents<Variable<Vector>>::Get("STRESS_ON_GP");
		const Variable<Vector>& rSTRESS_ON_NODE =
        	  KratosComponents<Variable<Vector>>::Get("STRESS_ON_NODE");		 

		if(rVariable == rSTRESS_DISP_DERIV_ON_GP) 
		{
       		this->CalculateStressDisplacementDerivative(rSTRESS_ON_GP, rOutput, rCurrentProcessInfo);
    	}
		else if(rVariable == rSTRESS_DISP_DERIV_ON_NODE)
		{
			this->CalculateStressDisplacementDerivative(rSTRESS_ON_NODE, rOutput, rCurrentProcessInfo);
		}
    	else if(rVariable == rSTRESS_DV_DERIV_ON_GP)
    	{
        	const Variable<std::string> & rDESIGN_VARIABLE_NAME =
           		KratosComponents<Variable<std::string> >::Get("DESIGN_VARIABLE_NAME");
        	std::string design_varible_name = this->GetValue( rDESIGN_VARIABLE_NAME );	

        	if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        	{
            	const Variable<double>& r_variable =
                	KratosComponents<Variable<double>>::Get(design_varible_name);
            	this->CalculateStressDesignVariableDerivative(r_variable, rSTRESS_ON_GP, rOutput, rCurrentProcessInfo);
        	}
        	else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        	{
            	const Variable<array_1d<double, 3>>& r_variable =
                	KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            	this->CalculateStressDesignVariableDerivative(r_variable, rSTRESS_ON_GP, rOutput, rCurrentProcessInfo);    
        	}      
    	}
		else if(rVariable == rSTRESS_DV_DERIV_ON_NODE)
		{
			const Variable<std::string> & rDESIGN_VARIABLE_NAME =
           		KratosComponents<Variable<std::string> >::Get("DESIGN_VARIABLE_NAME");
        	std::string design_varible_name = this->GetValue( rDESIGN_VARIABLE_NAME );	

        	if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        	{
            	const Variable<double>& r_variable =
                	KratosComponents<Variable<double>>::Get(design_varible_name);
            	this->CalculateStressDesignVariableDerivative(r_variable, rSTRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        	}
        	else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        	{
            	const Variable<array_1d<double, 3>>& r_variable =
                	KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            	this->CalculateStressDesignVariableDerivative(r_variable, rSTRESS_ON_NODE, rOutput, rCurrentProcessInfo);    
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
    	Vector stress_vector_undist;
    	Vector stress_vector_dist;
		ProcessInfo copy_process_info = rCurrentProcessInfo;	
		
		// Get disturbance measure
        double dist_measure = this->GetValue(DISTURBANCE_MEASURE); 
		//TODO: is here a correction of delta necessary???

    	this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);
	
		DofsVectorType element_dof_list;
    	CrBeamElement3D2N::GetDofList(element_dof_list, copy_process_info);
			
		unsigned int size_stress_vec = stress_vector_undist.size();
    	rOutput.resize(num_dofs, size_stress_vec);
    	for(int i = 0; i < num_dofs; i++)
    	{
        	element_dof_list[i]->GetSolutionStepValue() += dist_measure;

    		this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);
		
        	for(unsigned int j = 0; j < size_stress_vec; j++)
        	{
            	stress_vector_dist[j] -= stress_vector_undist[j];
            	stress_vector_dist[j] /= dist_measure;
            	rOutput(i,j) = stress_vector_dist[j];
        	}   

        	element_dof_list[i]->GetSolutionStepValue() -= dist_measure;
        	stress_vector_dist.clear();
    	}

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
		Matrix dummy_LHS;
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

			// Update stiffness matrix
			this->CalculateLeftHandSide(dummy_LHS, copy_process_info); 

			// Compute stress on GP after disturbance
		    this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

			// Compute derivative of stress w.r.t. design variable with finite differences
			stress_vector_dist  -= stress_vector_undist;
			stress_vector_dist  /= delta;

			for(int j = 0; j < stress_vector_size; j++)
			    rOutput(0, j) = stress_vector_dist[j];
		
            // Give element original properties back
            this->SetProperties(p_global_properties);

			// Update stiffness matrix
			this->CalculateLeftHandSide(dummy_LHS, copy_process_info); 
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
		Matrix dummy_LHS;
		ProcessInfo copy_process_info = rCurrentProcessInfo;	
		
		// Get disturbance measure
        double delta= this->GetValue(DISTURBANCE_MEASURE);
		double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
		delta *= correction_factor;	 	

		if(rDesignVariable == SHAPE_SENSITIVITY) 
		{
			const int number_of_nodes = GetGeometry().PointsNumber();
			const int dimension = this->GetGeometry().WorkingSpaceDimension();

			// Compute stress on GP before disturbance
	    	this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

			const int stress_vector_size = stress_vector_undist.size();
 
			rOutput.resize(dimension * number_of_nodes, stress_vector_size);
    
        	//TODO: look that this works also for parallel computing
			for(int j = 0; j < number_of_nodes; j++)
			{
				//begin: derive w.r.t. x-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].X0() += delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS(); 
				// Update stiffness matrix
				this->CalculateLeftHandSide(dummy_LHS, copy_process_info); 

				// Compute stress on GP after disturbance
				this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

				// Compute derivative of stress w.r.t. design variable with finite differences
				stress_vector_dist  -= stress_vector_undist;
				stress_vector_dist  /= delta;

				for(int i = 0; i < stress_vector_size; i++)
					rOutput( (0 + j*dimension), i) = stress_vector_dist[i]; 

				// Reset pertubed vector
				stress_vector_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].X0() -= delta;

				//end: derive w.r.t. x-coordinate-----------------------------------------------------

				//begin: derive w.r.t. y-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].Y0() += delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS(); 
				// Update stiffness matrix
				this->CalculateLeftHandSide(dummy_LHS, copy_process_info); 

				// Compute stress on GP after disturbance
				this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

				// Compute derivative of stress w.r.t. design variable with finite differences
				stress_vector_dist  -= stress_vector_undist;
				stress_vector_dist  /= delta;

				for(int i = 0; i < stress_vector_size; i++)
					rOutput((1 + j*dimension),i) = stress_vector_dist[i]; 

				// Reset pertubed vector
				stress_vector_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].Y0() -= delta;

				//end: derive w.r.t. y-coordinate-----------------------------------------------------

				//begin: derive w.r.t. z-coordinate---------------------------------------------------

				// disturb the design variable
				this->GetGeometry()[j].Z0() += delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();
				// Update stiffness matrix
				this->CalculateLeftHandSide(dummy_LHS, copy_process_info);  

				// Compute stress on GP after disturbance
				this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

				// Compute derivative of stress w.r.t. design variable with finite differences
				stress_vector_dist  -= stress_vector_undist;
				stress_vector_dist  /= delta;

				for(int i = 0; i < stress_vector_size; i++)
					rOutput((2 + j*dimension),i) = stress_vector_dist[i];

				// Reset pertubed vector
				stress_vector_dist = Vector(0);

				// undisturb the design variable
				this->GetGeometry()[j].Z0() -= delta;
				// Update CS and transformation matrix after geometry change
				this->CalculateInitialLocalCS();
				// Update stiffness matrix
				this->CalculateLeftHandSide(dummy_LHS, copy_process_info);  

				//end: derive w.r.t. z-coordinate-----------------------------------------------------

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

		if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
		{
			KRATOS_ERROR <<
				"The beam element works only in 3D and with 2 noded elements" << ""
			<< std::endl;
		}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0) {
			KRATOS_ERROR <<
				"VELOCITY has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (DISPLACEMENT.Key() == 0) {
			KRATOS_ERROR <<
				"DISPLACEMENT has Key zero! (check if the application is correctly registered"<< ""
				<< std::endl;
		}
		if (ACCELERATION.Key() == 0) {
			KRATOS_ERROR <<
				"ACCELERATION has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (DENSITY.Key() == 0) {
			KRATOS_ERROR <<
				"DENSITY has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (CROSS_AREA.Key() == 0) {
			KRATOS_ERROR <<
				"CROSS_AREA has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}

		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_ERROR << "CROSS_AREA not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_ERROR << "YOUNG_MODULUS not provided for this element" << this->Id()
				<< std::endl;
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_ERROR << "DENSITY not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(POISSON_RATIO) == false)
		{
			KRATOS_ERROR << "POISSON_RATIO not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(TORSIONAL_INERTIA) == false)
		{
			KRATOS_ERROR << "TORSIONAL_INERTIA not provided for this element" << this->Id()
				<< std::endl;
		}
		if (this->GetProperties().Has(I22) == false)
		{
			KRATOS_ERROR << "I22 not provided for this element" << this->Id()
				<< std::endl;
		}
		if (this->GetProperties().Has(I33) == false)
		{
			KRATOS_ERROR << "I33 not provided for this element" << this->Id()
				<< std::endl;
		}		
		

		//##################################################################################################
		// Check for specific sensitivity analysis stuff
		//##################################################################################################
		if (ADJOINT_DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_DISPLACEMENT Key is 0. "
                    "Check if the application was correctly registered.","");

		if (ADJOINT_ROTATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,
                    "ADJOINT_ROTATION Key is 0. "
                    "Check if the application was correctly registered.","");

		 // Check if the nodes have adjoint dofs.
        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_DISPLACEMENT_X) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_DISPLACEMENT_Y) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_DISPLACEMENT component degree of freedom on node ",
                        this->GetGeometry()[iNode].Id());

			if (this->GetGeometry()[iNode].HasDofFor(ADJOINT_ROTATION_X) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_ROTATION_Y) == false
                    || this->GetGeometry()[iNode].HasDofFor(ADJOINT_ROTATION_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing ADJOINT_ROTATION component degree of freedom on node ",
                        this->GetGeometry()[iNode].Id());	

			if (this->GetGeometry()[iNode].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,
                        "missing DISPLACEMENT variable on solution step data for node ",
                        this->GetGeometry()[iNode].Id());					
        }

		return 0;

		KRATOS_CATCH("")
	}

	/*std::string CrBeamAdjointElement3D2N::Info() const
    {
		return "CrBeamAdjointElement3D2N";
		//fusseder TODO: seperate between linear and nonliner case!!!!
    }*/

	void CrBeamAdjointElement3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElement3D2N);
	}

	void CrBeamAdjointElement3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElement3D2N);
	}

} // namespace Kratos.



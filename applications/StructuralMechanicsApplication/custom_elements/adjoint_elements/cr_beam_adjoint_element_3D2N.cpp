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
			NewId, rGeom.Create(rThisNodes), pProperties, this->mIsLinearElement));
	}

	Element::Pointer CrBeamAdjointElement3D2N::Create(IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
        return Element::Pointer(
                new CrBeamAdjointElement3D2N(NewId, pGeom, pProperties, this->mIsLinearElement));
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

		//TODO: reduce code duplication

		if (rDesignVariable == IT || rDesignVariable == IY || rDesignVariable == IZ)
        { 
			// Save properties and its pointer
            Properties& r_global_property = this->GetProperties(); 
            Properties::Pointer p_global_properties = this->pGetProperties(); 

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

			// Check which entry of the inertia vector shall be treated as design variable
			Vector& inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
			double& design_variable = rDesignVariable==IT ? inertia[0] : rDesignVariable==IY ? inertia[1] : inertia[2];

			// Get disturbance measure
            double delta = design_variable * 0.001; // TODO: get this from outside!

			design_variable += delta;

			// Compute RHS after disturbance
			this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

			// Compute derivative of RHS w.r.t. design variable with finite differences
			RHS_dist -= RHS_undist;
			RHS_dist /= delta;
			for(unsigned int i = 0; i < RHS_dist.size(); i++)
				rOutput(0, i) = RHS_dist[i];
		
            // Give element original properties back
            this->SetProperties(p_global_properties);
        }
		else if ( this->GetProperties().Has(rDesignVariable) ) 
		{
       
			// Save properties and its pointer
            Properties& r_global_property = this->GetProperties(); 
            Properties::Pointer p_global_properties = this->pGetProperties(); 

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Get disturbance measure
            double delta = r_global_property[rDesignVariable] * 0.001; // TODO: get this from outside!

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
		double delta = 1e-3;	//TODO: get this from outside!

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
    	double dist_measure = 1e-6; //------------------>TODO: get this from outside
    	ProcessInfo copy_process_info = rCurrentProcessInfo;		

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
		ProcessInfo copy_process_info = rCurrentProcessInfo;	

        // Compute stress before disturbance
		this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

		const int stress_vector_size = stress_vector_undist.size();
        rOutput.resize(1, stress_vector_size);

		if (rDesignVariable == IT || rDesignVariable == IY || rDesignVariable == IZ)
        { 
			// Save properties and its pointer
            Properties& r_global_property = this->GetProperties(); 
            Properties::Pointer p_global_properties = this->pGetProperties(); 

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

			// Check which entry of the inertia vector shall be treated as design variable
			Vector& inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
			double& design_variable = rDesignVariable==IT ? inertia[0] : rDesignVariable==IY ? inertia[1] : inertia[2];
			
			// Get disturbance measure
            double delta = 1e-05;//design_variable * 0.01; // TODO: get this from outside!

			// Disturb Design Variable
			design_variable += delta;

			// Update stiffness matrix
			Matrix dummy_LHS;
			this->CalculateLeftHandSide(dummy_LHS, copy_process_info); 

			// Compute stress after disturbance
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
		else if( this->GetProperties().Has(rDesignVariable) ) 
		{
			// Save properties and its pointer
            Properties& r_global_property = this->GetProperties(); 
            Properties::Pointer p_global_properties = this->pGetProperties(); 

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Get disturbance measure
            double delta = r_global_property[rDesignVariable] * 0.001; // TODO: get this from outside!

			// Disturb the design variable
			const double current_property_value = this->GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

			// Update stiffness matrix
			Matrix dummy_LHS;
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
		ProcessInfo copy_process_info = rCurrentProcessInfo;	
		double delta = 1e-5;	//TODO: get this from outside!

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
				Matrix dummy_LHS;
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

	/*void CrBeamAdjointElement3D2N::Initialize() {

		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		if (this->mIterationCount == 0)
		{
			this->mNodalForces = ZeroVector(local_size * 2);
		}
		KRATOS_CATCH("")
	}*/

	/*Matrix CrBeamAdjointElement3D2N::CreateElementStiffnessMatrix_Material() {

		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension * 2;

		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = this->CalculateShearModulus();
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();

		Vector inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y) == true) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z) == true) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}
		const double Psi_y = this->CalculatePsi(Iy, Az);
		const double Psi_z = this->CalculatePsi(Iz, Ay);



		Matrix LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);
		const double L3 = L*L*L;
		const double L2 = L*L;


		LocalStiffnessMatrix(0, 0) = E*A / L;
		LocalStiffnessMatrix(6, 0) = -1.0 * LocalStiffnessMatrix(0, 0);
		LocalStiffnessMatrix(0, 6) = LocalStiffnessMatrix(6, 0);
		LocalStiffnessMatrix(6, 6) = LocalStiffnessMatrix(0, 0);

		LocalStiffnessMatrix(1, 1) = 12.0 * E * Iz * Psi_z / L3;
		LocalStiffnessMatrix(1, 7) = -1.0 * LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(1, 5) = 6.0 * E * Iz * Psi_z / L2;
		LocalStiffnessMatrix(1, 11) = LocalStiffnessMatrix(1, 5);

		LocalStiffnessMatrix(2, 2) = 12.0 * E *Iy * Psi_y / L3;
		LocalStiffnessMatrix(2, 8) = -1.0 * LocalStiffnessMatrix(2, 2);
		LocalStiffnessMatrix(2, 4) = -6.0 * E *Iy *Psi_y / L2;
		LocalStiffnessMatrix(2, 10) = LocalStiffnessMatrix(2, 4);

		LocalStiffnessMatrix(4, 2) = LocalStiffnessMatrix(2, 4);
		LocalStiffnessMatrix(5, 1) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(3, 3) = G*J / L;
		LocalStiffnessMatrix(4, 4) = E*Iy*(3.0 * Psi_y + 1.0) / L;
		LocalStiffnessMatrix(5, 5) = E*Iz*(3.0 * Psi_z + 1.0) / L;
		LocalStiffnessMatrix(4, 8) = -1.0 * LocalStiffnessMatrix(4, 2);
		LocalStiffnessMatrix(5, 7) = -1.0 * LocalStiffnessMatrix(5, 1);
		LocalStiffnessMatrix(3, 9) = -1.0 * LocalStiffnessMatrix(3, 3);
		LocalStiffnessMatrix(4, 10) = E*Iy*(3.0 * Psi_y - 1) / L;
		LocalStiffnessMatrix(5, 11) = E*Iz*(3.0 * Psi_z - 1) / L;

		LocalStiffnessMatrix(7, 1) = LocalStiffnessMatrix(1, 7);
		LocalStiffnessMatrix(7, 5) = LocalStiffnessMatrix(5, 7);
		LocalStiffnessMatrix(7, 7) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(7, 11) = LocalStiffnessMatrix(7, 5);

		LocalStiffnessMatrix(8, 2) = LocalStiffnessMatrix(2, 8);
		LocalStiffnessMatrix(8, 4) = LocalStiffnessMatrix(4, 8);
		LocalStiffnessMatrix(8, 8) = LocalStiffnessMatrix(2, 2);
		LocalStiffnessMatrix(8, 10) = LocalStiffnessMatrix(8, 4);

		LocalStiffnessMatrix(9, 3) = LocalStiffnessMatrix(3, 9);
		LocalStiffnessMatrix(9, 9) = LocalStiffnessMatrix(3, 3);

		LocalStiffnessMatrix(10, 2) = LocalStiffnessMatrix(2, 10);
		LocalStiffnessMatrix(10, 4) = LocalStiffnessMatrix(4, 10);
		LocalStiffnessMatrix(10, 8) = LocalStiffnessMatrix(8, 10);
		LocalStiffnessMatrix(10, 10) = LocalStiffnessMatrix(4, 4);

		LocalStiffnessMatrix(11, 1) = LocalStiffnessMatrix(1, 11);
		LocalStiffnessMatrix(11, 5) = LocalStiffnessMatrix(5, 11);
		LocalStiffnessMatrix(11, 7) = LocalStiffnessMatrix(7, 11);
		LocalStiffnessMatrix(11, 11) = LocalStiffnessMatrix(5, 5);

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}*/

	/*Matrix CrBeamAdjointElement3D2N::CreateElementStiffnessMatrix_Geometry(
		const Vector qe) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension * 2;

		const double N = qe[6];
		const double Mt = qe[9];
		const double my_A = qe[4];
		const double mz_A = qe[5];
		const double my_B = qe[10];
		const double mz_B = qe[11];

		const double L = this->CalculateCurrentLength();
		const double Qy = -1.00 * (mz_A + mz_B) / L;
		const double Qz = (my_A + my_B) / L;

		Matrix LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);

		LocalStiffnessMatrix(0, 1) = -Qy / L;
		LocalStiffnessMatrix(0, 2) = -Qz / L;
		LocalStiffnessMatrix(0, 7) = -1.0 * LocalStiffnessMatrix(0, 1);
		LocalStiffnessMatrix(0, 8) = -1.0 * LocalStiffnessMatrix(0, 2);

		LocalStiffnessMatrix(1, 0) = LocalStiffnessMatrix(0, 1);

		LocalStiffnessMatrix(1, 1) = 1.2 * N / L;

		LocalStiffnessMatrix(1, 3) = my_A / L;
		LocalStiffnessMatrix(1, 4) = Mt / L;

		LocalStiffnessMatrix(1, 5) = N / 10.0;

		LocalStiffnessMatrix(1, 6) = LocalStiffnessMatrix(0, 7);
		LocalStiffnessMatrix(1, 7) = -1.00 * LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(1, 9) = my_B / L;
		LocalStiffnessMatrix(1, 10) = -1.00 * LocalStiffnessMatrix(1, 4);
		LocalStiffnessMatrix(1, 11) = LocalStiffnessMatrix(1, 5);

		LocalStiffnessMatrix(2, 0) = LocalStiffnessMatrix(0, 2);
		LocalStiffnessMatrix(2, 2) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(2, 3) = mz_A / L;
		LocalStiffnessMatrix(2, 4) = -1.00 * LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(2, 5) = LocalStiffnessMatrix(1, 4);
		LocalStiffnessMatrix(2, 6) = LocalStiffnessMatrix(0, 8);
		LocalStiffnessMatrix(2, 8) = LocalStiffnessMatrix(1, 7);
		LocalStiffnessMatrix(2, 9) = mz_B / L;
		LocalStiffnessMatrix(2, 10) = LocalStiffnessMatrix(2, 4);
		LocalStiffnessMatrix(2, 11) = LocalStiffnessMatrix(1, 10);

		for (int i = 0; i < 3; ++i) {
			LocalStiffnessMatrix(3, i) = LocalStiffnessMatrix(i, 3);
		}
		LocalStiffnessMatrix(3, 4) = (-mz_A / 3.00) + (mz_B / 6.00);
		LocalStiffnessMatrix(3, 5) = (my_A / 3.00) - (my_B / 6.00);
		LocalStiffnessMatrix(3, 7) = -my_A / L;
		LocalStiffnessMatrix(3, 8) = -mz_A / L;
		LocalStiffnessMatrix(3, 10) = L*Qy / 6.00;
		LocalStiffnessMatrix(3, 11) = L*Qz / 6.00;

		for (int i = 0; i < 4; ++i) {
			LocalStiffnessMatrix(4, i) = LocalStiffnessMatrix(i, 4);
		}
		LocalStiffnessMatrix(4, 4) = 2.00 * L*N / 15.00;
		LocalStiffnessMatrix(4, 7) = -Mt / L;
		LocalStiffnessMatrix(4, 8) = N / 10.00;
		LocalStiffnessMatrix(4, 9) = LocalStiffnessMatrix(3, 10);
		LocalStiffnessMatrix(4, 10) = -L*N / 30.00;
		LocalStiffnessMatrix(4, 11) = Mt / 2.00;


		for (int i = 0; i < 5; ++i) {
			LocalStiffnessMatrix(5, i) = LocalStiffnessMatrix(i, 5);
		}
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(4, 4);
		LocalStiffnessMatrix(5, 7) = -N / 10.0;
		LocalStiffnessMatrix(5, 8) = -Mt / L;
		LocalStiffnessMatrix(5, 9) = LocalStiffnessMatrix(3, 11);
		LocalStiffnessMatrix(5, 10) = -1.00 * LocalStiffnessMatrix(4, 11);
		LocalStiffnessMatrix(5, 11) = LocalStiffnessMatrix(4, 10);

		for (int i = 0; i < 6; ++i) {
			LocalStiffnessMatrix(6, i) = LocalStiffnessMatrix(i, 6);
		}
		LocalStiffnessMatrix(6, 7) = LocalStiffnessMatrix(0, 1);
		LocalStiffnessMatrix(6, 8) = LocalStiffnessMatrix(0, 2);

		for (int i = 0; i < 7; ++i) {
			LocalStiffnessMatrix(7, i) = LocalStiffnessMatrix(i, 7);
		}
		LocalStiffnessMatrix(7, 7) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(7, 9) = -1.00 * LocalStiffnessMatrix(1, 9);
		LocalStiffnessMatrix(7, 10) = LocalStiffnessMatrix(4, 1);
		LocalStiffnessMatrix(7, 11) = LocalStiffnessMatrix(2, 4);

		for (int i = 0; i < 8; ++i) {
			LocalStiffnessMatrix(8, i) = LocalStiffnessMatrix(i, 8);
		}
		LocalStiffnessMatrix(8, 8) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(8, 9) = -1.00 * LocalStiffnessMatrix(2, 9);
		LocalStiffnessMatrix(8, 10) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(8, 11) = LocalStiffnessMatrix(1, 4);

		for (int i = 0; i < 9; ++i) {
			LocalStiffnessMatrix(9, i) = LocalStiffnessMatrix(i, 9);
		}
		LocalStiffnessMatrix(9, 10) = (mz_A / 6.00) - (mz_B / 3.00);
		LocalStiffnessMatrix(9, 11) = (-my_A / 6.00) + (my_B / 3.00);

		for (int i = 0; i < 10; ++i) {
			LocalStiffnessMatrix(10, i) = LocalStiffnessMatrix(i, 10);
		}
		LocalStiffnessMatrix(10, 10) = LocalStiffnessMatrix(4, 4);

		for (int i = 0; i < 11; ++i) {
			LocalStiffnessMatrix(11, i) = LocalStiffnessMatrix(i, 11);
		}
		LocalStiffnessMatrix(11, 11) = LocalStiffnessMatrix(4, 4);

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}*/

	/*Matrix CrBeamAdjointElement3D2N::CalculateDeformationStiffness() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		Matrix Kd = ZeroMatrix(local_size, local_size);
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = this->CalculateShearModulus();
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();

		Vector inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y) == true) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z) == true) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}
		const double Psi_y = this->CalculatePsi(Iy, Az);
		const double Psi_z = this->CalculatePsi(Iz, Ay);

		Kd(0, 0) = G * J / L;
		Kd(1, 1) = E * Iy / L;
		Kd(2, 2) = E * Iz / L;
		Kd(3, 3) = E * A / L;
		Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
		Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;


		//add geometric stiffness part
		if (this->mIsLinearElement == false)
		{
			const double l = this->CalculateCurrentLength();
			const double N = this->mNodalForces[6];

			const double Qy = -1.00 * (this->mNodalForces[5] +
				this->mNodalForces[11]) / l;

			const double Qz = 1.00 * (this->mNodalForces[4] +
				this->mNodalForces[10]) / l;

			const double N1 = l*N / 12.00;
			const double N2 = l*N / 20.00;
			const double Qy1 = -l*Qy / 6.00;
			const double Qz1 = -l*Qz / 6.00;

			Kd(1, 1) += N1;
			Kd(2, 2) += N1;
			Kd(4, 4) += N2;
			Kd(5, 5) += N2;

			Kd(0, 1) += Qy1;
			Kd(0, 2) += Qz1;
			Kd(1, 0) += Qy1;
			Kd(2, 0) += Qz1;

		}
		return Kd;
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateInitialLocalCS() {

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const unsigned int local_size = size * 2;

		array_1d<double, 3> DirectionVectorX = ZeroVector(dimension);
		array_1d<double, 3> DirectionVectorY = ZeroVector(dimension);
		array_1d<double, 3> DirectionVectorZ = ZeroVector(dimension);
		Vector ReferenceCoordinates = ZeroVector(size);

		ReferenceCoordinates[0] = this->GetGeometry()[0].X0();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X0();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z0();

		for (int i = 0; i < dimension; ++i)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + dimension]
				- ReferenceCoordinates[i]);
		}

		//use orientation class 1st constructor
		double theta_costum = 0.00;
		if (this->GetProperties().Has(ANG_ROT) == true) {
			theta_costum = this->GetProperties()[ANG_ROT];
		}


		Orientation element_axis(DirectionVectorX, theta_costum);
		element_axis.CalculateBasisVectors(DirectionVectorX, DirectionVectorY,
			DirectionVectorZ);
		//save them to update the local axis in every following iter. step
		this->mNX0 = DirectionVectorX;
		this->mNY0 = DirectionVectorY;
		this->mNZ0 = DirectionVectorZ;


		Matrix Temp = ZeroMatrix(dimension);
		this->mRotationMatrix0 = ZeroMatrix(local_size);
		element_axis.CalculateRotationMatrix(Temp);
		this->AssembleSmallInBigMatrix(Temp, this->mRotationMatrix0);
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateTransformationMatrix(Matrix& rRotationMatrix) {

		KRATOS_TRY
			//12x12
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const unsigned int MatSize = 2 * size;

		//initialize local CS
		if (this->mIterationCount == 0) this->CalculateInitialLocalCS();

		//update local CS
		Matrix AuxRotationMatrix = ZeroMatrix(dimension);
		AuxRotationMatrix = this->UpdateRotationMatrixLocal();

		if (rRotationMatrix.size1() != MatSize) {
			rRotationMatrix.resize(MatSize, MatSize, false);
		}

		rRotationMatrix = ZeroMatrix(MatSize);
		//Building the rotation matrix for the local element matrix
		this->AssembleSmallInBigMatrix(AuxRotationMatrix, rRotationMatrix);
		KRATOS_CATCH("")
	}*/

	/*Matrix CrBeamAdjointElement3D2N::CalculateTransformationS() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const unsigned int MatSize = 2 * size;

		const double L = this->CalculateCurrentLength();
		Matrix S = ZeroMatrix(MatSize, size);
		S(0, 3) = -1.00;
		S(1, 5) = 2.00 / L;
		S(2, 4) = -2.00 / L;
		S(3, 0) = -1.00;
		S(4, 1) = -1.00;
		S(4, 4) = 1.00;
		S(5, 2) = -1.00;
		S(5, 5) = 1.00;
		S(6, 3) = 1.00;
		S(7, 5) = -2.00 / L;
		S(8, 4) = 2.00 / L;
		S(9, 0) = 1.00;
		S(10, 1) = 1.00;
		S(10, 4) = 1.00;
		S(11, 2) = 1.00;
		S(11, 5) = 1.00;

		return S;
		KRATOS_CATCH("")
	}*/

	/*Matrix CrBeamAdjointElement3D2N::UpdateRotationMatrixLocal() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const unsigned int MatSize = 2 * size;

		Vector dPhiA = ZeroVector(dimension);
		Vector dPhiB = ZeroVector(dimension);
		Vector IncrementDeformation = ZeroVector(MatSize);
		IncrementDeformation = this->mIncrementDeformation;

		for (int i = 0; i < dimension; ++i) {
			dPhiA[i] = IncrementDeformation[i + 3];
			dPhiB[i] = IncrementDeformation[i + 9];
		}

		//calculating quaternions
		Vector drA_vec = ZeroVector(dimension);
		Vector drB_vec = ZeroVector(dimension);
		double drA_sca, drB_sca;

		drA_vec = 0.50 * dPhiA;
		drB_vec = 0.50 * dPhiB;

		drA_sca = 0.00;
		drB_sca = 0.00;
		for (int i = 0; i < dimension; ++i) {
			drA_sca += drA_vec[i] * drA_vec[i];
			drB_sca += drB_vec[i] * drB_vec[i];
		}
		drA_sca = 1.00 - drA_sca;
		drB_sca = 1.00 - drB_sca;

		drA_sca = sqrt(drA_sca);
		drB_sca = sqrt(drB_sca);


		//1st solution step
		if (mIterationCount == 0) {
			this->mQuaternionVEC_A = ZeroVector(dimension);
			this->mQuaternionVEC_B = ZeroVector(dimension);
			this->mQuaternionSCA_A = 1.00;
			this->mQuaternionSCA_B = 1.00;
		}

		Vector tempVec = ZeroVector(dimension);
		double tempSca = 0.00;

		//Node A
		tempVec = this->mQuaternionVEC_A;
		tempSca = this->mQuaternionSCA_A;

		this->mQuaternionSCA_A = drA_sca *tempSca;
		for (int i = 0; i < dimension; ++i) {
			this->mQuaternionSCA_A -= drA_vec[i] * tempVec[i];
		}
		this->mQuaternionVEC_A = drA_sca*tempVec;
		this->mQuaternionVEC_A += tempSca * drA_vec;
		this->mQuaternionVEC_A += MathUtils<double>::CrossProduct(drA_vec, tempVec);

		//Node B
		tempVec = this->mQuaternionVEC_B;
		tempSca = this->mQuaternionSCA_B;

		this->mQuaternionSCA_B = drB_sca *tempSca;
		for (int i = 0; i < dimension; ++i) {
			this->mQuaternionSCA_B -= drB_vec[i] * tempVec[i];
		}

		this->mQuaternionVEC_B = drB_sca*tempVec;
		this->mQuaternionVEC_B += tempSca * drB_vec;
		this->mQuaternionVEC_B += MathUtils<double>::CrossProduct(drB_vec, tempVec);


		//scalar part of difference quaternion
		double scalar_diff;
		scalar_diff = (this->mQuaternionSCA_A + this->mQuaternionSCA_B) *
			(this->mQuaternionSCA_A + this->mQuaternionSCA_B);

		tempVec = this->mQuaternionVEC_A + this->mQuaternionVEC_B;
		scalar_diff += MathUtils<double>::Norm(tempVec) *
			MathUtils<double>::Norm(tempVec);

		scalar_diff = 0.50 * sqrt(scalar_diff);

		//mean rotation quaternion
		double meanRotationScalar;
		meanRotationScalar = (this->mQuaternionSCA_A + this->mQuaternionSCA_B) * 0.50;
		meanRotationScalar = meanRotationScalar / scalar_diff;

		Vector meanRotationVector = ZeroVector(dimension);
		meanRotationVector = (this->mQuaternionVEC_A + this->mQuaternionVEC_B) * 0.50;
		meanRotationVector = meanRotationVector / scalar_diff;

		//vector part of difference quaternion
		Vector vector_diff = ZeroVector(dimension);
		vector_diff = this->mQuaternionSCA_A * this->mQuaternionVEC_B;
		vector_diff -= this->mQuaternionSCA_B * this->mQuaternionVEC_A;
		vector_diff += MathUtils<double>::CrossProduct(this->mQuaternionVEC_A,
			this->mQuaternionVEC_B);

		vector_diff = 0.50 * vector_diff / scalar_diff;

		//rotate inital element basis
		const double r0 = meanRotationScalar;
		const double r1 = meanRotationVector[0];
		const double r2 = meanRotationVector[1];
		const double r3 = meanRotationVector[2];

		Quaternion<double> q(r0, r1, r2, r3);
		Vector rotatedNX0 = this->mNX0;
		Vector rotatedNY0 = this->mNY0;
		Vector rotatedNZ0 = this->mNZ0;
		q.RotateVector3(rotatedNX0);
		q.RotateVector3(rotatedNY0);
		q.RotateVector3(rotatedNZ0);

		Matrix RotatedCS = ZeroMatrix(dimension, dimension);
		for (int i = 0; i < dimension; ++i) {
			RotatedCS(i, 0) = rotatedNX0[i];
			RotatedCS(i, 1) = rotatedNY0[i];
			RotatedCS(i, 2) = rotatedNZ0[i];
		}

		//rotate basis to element axis + redefine R
		Vector n_bisectrix = ZeroVector(dimension);
		Vector deltaX = ZeroVector(dimension);
		double VectorNorm;

		deltaX[0] = this->mTotalNodalPosistion[3] - this->mTotalNodalPosistion[0];
		deltaX[1] = this->mTotalNodalPosistion[4] - this->mTotalNodalPosistion[1];
		deltaX[2] = this->mTotalNodalPosistion[5] - this->mTotalNodalPosistion[2];


		VectorNorm = MathUtils<double>::Norm(deltaX);
		if (VectorNorm != 0.00) deltaX /= VectorNorm;


		n_bisectrix = rotatedNX0 + deltaX;
		VectorNorm = MathUtils<double>::Norm(n_bisectrix);
		if (VectorNorm != 0.00) n_bisectrix /= VectorNorm;

		Matrix n_xyz = ZeroMatrix(dimension);
		for (int i = 0; i < dimension; ++i) {
			n_xyz(i, 0) = -1.0 * RotatedCS(i, 0);
			n_xyz(i, 1) = 1.0 * RotatedCS(i, 1);
			n_xyz(i, 2) = 1.0 * RotatedCS(i, 2);
		}

		Matrix Identity = ZeroMatrix(dimension);
		for (int i = 0; i < dimension; ++i) Identity(i, i) = 1.0;
		Identity -= 2.0 * outer_prod(n_bisectrix, n_bisectrix);
		n_xyz = prod(Identity, n_xyz);


		//save current CS for GID OUTPUT
		this->mNX = ZeroVector(dimension);
		this->mNY = ZeroVector(dimension);
		this->mNZ = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i)
		{
			this->mNX[i] = n_xyz(i, 0);
			this->mNY[i] = n_xyz(i, 1);
			this->mNZ[i] = n_xyz(i, 2);
		}

		//calculating deformation modes
		this->mPhiS = ZeroVector(dimension);
		this->mPhiA = ZeroVector(dimension);
		this->mPhiS = prod(Matrix(trans(n_xyz)), vector_diff);
		this->mPhiS *= 4.00;

		rotatedNX0 = ZeroVector(dimension);
		tempVec = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i) rotatedNX0[i] = n_xyz(i, 0);
		tempVec = MathUtils<double>::CrossProduct(rotatedNX0, n_bisectrix);
		this->mPhiA = prod(Matrix(trans(n_xyz)), tempVec);
		this->mPhiA *= 4.00;

		if (this->mIterationCount == 0)
		{
			this->mPhiS = ZeroVector(dimension);
			this->mPhiA = ZeroVector(dimension);
		}
		return n_xyz;
		KRATOS_CATCH("")
	}*/

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

	/*void CrBeamAdjointElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension * 2;
			rValues[index] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Z, Step);

			rValues[index + 3] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_X, Step);
			rValues[index + 4] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_Y, Step);
			rValues[index + 5] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_Z, Step);
		}

		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension * 2;

			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Z, Step);

			rValues[index + 3] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_X, Step);
			rValues[index + 4] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_Y, Step);
			rValues[index + 5] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z, Step);
		}
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);



		if (this->GetProperties().Has(LUMPED_MASS_MATRIX) == true) {
			this->mIsLumpedMassMatrix = GetProperties()[LUMPED_MASS_MATRIX];
		}
		else this->mIsLumpedMassMatrix = false;



		if (this->mIsLumpedMassMatrix == true)
		{
			this->CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
		}
		else
		{
			this->CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);

			Matrix RotationMatrix = ZeroMatrix(MatSize);
			Matrix aux_matrix = ZeroMatrix(MatSize);

			RotationMatrix = this->mRotationMatrix;
			aux_matrix = prod(RotationMatrix, rMassMatrix);
			rMassMatrix = prod(aux_matrix,
				Matrix(trans(RotationMatrix)));
		}
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension * 2;

		if (rDampingMatrix.size1() != MatSize)
		{
			rDampingMatrix.resize(MatSize, MatSize, false);
		}

		rDampingMatrix = ZeroMatrix(MatSize, MatSize);

		Matrix StiffnessMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		Matrix MassMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

		double alpha = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_ALPHA))
		{
			alpha = this->GetProperties()[RAYLEIGH_ALPHA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
		{
			alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
		}

		double beta = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_BETA))
		{
			beta = this->GetProperties()[RAYLEIGH_BETA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
		{
			beta = rCurrentProcessInfo[RAYLEIGH_BETA];
		}

		rDampingMatrix += alpha * MassMatrix;
		rDampingMatrix += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}*/


	/*Vector CrBeamAdjointElement3D2N::CalculateBodyForces()
	{
		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int localSize = number_of_nodes * dimension;
		const unsigned int MatSize = number_of_nodes * dimension * 2;

		//getting shapefunctionvalues for linear SF
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		Vector EquivalentLineLoad = ZeroVector(dimension);
		Vector BodyForcesGlobal = ZeroVector(MatSize);

		const double A = this->GetProperties()[CROSS_AREA];
		const double l = this->CalculateCurrentLength();
		const double rho = this->GetProperties()[DENSITY];

		//calculating equivalent line load
		for (int i = 0; i < number_of_nodes; ++i)
		{
			EquivalentLineLoad += A * rho*
				this->GetGeometry()[i].
				FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0, i);
		}


		// adding the nodal forces
		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i*localSize;
			for (int j = 0; j < dimension; ++j)
			{
				BodyForcesGlobal[j + index] =
					EquivalentLineLoad[j] * Ncontainer(0, i) * l;
			}
		}

		// adding the nodal moments
		this->CalculateAndAddWorkEquivalentNodalForcesLineLoad
			(EquivalentLineLoad, BodyForcesGlobal, l);


		// return the total ForceVector
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
		const Vector ForceInput, VectorType& rRightHandSideVector,
		const double GeometryLength)
	{
		KRATOS_TRY;
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		//calculate orthogonal load vector
		Vector GeometricOrientation = ZeroVector(dimension);
		GeometricOrientation[0] = this->GetGeometry()[1].X()
			- this->GetGeometry()[0].X();
		GeometricOrientation[1] = this->GetGeometry()[1].Y()
			- this->GetGeometry()[0].Y();
		if (dimension == 3)
		{
			GeometricOrientation[2] = this->GetGeometry()[1].Z()
				- this->GetGeometry()[0].Z();
		}

		const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
		if (VectorNormA != 0.00) GeometricOrientation /= VectorNormA;

		Vector LineLoadDir = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i)
		{
			LineLoadDir[i] = ForceInput[i];
		}

		const double VectorNormB = MathUtils<double>::Norm(LineLoadDir);
		if (VectorNormB != 0.00) LineLoadDir /= VectorNormB;

		double cosAngle = 0.00;
		for (int i = 0; i < dimension; ++i)
		{
			cosAngle += LineLoadDir[i] * GeometricOrientation[i];
		}

		const double sinAngle = sqrt(1.00 - (cosAngle*cosAngle));
		const double NormForceVectorOrth = sinAngle * VectorNormB;


		Vector NodeA = ZeroVector(dimension);
		NodeA[0] = this->GetGeometry()[0].X();
		NodeA[1] = this->GetGeometry()[0].Y();
		if (dimension == 3)	NodeA[2] = this->GetGeometry()[0].Z();

		Vector NodeB = ZeroVector(dimension);
		NodeB = NodeA + LineLoadDir;

		Vector NodeC = ZeroVector(dimension);
		NodeC = NodeA + (GeometricOrientation*cosAngle);

		Vector LoadOrthogonalDir = ZeroVector(dimension);
		LoadOrthogonalDir = NodeB - NodeC;
		const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
		if (VectorNormC != 0.00) LoadOrthogonalDir /= VectorNormC;



		// now caluclate respective work equivilent nodal moments

		const double CustomMoment = NormForceVectorOrth *
			GeometryLength*GeometryLength / 12.00;

		Vector MomentNodeA = ZeroVector(dimension);
		MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
			LoadOrthogonalDir);
		MomentNodeA *= CustomMoment;

		for (int i = 0; i < dimension; ++i)
		{
			rRightHandSideVector[(1 * dimension) + i] += MomentNodeA[i];
			rRightHandSideVector[(3 * dimension) + i] -= MomentNodeA[i];
		}

		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension * 2;

		this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

		//Nodal element forces global
		Vector nodalForcesGlobal_q = ZeroVector(LocalSize);
		nodalForcesGlobal_q = prod(this->mRotationMatrix,
			this->mNodalForces);

		//create+compute RHS
		//update Residual
		rRightHandSideVector = ZeroVector(LocalSize);
		rRightHandSideVector -= nodalForcesGlobal_q;


		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
		}
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		this->mIterationCount++;
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY;
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;
		rRightHandSideVector = ZeroVector(LocalSize);


		if (this->mIsLinearElement == false)
		{
			this->UpdateIncrementDeformation();
			Matrix TransformationMatrix = ZeroMatrix(LocalSize);
			this->CalculateTransformationMatrix(TransformationMatrix);
			Vector elementForces_t = ZeroVector(size);
			elementForces_t = this->CalculateElementForces();
			Vector nodalForcesLocal_qe = ZeroVector(LocalSize);
			Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
			TransformationMatrixS = this->CalculateTransformationS();
			nodalForcesLocal_qe = prod(TransformationMatrixS,
				elementForces_t);
			//save local nodal forces
			this->mNodalForces = ZeroVector(LocalSize);
			this->mNodalForces = nodalForcesLocal_qe;

			Vector nodalForcesGlobal_q = ZeroVector(LocalSize);
			nodalForcesGlobal_q = prod(TransformationMatrix, nodalForcesLocal_qe);
			rRightHandSideVector -= nodalForcesGlobal_q;
		}

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
		}

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		KRATOS_CATCH("")

	}*/

	/*void CrBeamAdjointElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		//update displacement_delta
		this->UpdateIncrementDeformation();

		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		this->mRotationMatrix = ZeroMatrix(LocalSize);
		this->mRotationMatrix = TransformationMatrix;

		//deformation modes
		Vector elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();

		//Nodal element forces local
		Vector nodalForcesLocal_qe = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		nodalForcesLocal_qe = prod(TransformationMatrixS, elementForces_t);

		//save local nodal forces
		this->mNodalForces = ZeroVector(LocalSize);
		this->mNodalForces = nodalForcesLocal_qe;

		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Material();
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Geometry(nodalForcesLocal_qe);


		Matrix aux_matrix = ZeroMatrix(LocalSize);
		aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
		rLeftHandSideMatrix = prod(aux_matrix,
			Matrix(trans(TransformationMatrix)));

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			TransformationMatrix = this->mRotationMatrix0;
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			rLeftHandSideMatrix +=
				this->CreateElementStiffnessMatrix_Material();
			aux_matrix = ZeroMatrix(LocalSize);
			aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
			rLeftHandSideMatrix = prod(aux_matrix,
				Matrix(trans(TransformationMatrix)));
		}
		//assign global element variables
		this->mLHS = rLeftHandSideMatrix;
		KRATOS_CATCH("")
	}*/

	/*Vector CrBeamAdjointElement3D2N::CalculateElementForces() {

		KRATOS_TRY;
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		Vector deformation_modes_total_V = ZeroVector(LocalSize);
		const double L = this->CalculateReferenceLength();
		const double l = this->CalculateCurrentLength();

		deformation_modes_total_V[3] = l - L;
		for (int i = 0; i < 3; ++i) deformation_modes_total_V[i] = this->mPhiS[i];
		for (int i = 0; i < 2; ++i) deformation_modes_total_V[i + 4] = this->mPhiA[i + 1];
		//calculate element forces
		Vector element_forces_t = ZeroVector(LocalSize);
		Matrix deformation_stiffness_Kd = ZeroMatrix(LocalSize);

		deformation_stiffness_Kd = this->CalculateDeformationStiffness();
		element_forces_t = prod(deformation_stiffness_Kd,
			deformation_modes_total_V);

		return element_forces_t;
		KRATOS_CATCH("")
	}*/

	/*double CrBeamAdjointElement3D2N::CalculateCurrentLength() {

		KRATOS_TRY;
		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) +
			(dw + dz)*(dw + dz));
		return l;
		KRATOS_CATCH("")

	}*/

	/*double CrBeamAdjointElement3D2N::CalculatePsi(const double I, const double A_eff) {

		KRATOS_TRY;
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double L = this->CalculateCurrentLength();
		const double G = this->CalculateShearModulus();

		const double phi = (12.0 * E * I) / (L*L * G*A_eff);
		double psi;
		//interpret input A_eff == 0 as shearstiff -> psi = 1.0
		if (A_eff == 0.00) psi = 1.00;
		else psi = 1.0 / (1.0 + phi);

		return psi;
		KRATOS_CATCH("")
	}*/

	/*double CrBeamAdjointElement3D2N::CalculateReferenceLength() {

		KRATOS_TRY;
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::UpdateIncrementDeformation() {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		Vector actualDeformation = ZeroVector(LocalSize);
		Vector total_nodal_def = ZeroVector(LocalSize);
		Vector total_nodal_pos = ZeroVector(size);
		this->mIncrementDeformation = ZeroVector(LocalSize);

		if (mIterationCount == 0) this->mTotalNodalDeformation = ZeroVector(LocalSize);
		this->GetValuesVector(actualDeformation, 0);

		this->mIncrementDeformation = actualDeformation
			- this->mTotalNodalDeformation;

		this->mTotalNodalDeformation = ZeroVector(LocalSize);
		this->mTotalNodalDeformation = actualDeformation;

		this->mTotalNodalPosistion = ZeroVector(size);
		this->mTotalNodalPosistion[0] = this->GetGeometry()[0].X0()
			+ actualDeformation[0];
		this->mTotalNodalPosistion[1] = this->GetGeometry()[0].Y0()
			+ actualDeformation[1];
		this->mTotalNodalPosistion[2] = this->GetGeometry()[0].Z0()
			+ actualDeformation[2];

		this->mTotalNodalPosistion[3] = this->GetGeometry()[1].X0()
			+ actualDeformation[6];
		this->mTotalNodalPosistion[4] = this->GetGeometry()[1].Y0()
			+ actualDeformation[7];
		this->mTotalNodalPosistion[5] = this->GetGeometry()[1].Z0()
			+ actualDeformation[8];
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;


		//element with two nodes can only represent results at one node 
		const unsigned int&  write_points_number = GetGeometry()
			.IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
		if (rOutput.size() != write_points_number) {
			rOutput.resize(write_points_number);
		}


		this->UpdateIncrementDeformation();
		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		//deformation modes
		Vector elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();
		Vector Stress = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		Stress = prod(TransformationMatrixS, elementForces_t);

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			LeftHandSideMatrix = this->mLHS;

			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			Stress = ZeroVector(LocalSize);
			Stress = prod(LeftHandSideMatrix, NodalDeformation);
			Matrix TransformationMatrix = ZeroMatrix(LocalSize);
			TransformationMatrix = this->mRotationMatrix;
			Stress = prod(Matrix(trans(TransformationMatrix)), Stress);
		}


		//rOutput[GP 1,2,3][x,y,z]

		if (rVariable == MOMENT)
		{
			rOutput[0][0] = -1.0 *Stress[3] * 0.75 + Stress[9] * 0.25;
			rOutput[1][0] = -1.0 *Stress[3] * 0.50 + Stress[9] * 0.50;
			rOutput[2][0] = -1.0 *Stress[3] * 0.25 + Stress[9] * 0.75;

			rOutput[0][1] = -1.0 *Stress[4] * 0.75 + Stress[10] * 0.25;
			rOutput[1][1] = -1.0 *Stress[4] * 0.50 + Stress[10] * 0.50;
			rOutput[2][1] = -1.0 *Stress[4] * 0.25 + Stress[10] * 0.75;

			rOutput[0][2] = -1.0 *Stress[5] * 0.75 + Stress[11] * 0.25;
			rOutput[1][2] = -1.0 *Stress[5] * 0.50 + Stress[11] * 0.50;
			rOutput[2][2] = -1.0 *Stress[5] * 0.25 + Stress[11] * 0.75;

		}
		if (rVariable == FORCE)
		{
			rOutput[0][0] = -1.0 * Stress[0] * 0.75 + Stress[6] * 0.25;
			rOutput[1][0] = -1.0 * Stress[0] * 0.50 + Stress[6] * 0.50;
			rOutput[2][0] = -1.0 * Stress[0] * 0.25 + Stress[6] * 0.75;

			rOutput[0][1] = -1.0 * Stress[1] * 0.75 + Stress[7] * 0.25;
			rOutput[1][1] = -1.0 *Stress[1] * 0.50 + Stress[7] * 0.50;
			rOutput[2][1] = -1.0 *Stress[1] * 0.25 + Stress[7] * 0.75;

			rOutput[0][2] = -1.0 *Stress[2] * 0.75 + Stress[8] * 0.25;
			rOutput[1][2] = -1.0 *Stress[2] * 0.50 + Stress[8] * 0.50;
			rOutput[2][2] = -1.0 *Stress[2] * 0.25 + Stress[8] * 0.75;

		}

		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::GetValueOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}*/



	/*void CrBeamAdjointElement3D2N::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
		std::vector< Vector >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		if (rVariable == LOCAL_AXES_VECTOR)
		{
			rOutput.resize(3);
			for (int i = 0; i < 3; ++i) rOutput[i] = ZeroVector(3);

			if (this->mIsLinearElement == true)
			{
				rOutput[0] = this->mNX0;
				rOutput[1] = this->mNY0;
				rOutput[2] = this->mNZ0;
			}
			else
			{
				rOutput[0] = this->mNX;
				rOutput[1] = this->mNY;
				rOutput[2] = this->mNZ;
			}
		}

		KRATOS_CATCH("");
	}*/

	/*void CrBeamAdjointElement3D2N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::AssembleSmallInBigMatrix(Matrix SmallMatrix,
		Matrix& BigMatrix) {

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int size = number_of_nodes * dimension;
		const unsigned int MatSize = 2 * size;

		for (unsigned int kk = 0; kk < MatSize; kk += dimension)
		{
			for (int i = 0; i<dimension; ++i)
			{
				for (int j = 0; j<dimension; ++j)
				{
					BigMatrix(i + kk, j + kk) = SmallMatrix(i, j);
				}
			}
		}
		KRATOS_CATCH("")
	}*/


	/*void CrBeamAdjointElement3D2N::BuildSingleMassMatrix(MatrixType& rMassMatrix,
		double Phi, double CT, double CR, double L)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int MatSize = number_of_nodes * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);
		Matrix TempMassMatrix = ZeroMatrix(MatSize, MatSize);
		const double Phi2 = Phi * Phi;
		const double L2 = L*L;


		TempMassMatrix(0, 0) = (13.00 / 35.00) + (7.00 / 10.00)*Phi
			+ (1.00 / 3.00)*Phi2;
		TempMassMatrix(0, 1) = ((11.00 / 210.00) + (11.00 / 210.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(0, 2) = (9.00 / 70.00) + (3.00 / 10.00)*Phi
			+ (1.00 / 6.00)*Phi2;
		TempMassMatrix(0, 3) = -((13.00 / 420.00) + (3.00 / 40.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(1, 0) = TempMassMatrix(0, 1);
		TempMassMatrix(1, 1) = ((1.00 / 105.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;
		TempMassMatrix(1, 2) = ((13.00 / 420.00) + (3.00 / 40.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(1, 3) = -((1.00 / 140.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;
		TempMassMatrix(2, 0) = TempMassMatrix(0, 2);
		TempMassMatrix(2, 1) = TempMassMatrix(1, 2);
		TempMassMatrix(2, 2) = (13.00 / 35.00) + (7.00 / 10.00)*Phi
			+ (1.00 / 3.00)*Phi2;
		TempMassMatrix(2, 3) = -((11.00 / 210.00) + (11.00 / 210.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(3, 0) = TempMassMatrix(0, 3);
		TempMassMatrix(3, 1) = TempMassMatrix(1, 3);
		TempMassMatrix(3, 2) = TempMassMatrix(2, 3);
		TempMassMatrix(3, 3) = ((1.00 / 105.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;

		TempMassMatrix *= CT;
		rMassMatrix += TempMassMatrix;


		TempMassMatrix = ZeroMatrix(MatSize, MatSize);

		TempMassMatrix(0, 0) = 6.00 / 5.00;
		TempMassMatrix(0, 1) = ((1.00 / 10.00) - (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(0, 2) = -6.00 / 5.00;
		TempMassMatrix(0, 3) = ((1.00 / 10.00) - (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(1, 0) = TempMassMatrix(0, 1);
		TempMassMatrix(1, 1) = ((2.00 / 15.00) + (1.00 / 6.00)*Phi
			+ (1.00 / 3.00)*Phi2)*L2;
		TempMassMatrix(1, 2) = ((-1.00 / 10.00) + (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(1, 3) = -((1.00 / 30.00) + (1.00 / 6.00)*Phi
			- (1.00 / 6.00)*Phi2)*L2;
		TempMassMatrix(2, 0) = TempMassMatrix(0, 2);
		TempMassMatrix(2, 1) = TempMassMatrix(1, 2);
		TempMassMatrix(2, 2) = 6.00 / 5.00;
		TempMassMatrix(2, 3) = ((-1.00 / 10.00) + (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(3, 0) = TempMassMatrix(0, 3);
		TempMassMatrix(3, 1) = TempMassMatrix(1, 3);
		TempMassMatrix(3, 2) = TempMassMatrix(2, 3);
		TempMassMatrix(3, 3) = ((2.00 / 15.00) + (1.00 / 6.00)*Phi
			+ (1.00 / 3.00)*Phi2)*L2;

		TempMassMatrix *= CR;
		rMassMatrix += TempMassMatrix;
		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateConsistentMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const int smallMatSize = number_of_nodes * 2;
		const unsigned int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double L = this->CalculateReferenceLength();
		const double L2 = L * L;
		const double rho = this->GetProperties()[DENSITY];
		const double A = this->GetProperties()[CROSS_AREA];
		const double E = this->GetProperties()[YOUNG_MODULUS];
		Vector inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];
		const double G = this->CalculateShearModulus();

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y) == true) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z) == true) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}

		double IRy = Iy;
		if (this->GetProperties().Has(INERTIA_ROT_Y) == true) {
			IRy = GetProperties()[INERTIA_ROT_Y];
		}

		double IRz = Iz;
		if (this->GetProperties().Has(INERTIA_ROT_Y) == true) {
			IRz = GetProperties()[INERTIA_ROT_Z];
		}

		double Phiy = 0.00;
		double Phiz = 0.00;

		if (Ay != 0.00) Phiz = (12.00 * E * Iz) / (L2*G*Ay);
		if (Az != 0.00) Phiy = (12.00 * E * Iy) / (L2*G*Az);

		const double CTy = (rho * A * L) / ((1 + Phiy)*(1 + Phiy));
		const double CTz = (rho * A * L) / ((1 + Phiz)*(1 + Phiz));

		const double CRy = (rho*IRy) / ((1 + Phiy)*(1 + Phiy)*L);
		const double CRz = (rho*IRz) / ((1 + Phiz)*(1 + Phiz)*L);

		//longitudinal forces + torsional moment
		const double M00 = (1.00 / 3.00)*A*rho*L;
		const double M06 = M00 / 2.00;
		const double M33 = (J*L*rho) / 3.00;
		const double M39 = M33 / 2.00;

		rMassMatrix(0, 0) = M00;
		rMassMatrix(0, 6) = M06;
		rMassMatrix(6, 6) = M00;
		rMassMatrix(3, 3) = M33;
		rMassMatrix(3, 9) = M39;
		rMassMatrix(9, 9) = M33;

		Matrix TempBendingMassMatrix = ZeroMatrix(smallMatSize, smallMatSize);
		this->BuildSingleMassMatrix(TempBendingMassMatrix, Phiz, CTz, CRz, L);

		rMassMatrix(1, 1) = TempBendingMassMatrix(0, 0);
		rMassMatrix(1, 5) = TempBendingMassMatrix(0, 1);
		rMassMatrix(1, 7) = TempBendingMassMatrix(0, 2);
		rMassMatrix(1, 11) = TempBendingMassMatrix(0, 3);
		rMassMatrix(5, 5) = TempBendingMassMatrix(1, 1);
		rMassMatrix(5, 7) = TempBendingMassMatrix(1, 2);
		rMassMatrix(5, 11) = TempBendingMassMatrix(1, 3);
		rMassMatrix(7, 7) = TempBendingMassMatrix(2, 2);
		rMassMatrix(7, 11) = TempBendingMassMatrix(2, 3);
		rMassMatrix(11, 11) = TempBendingMassMatrix(3, 3);

		TempBendingMassMatrix = ZeroMatrix(smallMatSize, smallMatSize);
		this->BuildSingleMassMatrix(TempBendingMassMatrix, Phiy, CTy, CRy, L);

		rMassMatrix(2, 2) = TempBendingMassMatrix(0, 0);
		rMassMatrix(2, 4) = TempBendingMassMatrix(0, 1);
		rMassMatrix(2, 8) = TempBendingMassMatrix(0, 2);
		rMassMatrix(2, 10) = TempBendingMassMatrix(0, 3);
		rMassMatrix(4, 4) = TempBendingMassMatrix(1, 1);
		rMassMatrix(4, 8) = TempBendingMassMatrix(1, 2);
		rMassMatrix(4, 10) = TempBendingMassMatrix(1, 3);
		rMassMatrix(8, 8) = TempBendingMassMatrix(2, 2);
		rMassMatrix(8, 10) = TempBendingMassMatrix(2, 3);
		rMassMatrix(10, 10) = TempBendingMassMatrix(3, 3);


		for (int j = 1; j < 12; ++j)
		{
			for (int i = 0; i < j; ++i)
			{
				rMassMatrix(j, i) = rMassMatrix(i, j);
			}
		}

		KRATOS_CATCH("")
	}*/

	/*void CrBeamAdjointElement3D2N::CalculateLumpedMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		const double TotalMass = A * L * rho;
		const double temp = 0.50 * TotalMass;

		//translatonal mass	
		for (int i = 0; i < number_of_nodes; ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				int index = i * (dimension * 2) + j;
				rMassMatrix(index, index) = temp;
			}
		}
		//rotaional mass neglected alpha = 0
		KRATOS_CATCH("")
	}*/

	/*CrBeamAdjointElement3D2N::IntegrationMethod
		CrBeamAdjointElement3D2N::GetIntegrationMethod() const
	{
		//do this to have 3GP as an output in GID
		return Kratos::GeometryData::GI_GAUSS_3;
	}*/

	/*void CrBeamAdjointElement3D2N::AddExplicitContribution(const VectorType& rRHSVector,
		const Variable<VectorType>& rRHSVariable,
		Variable<array_1d<double, 3> >& rDestinationVariable,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension;

		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
		{

			for (int i = 0; i< number_of_nodes; ++i)
			{
				int index = element_size * i;

				GetGeometry()[i].SetLock();

				array_1d<double, 3> &ForceResidual =
					GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

				for (int j = 0; j<dimension; ++j)
				{
					ForceResidual[j] += rRHSVector[index + j];
				}

				GetGeometry()[i].UnSetLock();
			}
		}


		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL)
		{

			for (int i = 0; i< number_of_nodes; ++i)
			{
				int index = (element_size * i) + dimension;

				GetGeometry()[i].SetLock();

				array_1d<double, 3> &MomentResidual =
					GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

				for (int j = 0; j<dimension; ++j)
				{
					MomentResidual[j] += rRHSVector[index + j];
				}

				GetGeometry()[i].UnSetLock();
			}
		}

		KRATOS_CATCH("")
	}*/

	/*double CrBeamAdjointElement3D2N::CalculateShearModulus() {
		KRATOS_TRY;
		const double nu = this->GetProperties()[POISSON_RATIO];
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = E / (2.0 * (1.0 + nu));
		return G;
		KRATOS_CATCH("")
	}*/

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

		if (this->GetProperties().Has(LOCAL_INERTIA_VECTOR) == false)
		{
			KRATOS_ERROR << "LOCAL_INERTIA_VECTOR not provided for this element" << this->Id()
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



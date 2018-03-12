// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

#include "shell_thin_adjoint_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "structural_mechanics_application_variables.h"

//#include "geometries/triangle_3d_3.h"

#include <string>
#include <iomanip>

#define OPT_NUM_NODES 3 //---------------------------------------------->TODO: check which of the stuff is here really needed
#define OPT_STRAIN_SIZE 6
#define OPT_NUM_DOFS 18

//----------------------------------------
// preprocessors for the integration
// method used by this element.

//#define OPT_1_POINT_INTEGRATION

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_1
#define OPT_NUM_GP 1
#else
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_2
#define OPT_NUM_GP 3
#endif // OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) Utilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

//#define OPT_AVARAGE_RESULTS

namespace Kratos
{

namespace Utilities
{
inline void InterpToStandardGaussPoints(double& v1, double& v2, double& v3)
{
    double vg1 = v1;
    double vg2 = v2;
    double vg3 = v3;
#ifdef OPT_AVARAGE_RESULTS
    v1 = (vg1+vg2+vg3)/3.0;
    v2 = (vg1+vg2+vg3)/3.0;
    v3 = (vg1+vg2+vg3)/3.0;
#else
    v1 = (2.0*vg1)/3.0 - vg2/3.0       + (2.0*vg3)/3.0;
    v2 = (2.0*vg1)/3.0 + (2.0*vg2)/3.0 - vg3/3.0;
    v3 = (2.0*vg2)/3.0 - vg1/3.0       + (2.0*vg3)/3.0;
#endif // OPT_AVARAGE_RESULTS
}

inline void InterpToStandardGaussPoints(std::vector< double >& v)
{
    if(v.size() != 3) return;
    InterpToStandardGaussPoints(v[0], v[1], v[2]);
}

inline void InterpToStandardGaussPoints(std::vector< array_1d<double,3> >& v)
{
    if(v.size() != 3) return;
    for(size_t i = 0; i < 3; i++)
        InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< array_1d<double,6> >& v)
{
    if(v.size() != 3) return;
    for(size_t i = 0; i < 6; i++)
        InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< Vector >& v)
{
    if(v.size() != 3) return;
    size_t ncomp = v[0].size();
    for(int i = 1; i < 3; i++)
        if(v[i].size() != ncomp)
            return;
    for(size_t i = 0; i < ncomp; i++)
        InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< Matrix >& v)
{
    if(v.size() != 3) return;
    size_t nrows = v[0].size1();
    size_t ncols = v[0].size2();
    for(int i = 1; i < 3; i++)
        if(v[i].size1() != nrows || v[i].size2() != ncols)
            return;
    for(size_t i = 0; i < nrows; i++)
        for(size_t j = 0; j < ncols; j++)
            InterpToStandardGaussPoints(v[0](i,j), v[1](i,j), v[2](i,j));
}

}


// =====================================================================================
//
// Class ShellThinAdjointElement3D3N
//
// =====================================================================================

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        bool NLGeom)
    : ShellThinElement3D3N(NewId, pGeometry, NLGeom)
{

}

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        bool NLGeom)
    : ShellThinElement3D3N(NewId, pGeometry, pProperties, NLGeom)
{

}

ShellThinAdjointElement3D3N::ShellThinAdjointElement3D3N(IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        CoordinateTransformationBasePointerType pCoordinateTransformation)
    : ShellThinElement3D3N(NewId, pGeometry, pProperties, pCoordinateTransformation)
{

}

ShellThinAdjointElement3D3N::~ShellThinAdjointElement3D3N()
{
}

Element::Pointer ShellThinAdjointElement3D3N::Create(IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const 
{
    KRATOS_TRY
    bool NLGeom = false; //-----------------------------------------------------------> hard coded linar shell element!!!
    return Element::Pointer(
                new ShellThinAdjointElement3D3N(NewId, pGeom, pProperties, NLGeom));
    
    KRATOS_CATCH("")
}

Element::Pointer ShellThinAdjointElement3D3N::Create(IndexType NewId, 
                    NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
{
    bool NLGeom = false; //-----------------------------------------------------------> hard coded linar shell element!!!
    GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
    return boost::make_shared< ShellThinAdjointElement3D3N >(NewId, newGeom, pProperties, NLGeom);

}

void ShellThinAdjointElement3D3N::InitializeAgainAfterDisturbing()
{
    KRATOS_TRY

    //ShellThinElement3D3N::mSections.clear();
    ShellThinElement3D3N::Initialize();

    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if(rResult.size() != OPT_NUM_DOFS)
        rResult.resize(OPT_NUM_DOFS, false);

    GeometryType & geom = this->GetGeometry();

    for(SizeType i = 0; i < geom.size(); i++)
    {
        int index = i * 6;
        NodeType & iNode = geom[i];

        rResult[index]     = iNode.GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = iNode.GetDof(ADJOINT_ROTATION_X).EquationId();
        rResult[index + 4] = iNode.GetDof(ADJOINT_ROTATION_Y).EquationId();
        rResult[index + 5] = iNode.GetDof(ADJOINT_ROTATION_Z).EquationId();
    }
}

void ShellThinAdjointElement3D3N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    ElementalDofList.resize(0);
    ElementalDofList.reserve(OPT_NUM_DOFS);

    GeometryType & geom = this->GetGeometry();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        NodeType & iNode = geom[i];

        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_X));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Y));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Z));

        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_X));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Y));
        ElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Z));
    }
}

int ShellThinAdjointElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& geom = GetGeometry();

   // verify that the variables are correctly initialized
    if(DISPLACEMENT.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
    if(ROTATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered","");
    if(SHELL_CROSS_SECTION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"SHELL_CROSS_SECTION has Key zero! (check if the application is correctly registered","");
    if(THICKNESS.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
    if(CONSTITUTIVE_LAW.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");


    // check properties
    if(this->pGetProperties() == NULL)
        KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

    const PropertiesType & props = this->GetProperties();    

    if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
    {
        const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
        if(section == NULL)
            KRATOS_THROW_ERROR(std::logic_error,"SHELL_CROSS_SECTION not provided for element ",this->Id());

        section->Check(props, geom, rCurrentProcessInfo);
    }
    else // ... allow the automatic creation of a homogeneous section from a material and a thickness
    {
        if(!props.Has(CONSTITUTIVE_LAW))
            KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());
        const ConstitutiveLaw::Pointer& claw = props[CONSTITUTIVE_LAW];
        if(claw == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());

        if(!props.Has(THICKNESS))
            KRATOS_THROW_ERROR(std::logic_error, "THICKNESS not provided for element ", this->Id());
        if(props[THICKNESS] <= 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "wrong THICKNESS value provided for element ", this->Id());

        ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thin);
        dummySection->Check(props, geom, rCurrentProcessInfo);
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

void ShellThinAdjointElement3D3N::GetValuesVector(Vector& values, int Step)
{
    if(values.size() != OPT_NUM_DOFS)
        values.resize(OPT_NUM_DOFS, false); 

    const GeometryType & geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        const NodeType & iNode = geom[i];
        const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ADJOINT_ROTATION, Step);

        int index = i*6;
        values[index]     = disp[0];
        values[index + 1] = disp[1];
        values[index + 2] = disp[2];

        values[index + 3] = rot[0];
        values[index + 4] = rot[1];
        values[index + 5] = rot[2];
    }
}


double ShellThinAdjointElement3D3N::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
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

double ShellThinAdjointElement3D3N::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
	KRATOS_TRY;

	if(rDesignVariable == SHAPE_SENSITIVITY) 
	{
        double dx, dy, dz, L = 0.0;
   
        dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = this->GetGeometry()[2].X0() - this->GetGeometry()[1].X0();
		dy = this->GetGeometry()[2].Y0() - this->GetGeometry()[1].Y0();
		dz = this->GetGeometry()[2].Z0() - this->GetGeometry()[1].Z0();
		L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = this->GetGeometry()[2].X0() - this->GetGeometry()[0].X0();
		dy = this->GetGeometry()[2].Y0() - this->GetGeometry()[0].Y0();
		dz = this->GetGeometry()[2].Z0() - this->GetGeometry()[0].Z0();
		L += sqrt(dx*dx + dy*dy + dz*dz);
        L /= 3.0;
        
		return L;
	}
	else
		return 1.0;

	KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, 
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
     
        ShellThinElement3D3N::ResetSections();
        ShellThinElement3D3N::Initialize();

		// Compute RHS after disturbance
		this->CalculateRightHandSide(RHS_dist, testProcessInfo); 

		// Compute derivative of RHS w.r.t. design variable with finite differences
		RHS_dist -= RHS_undist;
		RHS_dist /= delta;
		for(unsigned int i = 0; i < RHS_dist.size(); i++)
			rOutput(0, i) = RHS_dist[i];
	
        // Give element original properties back
        this->SetProperties(p_global_properties);
        ShellThinElement3D3N::ResetSections();
        ShellThinElement3D3N::Initialize();
        this->CalculateRightHandSide(RHS_dist, testProcessInfo); 
       	
	}
    else
        rOutput.clear();
	
	KRATOS_CATCH("")

}                                            
	
void ShellThinAdjointElement3D3N::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput, 
											const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

		// define working variables
		Vector RHS_undist;
		Vector RHS_dist;
		ProcessInfo testProcessInfo = rCurrentProcessInfo;

		// Get disturbance measure
        double delta= this->GetValue(DISTURBANCE_MEASURE); 	
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
				//end: derive w.r.t. z-coordinate-----------------------------------------------------

                this->CalculateRightHandSide(RHS_dist, testProcessInfo);

			}// end loop over element nodes
		}
        else
			KRATOS_ERROR << "Unsupported design variable!" << std::endl;  

		KRATOS_CATCH("")

}

void ShellThinAdjointElement3D3N::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

	if(rVariable == STRESS_ON_GP)
	{
	    std::string traced_stress_type = this->GetValue(TRACED_STRESS_TYPE);

        const char item_1 = traced_stress_type.at(0);
        const char item_2 = traced_stress_type.at(1);
        const char item_3 = traced_stress_type.at(2);
        int direction_1 = 0;
        int direction_2 = 0;   
        std::vector<Matrix> stress_vector;
  
        if(item_1 == 'M') 
            ShellThinElement3D3N::GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
        else if(item_1 == 'F') 
            ShellThinElement3D3N::GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);
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

        if(item_3 == 'X')  
            direction_2 = 0; 
        else if(item_3 == 'Y')  
            direction_2 = 1; 
        else if(item_3 == 'Z')  
            direction_2 = 2;   
        else 
            KRATOS_ERROR << "Invalid stress type! " << traced_stress_type << (" is not supported!")  << std::endl;        

        rOutput.resize(OPT_NUM_GP);   
        for(size_t i = 0; i < OPT_NUM_GP; i++)
        {
            rOutput(i) = stress_vector[i](direction_1, direction_2);
        }

    }
    else
    {
        rOutput.resize(OPT_NUM_GP);  
        rOutput.clear(); 
    }

    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput, 
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;                           
           
	if(rVariable == STRESS_DISP_DERIV_ON_GP)   
	{
       this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DV_DERIV_ON_GP)
    {
        std::string design_variable_name = this->GetValue( DESIGN_VARIABLE_NAME );	

        if (KratosComponents<Variable<double>>::Has(design_variable_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);    
        }      
    }
    else
	{
		rOutput.clear();
	}
      
    KRATOS_CATCH("")
}

void ShellThinAdjointElement3D3N::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable, 
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    Vector stress_vector_undist;
    Vector stress_vector_dist;
    ProcessInfo copy_process_info = rCurrentProcessInfo;
    DofsVectorType element_dof_list;
    double dist_measure  = 0.0;
    double original_value = 0.0;

	// Get disturbance measure
    double eta = this->GetValue(DISTURBANCE_MEASURE); 	

    ShellThinElement3D3N::GetDofList(element_dof_list, copy_process_info);

    this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

    rOutput.resize(OPT_NUM_DOFS, OPT_NUM_GP);
    rOutput.clear();
    for(size_t i = 0; i < OPT_NUM_DOFS; i++)
    {
        dist_measure =  eta * element_dof_list[i]->GetSolutionStepValue();
 
        if(fabs(dist_measure) < 1e-20 )
            continue;

        original_value = element_dof_list[i]->GetSolutionStepValue();  
        element_dof_list[i]->GetSolutionStepValue() += dist_measure;

  
        this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

        for(size_t j = 0; j < OPT_NUM_GP; j++)
        {
            stress_vector_dist[j] -= stress_vector_undist[j];
            stress_vector_dist[j] /= dist_measure;
            rOutput(i,j) = stress_vector_dist[j];
        }   

        element_dof_list[i]->GetSolutionStepValue()  = original_value; 
        original_value = 0.0;
    
        stress_vector_dist.clear();
    }

    KRATOS_CATCH("")
}   

void ShellThinAdjointElement3D3N::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, 
                                                const Variable<Vector>& rStressVariable, Matrix& rOutput, 
											    const ProcessInfo& rCurrentProcessInfo) 
{
    KRATOS_TRY;

        // Define working variables
		Vector stress_vector_undist;
		Vector stress_vector_dist;

        // Compute stress on GP before disturbance
		this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        // Get disturbance measure
        double delta= this->GetValue(DISTURBANCE_MEASURE); 	
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
	    delta *= correction_factor;	

        rOutput.resize(1, OPT_NUM_GP);

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

            ShellThinElement3D3N::ResetSections();
            ShellThinElement3D3N::Initialize();

			// Compute stress on GP after disturbance
		    this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

			// Compute derivative of stress w.r.t. design variable with finite differences
			stress_vector_dist  -= stress_vector_undist;
			stress_vector_dist  /= delta;

			for(size_t j = 0; j < OPT_NUM_GP; j++)
			    rOutput(0, j) = stress_vector_dist[j];
		
            // Give element original properties back
            this->SetProperties(p_global_properties);

            ShellThinElement3D3N::ResetSections();
            ShellThinElement3D3N::Initialize();
          
		}
        else
         rOutput.clear();

    KRATOS_CATCH("")
} 

void ShellThinAdjointElement3D3N::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable, 
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // define working variables
	Vector stress_vector_undist;
	Vector stress_vector_dist;
	
    // Get disturbance measure
    double delta= this->GetValue(DISTURBANCE_MEASURE); 	
    double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
	delta *= correction_factor;	

	if(rDesignVariable == SHAPE_SENSITIVITY) 
	{
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
 
		rOutput.resize(dimension * number_of_nodes, OPT_NUM_GP);
     
		// Compute stress on GP before disturbance
	    this->Calculate(rStressVariable, stress_vector_undist, rCurrentProcessInfo);

        //TODO: look that this works also for parallel computing
		for(int j = 0; j < number_of_nodes; j++)
		{
			//begin: derive w.r.t. x-coordinate---------------------------------------------------
			// disturb the design variable
			this->GetGeometry()[j].X0() += delta;

			// Compute stress on GP after disturbance
			this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

			// Compute derivative of stress w.r.t. design variable with finite differences
			stress_vector_dist  -= stress_vector_undist;
			stress_vector_dist  /= delta;

			for(size_t i = 0; i < OPT_NUM_GP; i++)
				rOutput( (0 + j*dimension), i) = stress_vector_dist[i]; 

			// Reset pertubed vector
			stress_vector_dist = Vector(0);

			// undisturb the design variable
			this->GetGeometry()[j].X0() -= delta;
			//end: derive w.r.t. x-coordinate-----------------------------------------------------

			//begin: derive w.r.t. y-coordinate---------------------------------------------------
			// disturb the design variable
			this->GetGeometry()[j].Y0() += delta;

			// Compute stress on GP after disturbance
			this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

			// Compute derivative of stress w.r.t. design variable with finite differences
			stress_vector_dist  -= stress_vector_undist;
			stress_vector_dist  /= delta;

			for(size_t i = 0; i < OPT_NUM_GP; i++)
				rOutput((1 + j*dimension),i) = stress_vector_dist[i]; 

			// Reset pertubed vector
			stress_vector_dist = Vector(0);

			// undisturb the design variable
			this->GetGeometry()[j].Y0() -= delta;
			//end: derive w.r.t. y-coordinate-----------------------------------------------------

			//begin: derive w.r.t. z-coordinate---------------------------------------------------
			// disturb the design variable
			this->GetGeometry()[j].Z0() += delta;

			// Compute stress on GP after disturbance
			this->Calculate(rStressVariable, stress_vector_dist, rCurrentProcessInfo);

			// Compute derivative of stress w.r.t. design variable with finite differences
			stress_vector_dist  -= stress_vector_undist;
			stress_vector_dist  /= delta;

			for(size_t i = 0; i < OPT_NUM_GP; i++)
			    rOutput((2 + j*dimension),i) = stress_vector_dist[i]; 

			// Reset pertubed vector
			stress_vector_dist = Vector(0);

			// undisturb the design variable
			this->GetGeometry()[j].Z0() -= delta;
			//end: derive w.r.t. z-coordinate-----------------------------------------------------

		}// end loop over element nodes
	}
    else
		KRATOS_ERROR << "Unsupported design variable!" << std::endl;  

    KRATOS_CATCH("")
}                                            

void ShellThinAdjointElement3D3N::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, 
        ProcessInfo& rCurrentProcessInfo)
{
    Vector dummy;
    ShellThinElement3D3N::CalculateLocalSystem(rLeftHandSideMatrix, dummy, rCurrentProcessInfo);
}        
// =====================================================================================
//
// Class ShellThinAdjointElement3D3N - Results on Gauss Points
//
// =====================================================================================

void ShellThinAdjointElement3D3N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
					      std::vector<double>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY;
		
	if(this->Has(rVariable))
	{
		// Get result value for output
		double output_value = this->GetValue(rVariable);

		// Resize Output
        if(rOutput.size() != OPT_NUM_GP)
            rOutput.resize(OPT_NUM_GP);

		// Write scalar result value on all Gauss-Points
		for(int i = 0; i < OPT_NUM_GP; i++)
			rOutput[i] = output_value; 

        OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rOutput);    
	}
	else
        KRATOS_ERROR << "Unsupported output variable." << std::endl;



	KRATOS_CATCH("")

}

void ShellThinAdjointElement3D3N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
				    std::vector<double>& rValues,
					const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY;
	this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	KRATOS_CATCH("")
}


// =====================================================================================
//
// Class ShellThinAdjointElement3D3N - Serialization
//
// =====================================================================================

void ShellThinAdjointElement3D3N::save(Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  ShellThinElement3D3N );
}

void ShellThinAdjointElement3D3N::load(Serializer& rSerializer) 
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ShellThinElement3D3N );

}

}



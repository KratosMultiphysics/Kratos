// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
//
// ==============================================================================

// Application includes
#include "small_displacement_simp_element.hpp"

#include "topology_optimization_application.h"


namespace Kratos
{

// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry )
: SmallDisplacementElement( NewId, pGeometry )
{
	//DO NOT ADD DOFS HERE!!!
}


// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
: SmallDisplacementElement( NewId, pGeometry, pProperties )
{
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


// =============================================================================================================================================
// COPY CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( SmallDisplacementSIMPElement const& rOther)
:SmallDisplacementElement(rOther)
{
}


// =============================================================================================================================================
// OPERATIONS
// =============================================================================================================================================

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
	return Element::Pointer( new SmallDisplacementSIMPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


// =============================================================================================================================================
// CLONE
// =============================================================================================================================================

Element::Pointer SmallDisplacementSIMPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

	SmallDisplacementSIMPElement NewElement ( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

	//-----------//

	NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

	if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
	{
		NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

		if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
	}


	return Element::Pointer( new SmallDisplacementSIMPElement(NewElement) );
}

// =============================================================================================================================================
// DESTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::~SmallDisplacementSIMPElement()
{
}


// =============================================================================================================================================
// STARTING / ENDING METHODS
// =============================================================================================================================================

//************************************************************************************
//************************************************************************************

void SmallDisplacementSIMPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo )
{
	KRATOS_TRY

	// Additional part for post-processing of the topology optimized model part
	if (rVariable == X_PHYS)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);

	// From original SmallDisplacementElement
	else if (rVariable == VON_MISES_STRESS)
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	else {

		const unsigned int& integration_points_number = GetGeometry()
            				.IntegrationPointsNumber(mThisIntegrationMethod);

		if (rValues.size() != integration_points_number)
			rValues.resize(integration_points_number);

		for (unsigned int ii = 0; ii < integration_points_number; ii++) {
			rValues[ii] = 0.0;
			rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable,
					rValues[ii]);
		}

	}

	KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementSIMPElement::Calculate(const Variable<double> &rVariable, double &rOutput, const ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	if (rVariable == DCDX || rVariable == LOCAL_STRAIN_ENERGY)
	{
		// Get values
		double E_min     = this->GetValue(E_MIN);
		double E_initial = this->GetValue(E_0);
		double E_current = this->GetValue(YOUNG_MODULUS);
		double penalty   = this->GetValue(PENAL);
		double x_phys    = this->GetValue(X_PHYS);

		// Get element stiffness matrix and modify it with the factor that considers the adjusted Youngs Modulus according to the SIMP method
		// Note that Ke0 is computed based on the originally provided Youngs Modulus in the .mdpa-file
		MatrixType Ke0 = Matrix();
		this->CalculateLeftHandSide(Ke0, const_cast <ProcessInfo&>(rCurrentProcessInfo));
		double E_new     = (E_min + pow(x_phys, penalty) * (E_initial - E_min));
		double factor    = (1/E_current)*E_new;
		MatrixType Ke = Ke0 * factor;

		// Loop through nodes of elements and create elemental displacement vector "ue"
		Element::GeometryType& rGeom = this->GetGeometry();
		unsigned int NumNodes = rGeom.PointsNumber();  //NumNodes=8

		// Resize "ue" according to element type
		Vector ue;
		ue.resize(NumNodes * 3);

		// Set the displacement obtained from the FE-Analysis
		for (unsigned int node_i = 0; node_i < NumNodes; node_i++) {
			array_1d<double, 3> &CurrentDisplacement = rGeom[node_i].FastGetSolutionStepValue(DISPLACEMENT);
			ue[3 * node_i + 0] = CurrentDisplacement[0];
			ue[3 * node_i + 1] = CurrentDisplacement[1];
			ue[3 * node_i + 2] = CurrentDisplacement[2];
		}

		// Calculate trans(ue)*Ke0*ue
		Vector intermediateVector;
		intermediateVector.resize(NumNodes * 3);
		intermediateVector = prod(trans(ue), Ke0);
		double ue_Ke0_ue = inner_prod(intermediateVector, ue);

		if (rVariable == DCDX)
		{
			// Calculation of the compliance sensitivities DCDX
			double dcdx = (-penalty) * (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
			this->SetValue(DCDX, dcdx);
		}
		if (rVariable == LOCAL_STRAIN_ENERGY)
		{
			// Calculation of the local strain energy (requires Ke)
			double local_strain_energy = factor * ue_Ke0_ue;
			this->SetValue(LOCAL_STRAIN_ENERGY, local_strain_energy);
		}

	} else if (rVariable == DVDX) {
		// Calculation of the volume sensitivities DVDX
		this->SetValue(DVDX, 1.0);
	}

	KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementSIMPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	// From original SmallDisplacementElement
	const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

	if (rOutput.size() != integration_points_number)
		rOutput.resize(integration_points_number, false);

	if (rVariable == VON_MISES_STRESS) {
		//create and initialize element variables:
		GeneralVariables Variables;
		this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

		//create constitutive law parameters:
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(),
				rCurrentProcessInfo);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();

		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

		for (unsigned int PointNumber = 0;
				PointNumber < mConstitutiveLawVector.size(); PointNumber++) {
			//compute element kinematics B, F, DN_DX ...
			this->CalculateKinematics(Variables, PointNumber);

			//set general variables to constitutivelaw parameters
			this->SetGeneralVariables(Variables, Values, PointNumber);

			//call the constitutive law to update material variables
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(
					Values);

			ComparisonUtilities EquivalentStress;
			rOutput[PointNumber] = EquivalentStress.CalculateVonMises(
					Variables.StressVector);
		}
	}

	// Additional part for post-processing of the topology optimized model part
	else if (rVariable == X_PHYS)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
			rOutput[PointNumber] = this->GetValue(X_PHYS);
	}
	else
	{
		for (unsigned int ii = 0; ii < integration_points_number; ii++)
			rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
	}

	KRATOS_CATCH( "" )
}


// =============================================================================================================================================
// =============================================================================================================================================


void SmallDisplacementSIMPElement::save( Serializer& rSerializer ) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
}

void SmallDisplacementSIMPElement::load( Serializer& rSerializer )
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
}

} // Namespace Kratos

// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

// Application includes

#include "small_displacement_simp_element.h"
#include "comparison_utilities.h" ///nicht mehr vorhanden auf momentanem Stand von Kratos --> muss neue Verknüpfung gefunden werden, vorerst wird diese verwendet. Auch soli_mechanics_math_utilities
#include "topology_optimization_application.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseType( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const 
{
	return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
	return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSIMPElement::~SmallDisplacementSIMPElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementSIMPElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementSIMPElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementSIMPElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


// =============================================================================================================================================
// STARTING / ENDING METHODS
// =============================================================================================================================================

//************************************************************************************
//************************************************************************************


/*  void SmallDisplacementSIMPElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo ) */
	
void SmallDisplacementSIMPElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	// Additional part for post-processing of the topology optimized model part
	if (rVariable == X_PHYS)
		//CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		;
	// From original SmallDisplacementElement
	else if (rVariable == VON_MISES_STRESS)
		//CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		;
	else {

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		if (rValues.size() != integration_points.size())
			rValues.resize(integration_points.size());

		for ( SizeType ii = 0; ii < integration_points.size(); ii++ )
      	{
		rValues[ii] = 0.0;
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
      	}


	}

	// From original SmallDisplacementElement
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);


	if (rValues.size() != integration_points.size()) /// Hier eventuell .size() einbauen
		rValues.resize(integration_points.size(), false);

	if (rVariable == VON_MISES_STRESS) {

		//create and initialize element variables:
		const SizeType number_of_nodes = GetGeometry().size();
		const SizeType dimension = GetGeometry().WorkingSpaceDimension();
		const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

		KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
		ConstitutiveVariables this_constitutive_variables(strain_size);

		//create constitutive law parameters:
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);


		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		for (unsigned int PointNumber = 0;
				PointNumber < mConstitutiveLawVector.size(); PointNumber++) {
			//compute element kinematics B, F, DN_DX ...
			this->CalculateKinematicVariables( this_kinematic_variables, PointNumber,  this->GetIntegrationMethod());
			//set general variables to constitutivelaw parameters
			this->SetElementData(this_kinematic_variables, Values, PointNumber); 
			//call the constitutive law to update material variables
			mConstitutiveLawVector[PointNumber]->InitializeMaterialResponse(
					Values, GetStressMeasure());
			ComparisonUtilities EquivalentStress;
			rValues[PointNumber] = EquivalentStress.CalculateVonMises(
					this_constitutive_variables.StressVector);
		}
	} 

	// Additional part for post-processing of the topology optimized model part
	else if (rVariable == X_PHYS)
	{
		for (SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
			rValues[PointNumber] = this->GetValue(X_PHYS);
	}
	else
	{
		for (SizeType ii = 0; ii < integration_points.size(); ii++ )
			rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rValues[ii]);
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

		///Normalize the youngs modulus
		double factor    = (1/E_current)*E_new; 
		MatrixType Ke = Ke0 * factor;

		// Loop through nodes of elements and create elemental displacement vector "ue"
		Element::GeometryType& rGeom = this->GetGeometry();
		unsigned int NumNodes = rGeom.PointsNumber(); 

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
			// Do they have to be normalized?
			double dcdx = (-penalty)* (E_initial - E_min) * pow(x_phys, penalty - 1) * ue_Ke0_ue;
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
 
// =============================================================================================================================================
// =============================================================================================================================================

void SmallDisplacementSIMPElement::SetElementData(const KinematicVariables& rThisKinematicVariables,
                                              ConstitutiveLaw::Parameters& rValues,
                                              const int & rPointNumber)
{
    KRATOS_TRY
	
	const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

	ConstitutiveVariables this_constitutive_variables(strain_size);

    rValues.SetStrainVector(this_constitutive_variables.StrainVector);
    rValues.SetStressVector(this_constitutive_variables.StressVector);
    rValues.SetConstitutiveMatrix(this_constitutive_variables.D);
    rValues.SetShapeFunctionsDerivatives(rThisKinematicVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N);

    if(rThisKinematicVariables.detJ0<0)
      {
	KRATOS_ERROR << " (small displacement) ELEMENT INVERTED |J|<0 : " << rThisKinematicVariables.detJ0 << std::endl;
      }

    rValues.SetDeterminantF(rThisKinematicVariables.detF);
    rValues.SetDeformationGradientF(rThisKinematicVariables.F);

    KRATOS_CATCH( "" )
}

// =============================================================================================================================================
// =============================================================================================================================================


void SmallDisplacementSIMPElement::save( Serializer& rSerializer ) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

void SmallDisplacementSIMPElement::load( Serializer& rSerializer )
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos

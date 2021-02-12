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

#include "small_displacement_simp_element.hpp"
#include "compairson_utilities.hpp" ///nicht mehr vorhanden auf momentanem Stand von Kratos --> muss neue Verknüpfung gefunden werden, vorerst wird diese verwendet. Auch soli_mechanics_math_utilities
#include "topology_optimization_application.h"
#include "custom_elements/small_displacement.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry )
: SmallDisplacement( NewId, pGeometry )
{
	//DO NOT ADD DOFS HERE!!!
}


// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : SmallDisplacement( NewId, pGeometry, pProperties )
{
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


// =============================================================================================================================================
// COPY CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementSIMPElement::SmallDisplacementSIMPElement( SmallDisplacementSIMPElement const& rOther)
:SmallDisplacement(rOther)
{
}


// =============================================================================================================================================
// OPERATIONS
// =============================================================================================================================================

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

Element::Pointer SmallDisplacementSIMPElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementSIMPElement>( NewId, pGeom, pProperties );
}


// =============================================================================================================================================
// CLONE
// =============================================================================================================================================

Element::Pointer SmallDisplacementSIMPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    SmallDisplacementSIMPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR << " constitutive law not has the correct size small displacement element " << std::endl;
      }


    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(Flags(*this));
	NewElement.SetIntegrationMethod(BaseType::mThisIntegrationMethod);
	NewElement.SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return Kratos::make_intrusive< SmallDisplacementSIMPElement >(NewElement);
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
	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
	const std::size_t number_of_integration_points = integration_points.size();
    const auto& r_geometry = GetGeometry();

	 ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);




    if ( rOutput.size() != number_of_integration_points )
         rOutput.resize( number_of_integration_points );
	
	if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);

	} 
	if (rVariable == VON_MISES_STRESS) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( this_constitutive_variables.StressVector );

                double sigma_equivalent = 0.0;

                if (dimension == 2) {
                    sigma_equivalent = std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                                3*(stress_tensor(0,1) * stress_tensor(1,0));
                } else {
                    sigma_equivalent = 0.5*(std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                            std::pow((stress_tensor(1,1) - stress_tensor(2,2)), 2.0) +
                                            std::pow((stress_tensor(2,2) - stress_tensor(0,0)), 2.0) +
                                                    6*(stress_tensor(0,1) * stress_tensor(1,0) +
                                                        stress_tensor(1,2) * stress_tensor(2,1) +
                                                        stress_tensor(2,0) * stress_tensor(0,2)));
                }

                if( sigma_equivalent < 0.0 )
                    rOutput[point_number] = 0.0;
                else
                    rOutput[point_number] = std::sqrt(sigma_equivalent);
            }
        }
    


///	if (rVariable == VON_MISES_STRESS) {
		//create and initialize element variables:
		// ElementType Variables;
	/* 	this->InitializeElementData(Variables, rCurrentProcessInfo);
		//create constitutive law parameters:
		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(),
				rCurrentProcessInfo);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();

		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

		for (unsigned int PointNumber = 0;
				PointNumber < mConstitutiveLawVector.size(); PointNumber++) {
			//compute element kinematics B, F, DN_DX ...
			this->CalculateKinematicVariables( Variables, PointNumber);

			//set general variables to constitutivelaw parameters
			this->SetElementData(Variables, Values, PointNumber);

			//call the constitutive law to update material variables
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(
					Values);

			ComparisonUtilities EquivalentStress;
			rOutput[PointNumber] = EquivalentStress.CalculateVonMises(
					Variables.StressVector);
		}
	} */

	// Additional part for post-processing of the topology optimized model part
	if (rVariable == X_PHYS)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
			rOutput[PointNumber] = this->GetValue(X_PHYS);
	}
	else
	{
		for (unsigned int ii = 0; ii < number_of_integration_points; ii++)
			rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rOutput[ii]);
	}
	
	KRATOS_CATCH( "" )
}


// =============================================================================================================================================
// =============================================================================================================================================


void SmallDisplacementSIMPElement::save( Serializer& rSerializer ) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacement )
}

void SmallDisplacementSIMPElement::load( Serializer& rSerializer )
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement )
}

} // Namespace Kratos

//
//   Project Name:        
//   Last modified by:    $Author:          
//   Date:                $Date:            
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement() : SmallDisplacementElement() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement( IndexType NewId, GeometryType::Pointer pGeometry ) : SmallDisplacementElement( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallDisplacementThermoMechanicElement::SmallDisplacementThermoMechanicElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : SmallDisplacementElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallDisplacementThermoMechanicElement::~SmallDisplacementThermoMechanicElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallDisplacementThermoMechanicElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallDisplacementThermoMechanicElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == THERMAL_STRESS_TENSOR || rVariable == MECHANICAL_STRESS_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == THERMAL_STRAIN_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
        for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == THERMAL_STRESS_VECTOR  || rVariable == MECHANICAL_STRESS_VECTOR )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        
        if( rVariable == CAUCHY_STRESS_VECTOR)
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        else if(rVariable == THERMAL_STRESS_VECTOR) 
            ConstitutiveLawOptions.Set(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY);
        else if(rVariable == MECHANICAL_STRESS_VECTOR)  
            ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
                rOutput[PointNumber].resize( Variables.StressVector.size(), false );

            rOutput[PointNumber] = Variables.StressVector;
        }
    }
    else if ( rVariable == THERMAL_STRAIN_VECTOR )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
        
        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == THERMAL_STRESS_TENSOR || rVariable == MECHANICAL_STRESS_TENSOR  )
    {
        std::vector<Vector> StressVector;
	
        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
        else if ( rVariable == THERMAL_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( THERMAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( MECHANICAL_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
	    }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == THERMAL_STRAIN_TENSOR)
    {
        std::vector<Vector> StrainVector;
        
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else
            CalculateOnIntegrationPoints( THERMAL_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }

    KRATOS_CATCH( "" )
}

} // Namespace Kratos

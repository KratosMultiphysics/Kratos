//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//


/* Project includes */
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"


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

void SmallDisplacementThermoMechanicElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE); //Note: this is for nonlocal damage
    
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetElementData(Variables,Values,PointNumber);
        
        //call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    this->InitializeNonLinearIteration(rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{    
    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    //Extrapolation variables
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int& Dim  = Geom.WorkingSpaceDimension();
    const unsigned int NumGPoints = Geom.IntegrationPointsNumber( mThisIntegrationMethod );
    unsigned int VoigtSize = 6;
    if(Dim == 2) VoigtSize = 3;
    Matrix StressContainer(NumGPoints,VoigtSize);

    for ( unsigned int PointNumber = 0; PointNumber < NumGPoints; PointNumber++ )
    {

        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetElementData(Variables,Values,PointNumber);
        
        //call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy(Values);
        
        this->SaveGPStress(StressContainer,Variables.StressVector,VoigtSize,PointNumber);
    }
    
    this->ExtrapolateGPStress(StressContainer,Dim,VoigtSize);
}

//----------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::SaveGPStress(Matrix& rStressContainer, const Vector& StressVector, const unsigned int& VoigtSize, const unsigned int& GPoint)
{
    for(unsigned int i = 0; i < VoigtSize; i++)
    {
        rStressContainer(GPoint,i) = StressVector[i];
    }

    /* INFO: (Quadrilateral_2D_4 with GI_GAUSS_2)
     * 
     *                      |S0-0 S1-0 S2-0|
     * rStressContainer =   |S0-1 S1-1 S2-1|
     *                      |S0-2 S1-2 S2-2|
     *                      |S0-3 S1-3 S2-3|
     * 
     * S1-0 = S[1] at GP 0
    */
}

//----------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicElement::ExtrapolateGPStress(const Matrix& StressContainer, const unsigned int& Dim, const unsigned int& VoigtSize)
{
    GeometryType& rGeom = this->GetGeometry();
    //const unsigned int& Dim  = rGeom.WorkingSpaceDimension();
    const unsigned int& NumNodes = rGeom.size();
    const double& Area = rGeom.Area(); // In 3D this is Volume
    
    std::vector<Vector> NodalStressVector(NumNodes); //List with stresses at each node
    std::vector<Matrix> NodalStressTensor(NumNodes);
    
    for(unsigned int Node = 0; Node < NumNodes; Node ++)
    {
        NodalStressVector[Node].resize(VoigtSize);
        NodalStressTensor[Node].resize(Dim,Dim);
    }
    
    if (Dim == 2)
    {    
        if(NumNodes == 3)
        {
            // Triangle_2d_3 with GI_GAUSS_1
            for(unsigned int i = 0; i < 3; i++) //NumNodes
            {
                noalias(NodalStressVector[i]) = row(StressContainer,0)*Area;
                noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);
                
                rGeom[i].SetLock();
                noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
                rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
                rGeom[i].UnSetLock();
            }
        }
        else if(NumNodes == 4)
        {
            // Quadrilateral_2d_4 with GI_GAUSS_2
            BoundedMatrix<double,4,4> ExtrapolationMatrix;
            PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
            
            BoundedMatrix<double,4,3> AuxNodalStress;
            noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

            /* INFO:
             * 
             *                  |S0-0 S1-0 S2-0|
             * AuxNodalStress = |S0-1 S1-1 S2-1|
             *                  |S0-2 S1-2 S2-2|
             *                  |S0-3 S1-3 S2-3|
             * 
             * S1-0 = S[1] at node 0
            */

            for(unsigned int i = 0; i < 4; i++) //TNumNodes
            {
                noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
                noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);
                
                rGeom[i].SetLock();
                noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
                rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
                rGeom[i].UnSetLock();
            }
        }
    }
    else
    {
        if(NumNodes == 4)
        {
            // Tetrahedra_3d_4 with GI_GAUSS_1
            for(unsigned int i = 0; i < 4; i++) //NumNodes
            {
                noalias(NodalStressVector[i]) = row(StressContainer,0)*Area;
                noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);
                
                rGeom[i].SetLock();
                noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
                rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
                rGeom[i].UnSetLock();
            }
        }
        else if(NumNodes == 8)
        {
            // Hexahedra_3d_8 with GI_GAUSS_2
            BoundedMatrix<double,8,8> ExtrapolationMatrix;
            PoroElementUtilities::CalculateExtrapolationMatrix(ExtrapolationMatrix);
            
            BoundedMatrix<double,8,6> AuxNodalStress;
            noalias(AuxNodalStress) = prod(ExtrapolationMatrix,StressContainer);

            for(unsigned int i = 0; i < 8; i++) //TNumNodes
            {
                noalias(NodalStressVector[i]) = row(AuxNodalStress,i)*Area;
                noalias(NodalStressTensor[i]) = MathUtils<double>::StressVectorToTensor(NodalStressVector[i]);
                
                rGeom[i].SetLock();
                noalias(rGeom[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)) += NodalStressTensor[i];
                rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += Area;
                rGeom[i].UnSetLock();
            }
        }
    }
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
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        
        if( rVariable == CAUCHY_STRESS_VECTOR){
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	}
	else if(rVariable == THERMAL_STRESS_VECTOR){
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::THERMAL_RESPONSE_ONLY);
	}
	else if(rVariable == MECHANICAL_STRESS_VECTOR){
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::MECHANICAL_RESPONSE_ONLY);
	}

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

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
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::THERMAL_RESPONSE_ONLY);
	
        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

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
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

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

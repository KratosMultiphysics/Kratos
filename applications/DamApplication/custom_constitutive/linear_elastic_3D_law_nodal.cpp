//
//   Project Name:   
//   Last modified by:    $Author:     
//   Date:                $Date:     
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_constitutive/linear_elastic_3D_law_nodal.hpp"


namespace Kratos
{

//Default Constructor
LinearElastic3DLawNodal::LinearElastic3DLawNodal() : LinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic3DLawNodal::LinearElastic3DLawNodal(const LinearElastic3DLawNodal& rOther) : LinearElastic3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic3DLawNodal::~LinearElastic3DLawNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic3DLawNodal::Clone() const
{
    LinearElastic3DLawNodal::Pointer p_clone(new LinearElastic3DLawNodal(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLawNodal::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    KRATOS_TRY
    
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    Flags& Options = rValues.GetOptions();
    
    Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;

	ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
	ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());

    //1.- Lame constants
    double YoungModulus;
    this->CalculateNodalYoungModulus( ElasticVariables, YoungModulus);
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];


    if(Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
    {
		this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
		
		if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) //TOTAL STRESS
		{
			noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
		}
	}
    else if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) //TOTAL STRESS
    {        
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
    }


    KRATOS_CATCH( "" )
    
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  LinearElastic3DLawNodal::CalculateNodalYoungModulus (const MaterialResponseVariables & rElasticVariables, double & rYoungModulus)
{
    KRATOS_TRY
    
    //1.-Young Modulus from nodes 
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();
    
    rYoungModulus = 0.0;
    
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rYoungModulus += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(NODAL_YOUNG_MODULUS);
    }

    return rYoungModulus;
    
    KRATOS_CATCH( "" )
}


} // Namespace Kratos

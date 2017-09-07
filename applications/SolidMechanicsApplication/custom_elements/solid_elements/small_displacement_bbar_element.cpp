//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               MCaicedo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/solid_elements/small_displacement_bbar_element.hpp"
#include "includes/constitutive_law.h"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement()
  : SmallDisplacementElement()
{
    // DO NOT CALL IT: only needed for Register and Serialization!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(IndexType NewId,
                                                           GeometryType::Pointer pGeometry)
    : SmallDisplacementElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(IndexType NewId,
                                                           GeometryType::Pointer pGeometry,
                                                           PropertiesType::Pointer pProperties)
    : SmallDisplacementElement(NewId, pGeometry, pProperties)
{
    
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(SmallDisplacementBbarElement const& rOther)
    : SmallDisplacementElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SmallDisplacementBbarElement& SmallDisplacementBbarElement::operator=(SmallDisplacementBbarElement const& rOther)
{
    SmallDisplacementElement::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbarElement::Create(IndexType NewId,
                                                      NodesArrayType const& rThisNodes,
                                                      PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SmallDisplacementBbarElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbarElement::Clone(IndexType NewId,
                                                     NodesArrayType const& rThisNodes) const
{
    SmallDisplacementBbarElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if (NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
    {
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

        if (NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber())
            KRATOS_THROW_ERROR(std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size());
    }

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    //-----------//

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new SmallDisplacementBbarElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::~SmallDisplacementBbarElement()
{
}


//************************************************************************************
//************************************************************************************


void SmallDisplacementBbarElement::InitializeElementVariables(ElementVariables& rVariables,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = 4; // added component zz, necessary for plasticity.

    if (dimension == 3)
    {
        voigt_size = 6;
    }

    rVariables.Initialize(voigt_size, dimension, number_of_nodes);

    // needed parameters for consistency with the general constitutive law:
    // small displacements
    rVariables.detF = 1.0;
    rVariables.F = identity_matrix<double>(dimension);

    // set variables including all integration points values

    // reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod));

    // reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod));

    // calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian(rVariables.j, mThisIntegrationMethod);

    // Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    // calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian(rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition);

    // calculate hydrostatic deformation matrix
    CalculateHydrostaticDeformationMatrix(rVariables);

    KRATOS_CATCH( "" )    
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


void SmallDisplacementBbarElement::CalculateKinematics(ElementVariables& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    // Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    // Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix(rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    // Compute cartesian derivatives  [dN/dx_n]
    noalias(rVariables.DN_DX) = prod(DN_De[rPointNumber], InvJ);

    // Set Shape Functions Values for this integration point
    rVariables.N = row(Ncontainer, rPointNumber);

    // Compute the deformation matrix B_bar
    this->CalculateDeformationMatrixBbar(rVariables.B, rVariables.H, rVariables.DN_DX);
    
    // Compute infinitessimal B_bar strain
    this->CalculateInfinitesimalStrainBbar(rVariables.B, rVariables.StrainVector);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBbarElement::CalculateInfinitesimalStrainBbar(const Matrix& rB,
								    Vector& rStrainVector)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension == 2)
    {
        // Infinitesimal Bbar Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            array_1d<double, 3>& Displacement =
                GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rStrainVector[0] += Displacement[0] * rB(0, i * 2) +
                                Displacement[1] * rB(0, i * 2 + 1); // xx
            rStrainVector[1] += Displacement[0] * rB(1, i * 2) +
                                Displacement[1] * rB(1, i * 2 + 1); // yy
            rStrainVector[2] += Displacement[0] * rB(2, i * 2) +
                                Displacement[1] * rB(2, i * 2 + 1); // zz
            rStrainVector[3] += Displacement[0] * rB(3, i * 2) +
                                Displacement[1] * rB(3, i * 2 + 1); // xy
        }
    }
    else if (dimension == 3)
    {

        //Infinitesimal Bbar Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

      
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            array_1d<double, 3>& Displacement =
                GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rStrainVector[0] += Displacement[0] * rB(0, i * 3) +
                                Displacement[1] * rB(0, i * 3 + 1) +
                                Displacement[2] * rB(0, i * 3 + 2); // xx
            rStrainVector[1] += Displacement[0] * rB(1, i * 3) +
                                Displacement[1] * rB(1, i * 3 + 1) +
                                Displacement[2] * rB(1, i * 3 + 2); // yy
            rStrainVector[2] += Displacement[0] * rB(2, i * 3) +
                                Displacement[1] * rB(2, i * 3 + 1) +
                                Displacement[2] * rB(2, i * 3 + 2); // zz
            rStrainVector[3] += Displacement[0] * rB(3, i * 3) +
                                Displacement[1] * rB(3, i * 3 + 1) +
                                Displacement[2] * rB(3, i * 3 + 2); // xy
            rStrainVector[4] += Displacement[0] * rB(4, i * 3) +
                                Displacement[1] * rB(4, i * 3 + 1) +
                                Displacement[2] * rB(4, i * 3 + 2); // yz
            rStrainVector[5] += Displacement[0] * rB(5, i * 3) +
                                Displacement[1] * rB(5, i * 3 + 1) +
                                Displacement[2] * rB(5, i * 3 + 2); // xz
        }
    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument, "something is wrong with the dimension", "")
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
  
void SmallDisplacementBbarElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = 4; // added component zz, necessary for plasticity.

    if (dimension == 3)
    {
        voigt_size = 6;
    }
    
    if (rB.size1() != voigt_size || rB.size2() != dimension * number_of_nodes)
        rB.resize(voigt_size, dimension * number_of_nodes, false);
    
    if (dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 2 * i;

            rB(0,index + 0) = rDN_DX(i,0);
            rB(0,index + 1) = 0.0;

	    rB(1,index + 0) = 0.0;
            rB(1,index + 1) = rDN_DX(i,1);

	    rB(2,index + 0) = 0.0;
            rB(2,index + 1) = 0.0;

	    rB(3,index + 0) = rDN_DX(i,1);
            rB(3,index + 1) = rDN_DX(i,0);
        }
    }
    else if (dimension == 3)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 3 * i;

            rB(0,index + 0) = rDN_DX(i,0);
            rB(1,index + 1) = rDN_DX(i,1);
            rB(2,index + 2) = rDN_DX(i,2);

	    rB(3,index + 0) = rDN_DX(i,1);
            rB(3,index + 1) = rDN_DX(i,0);

	    rB(4,index + 1) = rDN_DX(i,2);
            rB(4,index + 2) = rDN_DX(i,1);
	    
            rB(5,index + 0) = rDN_DX(i,2);
            rB(5,index + 2) = rDN_DX(i,0);
        }
    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument, "something is wrong with the dimension", "")
    }
    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbarElement::CalculateDeformationMatrixBbar(Matrix& rB,
                                                                  const Matrix& rH,
                                                                  const Matrix& rDN_DX) // JLM debe entrar rBh
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = 4; // added component zz, necessary for plasticity.

    if (dimension == 3)
    {
        voigt_size = 6;
    }

    // Inicializar Auxiliar B matrix
    Matrix B;
    this->CalculateDeformationMatrix(B, rDN_DX);
    
    if ( rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes )
      rB.resize(voigt_size, dimension*number_of_nodes, false );

    double athird = 1.0/3.0;
    
    if (dimension == 2)
    {
        rB(0,0) =  2.0 * athird * B(0,0);
        rB(0,1) = -athird * B(1,1);

        rB(0,2) =  2.0 * athird * B(0,2);
        rB(0,3) = -athird * B(1,3);
        rB(0,4) =  2.0 * athird * B(0,4);
        rB(0,5) = -athird * B(1,5);
        rB(0,6) =  2.0 * athird * B(0,6);
        rB(0,7) = -athird * B(1,7);

        rB(1,0) = -athird * B(0,0);
        rB(1,1) =  2.0 * athird * B(1,1);

        rB(1,2) = -athird * B(0,2);
        rB(1,3) =  2.0 * athird * B(1,3);
        rB(1,4) = -athird * B(0,4);
        rB(1,5) =  2.0 * athird * B(1,5);
        rB(1,6) = -athird * B(0,6);
        rB(1,7) =  2.0 * athird * B(1,7);

        rB(2,0) = -athird * B(0,0);
        rB(2,1) = -athird * B(1,1);

        rB(2,2) = -athird * B(0,2);
        rB(2,3) = -athird * B(1,3);
        rB(2,4) = -athird * B(0,4);
        rB(2,5) = -athird * B(1,5);
        rB(2,6) = -athird * B(0,6);
        rB(2,7) = -athird * B(1,7);

        for (unsigned int i = 0; i < number_of_nodes * dimension; i++)
        {
            // unsigned int index = 2 * i;
   	    rB(0,i) += athird * rH(i,0);
            rB(1,i) += athird * rH(i,0);
            rB(2,i) += athird * rH(i,0);
            rB(3,i)  = B(3,i);
        }
    }
    else if (dimension == 3)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = 3 * i;

            rB(0,index + 0) = 2.0 * athird * B(0,index + 0);
            rB(1,index + 0) = -athird * B(0,index + 0);
            rB(2,index + 0) = -athird * B(0,index + 0);
            rB(0,index + 1) = -athird * B(1,index + 1);
            rB(1,index + 1) = 2.0 * athird * B(1,index + 1);
            rB(2,index + 1) = -athird * B(1,index + 1);
            rB(0,index + 2) = -athird * B(2,index + 2);
            rB(1,index + 2) = -athird * B(2,index + 2);
            rB(2,index + 2) = 2.0 * athird * B(2,index + 2);
        }
        for (unsigned int i = 0; i < number_of_nodes * dimension; i++)
        {
            rB(0,i) += athird * rH(i,0);
            rB(1,i) += athird * rH(i,0);
            rB(2,i) += athird * rH(i,0);
            rB(3,i) = B(3,i);
            rB(4,i) = B(4,i);
            rB(5,i) = B(5,i);
        }
    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "something is wrong with the dimension", "")
    }


    KRATOS_CATCH("")
}

void SmallDisplacementBbarElement::CalculateHydrostaticDeformationMatrix(ElementVariables& rVariables)
{
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // if (rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes)
    //    rB.resize(voigt_size, dimension*number_of_nodes, false);

    rVariables.H.resize(dimension * number_of_nodes, 1, false);
    rVariables.H.clear();

    // reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    double GeometrySize = 0.0;
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
      {
	const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
	    
	Matrix InvJ;
	MathUtils<double>::InvertMatrix(rVariables.J[PointNumber], InvJ, rVariables.detJ);
	noalias(rVariables.DN_DX) = prod(DN_De[PointNumber], InvJ);

	this->CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX);

	double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;
	GeometrySize += IntegrationWeight;

	unsigned int index = 0;
	for (unsigned int i = 0; i < number_of_nodes; i++)
	  {
	    for (unsigned int j= 0; j < dimension; j++)
	      {
		index = i * dimension + j;
		rVariables.H(index,0) += rVariables.B(j, index) * IntegrationWeight;
	      }
	  }
      }

    rVariables.H /= GeometrySize;

    KRATOS_CATCH("")
}

void SmallDisplacementBbarElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

void SmallDisplacementBbarElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

} // Namespace Kratos

/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Zhiming $
//   Date:                $Date: 2014-06 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "lagrange_multiplier_condition_0.h"
#include "utilities/math_utils.h"
#include "geometries/point_2d.h"
//#include "custom_elements/MLSparticle.h"
#include "custom_utilities/compute_mls_shape_functions_utility.h"



namespace Kratos
{
	
LagrangeMultiplierCondition2D0::LagrangeMultiplierCondition2D0(IndexType NewId, GeometryType::Pointer pGeometry)
	: Condition(NewId, pGeometry)
{
}

LagrangeMultiplierCondition2D0::LagrangeMultiplierCondition2D0(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: Condition(NewId, pGeometry, pProperties)
{
}
	

Condition::Pointer LagrangeMultiplierCondition2D0::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LagrangeMultiplierCondition2D0(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

LagrangeMultiplierCondition2D0::~LagrangeMultiplierCondition2D0()
{
}

void LagrangeMultiplierCondition2D0::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
   	KRATOS_TRY;

    const unsigned int number_of_nodes = this->GetGeometry().size();

    //Geometry< Node<3> >& geom = this->GetGeometry();
    //const unsigned int number_of_nodes = geom.size();

/*
    Matrix TempCoordinates(number_of_nodes,2);
    const double h = this->GetValue(EFFECTIVE_RADIUS);
    //const double& A = this->GetValue(GAUSS_AREA);


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {

        row(TempCoordinates,i) = this->GetGeometry()[i].Coordinates();
        //row(TempCoordinates,i) = this->GetGeometry()[i].GetInitialPosition()
               // + this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

    }

     //const array_1d<double,2>& xg = row(TempCoordinates,0);
     //const array_1d<double,3>& xg = this->GetValue(TEMP_POS);

     const array_1d<double,3> xg = row(TempCoordinates,0);
     Vector Nc( number_of_nodes);
     Matrix DN_DXc( number_of_nodes, 2 );
     LinearMLSKernel::ComputeMLSKernel(Nc, DN_DXc, TempCoordinates, xg, h);
*/

     Vector& Nc = this->GetValue(SHAPE_FUNCTIONS);


    unsigned int dim = 2;
    unsigned int matsize = number_of_nodes*dim+2;

    if(rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize,matsize,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize,matsize);
    if(rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false);
    noalias(rRightHandSideVector) = ZeroVector(matsize);

    Vector currentValues(matsize);



    for ( unsigned int node = 0; node < number_of_nodes; node++ )
    {

        int index = node * dim;
        currentValues(index) = GetGeometry()[node].FastGetSolutionStepValue(DISPLACEMENT_X);
        currentValues(index + 1) = GetGeometry()[node].FastGetSolutionStepValue(DISPLACEMENT_Y);
     }
    currentValues(number_of_nodes*dim) = GetGeometry()[0].FastGetSolutionStepValue(LAGRANGE_DISPLACEMENT_X);
    currentValues(number_of_nodes*dim + 1) = GetGeometry()[0].FastGetSolutionStepValue(LAGRANGE_DISPLACEMENT_Y);






    for( unsigned int n=0; n<number_of_nodes; n++ )
    {
        //filling the column
        rLeftHandSideMatrix(n*dim, number_of_nodes*dim) = Nc[n]; //1.0;
        rLeftHandSideMatrix(n*dim+1, number_of_nodes*dim+1) = Nc[n];

        //filling the row
        rLeftHandSideMatrix(number_of_nodes*dim, n*dim ) = Nc[n];
        rLeftHandSideMatrix(number_of_nodes*dim+1, n*dim+1) = Nc[n];


    }

    //rLeftHandSideMatrix *= Nc[0];

     //KRATOS_WATCH(GetGeometry()[0].X());
     //KRATOS_WATCH(GetGeometry()[0].Y());

    rRightHandSideVector[number_of_nodes*dim] = 0.0;
    rRightHandSideVector[number_of_nodes*dim + 1] =0.0;
    //this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // form residual

    noalias(rRightHandSideVector) -= prod( rLeftHandSideMatrix, currentValues );

	
    KRATOS_CATCH("");
       
       
}

void LagrangeMultiplierCondition2D0::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

            const unsigned int number_of_nodes = this->GetGeometry().size();

          //Geometry< Node<3> >& geom = this->GetGeometry();
           //onst unsigned int number_of_nodes = geom.size();

            unsigned int dim = 2;

            unsigned int matsize = number_of_nodes*dim+2;

            if(rRightHandSideVector.size() != matsize)
                rRightHandSideVector.resize(matsize,false);
            noalias(rRightHandSideVector) = ZeroVector(matsize);


            //MatrixType temp = Matrix();
            //CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);


    KRATOS_CATCH("");
}

void LagrangeMultiplierCondition2D0::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;

            unsigned int number_of_nodes = this->GetGeometry().size();
            //Geometry< Node<3> >& geom = this->GetGeometry();
            //const unsigned int number_of_nodes = geom.size();

            int dim = 2;
            unsigned int dim2 = number_of_nodes * dim + 2;

            if ( rResult.size() != dim2 )
                rResult.resize( dim2, false );

            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                int index = i * dim;
                rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

            }
            rResult[number_of_nodes * dim ] = GetGeometry()[0].GetDof( LAGRANGE_DISPLACEMENT_X ).EquationId();
            rResult[number_of_nodes * dim + 1] = GetGeometry()[0].GetDof( LAGRANGE_DISPLACEMENT_Y ).EquationId();




    KRATOS_CATCH("");
}

void LagrangeMultiplierCondition2D0::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;

            ConditionalDofList.resize( 0 );


            unsigned int number_of_nodes = this->GetGeometry().size();
            //Geometry< Node<3> >& geom = this->GetGeometry();
            //const unsigned int number_of_nodes = geom.size();
             //KRATOS_WATCH(TempNeighbours)


            //for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
                ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

            }

            ConditionalDofList.push_back( GetGeometry()[0].pGetDof( LAGRANGE_DISPLACEMENT_X ) );
            ConditionalDofList.push_back( GetGeometry()[0].pGetDof( LAGRANGE_DISPLACEMENT_Y ) );




	KRATOS_CATCH("")
  
}

int LagrangeMultiplierCondition2D0::Check(const ProcessInfo& rCurrentProcessInfo)
{

	KRATOS_TRY

	

    KRATOS_CATCH("");
	return 0;
}

}




/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author: kkazem $
//   Date:                $Date: 2008-02-14 09:41:09 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/environment_contact_3d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
EnvironmentContact3D::EnvironmentContact3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
EnvironmentContact3D::EnvironmentContact3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

}

Condition::Pointer EnvironmentContact3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new EnvironmentContact3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

EnvironmentContact3D::~EnvironmentContact3D()
{
}


//************************************************************************************
//************************************************************************************
void EnvironmentContact3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    MatrixType temp = Matrix();
// 		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

}

//************************************************************************************
//************************************************************************************
void EnvironmentContact3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    unsigned int nodes_number = GetGeometry().PointsNumber();
    unsigned int Msize = nodes_number;
    
    if(rLeftHandSideMatrix.size1() != Msize)
    {
        rLeftHandSideMatrix.resize(Msize,Msize,false);
        rRightHandSideVector.resize(Msize,false);

    }

    noalias(rRightHandSideVector) = ZeroVector(Msize);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(Msize, Msize);
    //calculate area
    double area = GetGeometry().Area();

    //take thermal properties
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();

    array_1d<double,3> NN;
    double amb_T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];
	
    //3 Gauss points coordinates
    bounded_matrix<double, 3, 3 > GaussCrd = ZeroMatrix(3, 3);
    GaussCrd(0,0) = 2.0/3.0; GaussCrd(0,1) = 1.0/6.0; GaussCrd(0,2) = 1.0/6.0; 
    GaussCrd(1,0) = 1.0/6.0; GaussCrd(1,1) = 2.0/3.0; GaussCrd(1,2) = 1.0/6.0; 	
    GaussCrd(2,0) = 1.0/6.0; GaussCrd(2,1) = 1.0/6.0; GaussCrd(2,2) = 2.0/3.0; 
    
    double wgauss = 0.33333333333333333333333333333333333333333;    
    for(unsigned int gp = 0; gp<3; gp++)
    {
      NN[0] = GaussCrd(gp,0); NN[1] = GaussCrd(gp,1); NN[2] = GaussCrd(gp,2);
      
      double HTC_Alpha = NN[0] * GetGeometry()[0].FastGetSolutionStepValue(rTransferCoefficientVar);
      HTC_Alpha += NN[1] * GetGeometry()[1].FastGetSolutionStepValue(rTransferCoefficientVar);
      HTC_Alpha += NN[2] * GetGeometry()[2].FastGetSolutionStepValue(rTransferCoefficientVar);      
      
      for(unsigned int ii=0; ii<nodes_number; ii++)
	for(unsigned int jj=0; jj<nodes_number; jj++)
	  rLeftHandSideMatrix(ii,jj) += HTC_Alpha * NN[ii] * NN[jj];      
    }
    
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    array_1d<double, 3 > unknown_vec; 
    for(unsigned int kk = 0; kk<nodes_number; kk++)
      unknown_vec[kk] = GetGeometry()[kk].FastGetSolutionStepValue(rUnknownVar) - amb_T;

     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, unknown_vec);
    
    rLeftHandSideMatrix *=  (area * wgauss);
    rRightHandSideVector *=  (area * wgauss);    

}

//************************************************************************************
//************************************************************************************

void EnvironmentContact3D::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

// 	int nodes_number = 4;
// 	int dim = 2;
// 	unsigned int matsize = nodes_number*(dim);
//
// 	if(rDampingMatrix.size1() != matsize)
// 			rDampingMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!




    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void EnvironmentContact3D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
{
    KRATOS_TRY


    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************
void EnvironmentContact3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void EnvironmentContact3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    }
    KRATOS_CATCH("");
}




} // Namespace Kratos



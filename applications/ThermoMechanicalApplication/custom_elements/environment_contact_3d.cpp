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
//#include "math.h"

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
    //stage = 0 is a pure convection step - this condition does nothing
    //stage = 1 is a diffusion step
    const unsigned int stage = rCurrentProcessInfo[FRACTIONAL_STEP];
    
    if(stage == 0) //pure convection step
    {
            if(rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0,0,false);
    
        if(rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0,false);

    }
    else if(stage == 1) //diffusion step
    {
        unsigned int nodes_number = GetGeometry().PointsNumber();

        if(rLeftHandSideMatrix.size1() != nodes_number)
            rLeftHandSideMatrix.resize(nodes_number,nodes_number,false);

        if(rRightHandSideVector.size() != nodes_number)
            rRightHandSideVector.resize(nodes_number,false);

        noalias(rRightHandSideVector) = ZeroVector(nodes_number);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nodes_number, nodes_number);

        //calculate area
        const double area = GetGeometry().Area();
        const double nodal_area = area/3.0;

        //take thermal properties
  //      ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
//        const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();
//        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

//        const double ambient_T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];

        //loop over integration poitns. Nodeal integration is assumed!!
        for(unsigned int igauss = 0; igauss<3; igauss++)
        {
               double HTC_Alpha, node_heat_flux;
               ComputeGaussHeatFluxAndHTC( igauss, HTC_Alpha, node_heat_flux, rCurrentProcessInfo);
               
               rLeftHandSideMatrix(igauss,igauss) += HTC_Alpha * nodal_area;  
               rRightHandSideVector[igauss] = node_heat_flux *  nodal_area;
               
//             const double HTC_Alpha = GetGeometry()[igauss].FastGetSolutionStepValue(rTransferCoefficientVar);
//             rLeftHandSideMatrix(igauss,igauss) += HTC_Alpha * nodal_area;  
//             rRightHandSideVector[igauss] = HTC_Alpha *  nodal_area * ( ambient_T - GetGeometry()[igauss].FastGetSolutionStepValue(rUnknownVar) );

        }
    }
}

void EnvironmentContact3D::ComputeGaussHeatFluxAndHTC( unsigned int igauss, 
                                             double& HTC_Alpha, 
                                             double& heat_flux, 
                                             ProcessInfo& rCurrentProcessInfo,
                                             bool save_internal_variables)
{
    const double ambient_T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        
    HTC_Alpha = GetGeometry()[igauss].FastGetSolutionStepValue(rTransferCoefficientVar);
    heat_flux = HTC_Alpha  * ( ambient_T - GetGeometry()[igauss].FastGetSolutionStepValue(rUnknownVar) );
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
        
        //stage = 0 is a pure convection step - this condition does nothing
    //stage = 1 is a diffusion step
    const unsigned int stage = CurrentProcessInfo[FRACTIONAL_STEP];
    
    if(stage == 0) //pure convection step
    {
        if(rResult.size() != 0)
            rResult.resize(0,false);
    }
    else if(stage == 1) //diffusion step
    {
        if(rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes,false);
        for (unsigned int i=0; i<number_of_nodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }
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
        
        //stage = 0 is a pure convection step - this condition does nothing
    //stage = 1 is a diffusion step
    const unsigned int stage = CurrentProcessInfo[FRACTIONAL_STEP];
    
    if(stage == 0) //pure convection step
    {
        if(ElementalDofList.size() != 0)
            ElementalDofList.resize(0);
    }
    else if(stage == 1) //diffusion step
    {
        if(ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        for (unsigned int i=0; i<number_of_nodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

        }
    }
    KRATOS_CATCH("");
}
//************************************************************************************
//************************************************************************************
// double EnvironmentContact3D::CalcTempSemiInfiniteWall(ProcessInfo& rCurrentProcessInfo, const double T_0)
// {
// 	// This function comutes the new mold temperature T(0,t) from the 
// 	// formulation presented in "BASIC HEAT and MASS TRANSFER" second edition page 170 for 
// 	// Convection heat transfer to the surface. In this formulation position, X, is put zero to 
// 	// obtain temperature on the wall at each time and ALPHA is compuated as k/ro*C
// 	//ERFC function is implemented as explained in Appendix B4.
//     ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
// 
//     const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
//     const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
//     const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();
// 	const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
// 
//     double wgauss = 0.33333333333333333333333333333333333333333;    
// 
// 	double ro =  GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
//  	double KK =  GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
// 	double CC =  GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
// 	double H = GetGeometry()[0].FastGetSolutionStepValue(rTransferCoefficientVar);
// 	double T_e = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
// 
// 	for (unsigned int i=1; i<3; i++)
//       {
// 		  ro += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
//           KK +=  GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
// 		  CC +=  GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
// 		  H += GetGeometry()[i].FastGetSolutionStepValue(rTransferCoefficientVar);
// 		  T_e += GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
// 	  }
// 
// 	ro *= wgauss;
// 	KK *= wgauss;
// 	CC *= wgauss;
// 	H *= wgauss;
// 	T_e *= wgauss;
// 
// 	double time = rCurrentProcessInfo[TIME];
// 	double aa = H/KK;
// 	double alpha = KK/(ro*CC);
// 
// 	double AA = exp(aa * aa * alpha * time);
// 
// 	double bb = aa * sqrt(alpha * time);
// 	double BB = ERFC(bb);
// 
// 	//T_e = rCurrentProcessInfo[FLUID_TEMPERATURE];
// 
// 	//(T-T_0)/(T_e-T_0) = 1 - AA * BB
// 	double new_T =  T_0;
// 	if( (AA * BB) >= 0.0 && (AA * BB) < 1.0)
// 	  new_T =  T_0 + (1.0 - AA * BB) * (T_e - T_0);
// 
// 	return new_T;
// 
// }

//************************************************************************************
//************************************************************************************
// double EnvironmentContact3D::ERFC( const double etta)
// {
// 
// 	double a_1 = 0.3480242;
// 	double a_2 = -0.0958798;
// 	double a_3 = 0.7478556;
// 	double p = 0.47047;
// 
// 	double t = 1.0/(1.0 + p*etta);
// 
// 	double erfc = (a_1 * t + a_2 * t*t + a_3 * t*t*t) * exp(-etta*etta);
// 
// 	return erfc;
// 
// }


} // Namespace Kratos



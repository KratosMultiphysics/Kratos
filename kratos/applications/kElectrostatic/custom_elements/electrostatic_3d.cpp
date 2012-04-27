/*
==============================================================================
KratosR1ElectrostaticApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
Pooyan Dadvand, Riccardo Rossi, Javier Mora
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mora@cimne.upc.edu
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
//   Last modified by:    $Author: rrossi  jmora $
//   Date:                $Date: 2010-02-02 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/electrostatic_3d.h"
#include "kElectrostatic.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Electrostatic3Dauxiliaries
{
//static variables
//boost::numeric::ublas::bounded_matrix<double,4,3> Electrostatic3D::msDN_DX;
boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3);
#pragma omp threadprivate(msDN_DX)
//boost::numeric::ublas::bounded_matrix<double,3,4> Electrostatic3D::msB;
boost::numeric::ublas::bounded_matrix<double,3,4> msB = ZeroMatrix(3,4);
#pragma omp threadprivate(msB)
//boost::numeric::ublas::bounded_matrix<double,3,3> Electrostatic3D::msD;
boost::numeric::ublas::bounded_matrix<double,3,3> msD = ZeroMatrix(3,3);
#pragma omp threadprivate(msD)

//array_1d<double,4> Electrostatic3D::msN; //dimension = number of nodes
array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
#pragma omp threadprivate(msN)
//array_1d<double,4> Electrostatic3D::ms_temp; //dimension = number of nodes
array_1d<double,4> ms_temp = ZeroVector(4); //dimension = number of nodes
#pragma omp threadprivate(ms_temp)
//array_1d<double,4> Electrostatic3D::point_sources; //dimension = number of nodes
array_1d<double,4> point_sources = ZeroVector(4); //dimension = number of nodes
#pragma omp threadprivate(point_sources)

array_1d<double,4> volume_sources = ZeroVector(4); //dimension = number of nodes
#pragma omp threadprivate(volume_sources)
}
using  namespace Electrostatic3Dauxiliaries;


//************************************************************************************
//************************************************************************************
Electrostatic3D::Electrostatic3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
Electrostatic3D::Electrostatic3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{

}

Element::Pointer Electrostatic3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new Electrostatic3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Electrostatic3D::~Electrostatic3D()
{
}

//************************************************************************************
//************************************************************************************
void Electrostatic3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_points = GetGeometry().size();

    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

    //reading properties and conditions
    array_1d<double,3> permittivity = GetProperties()[ELECTRICAL_PERMITTIVITY];
    msD(0,0)=permittivity[0];
    msD(1,1)=permittivity[1];
    msD(2,2)=permittivity[2];
    msD(0,1)=0.0;
    msD(0,2)=0.0;
    msD(1,0)=0.0;
    msD(1,2)=0.0;
    msD(2,0)=0.0;
    msD(2,1)=0.0;

    //double volume_sources = (this)->GetValue(ELECTROSTATIC_SURFACE_CHARGE);
    volume_sources[0] = (this)->GetValue(ELECTROSTATIC_VOLUME_CHARGE);
    volume_sources[1] = (this)->GetValue(ELECTROSTATIC_VOLUME_CHARGE);
    volume_sources[2] = (this)->GetValue(ELECTROSTATIC_VOLUME_CHARGE);
    volume_sources[3] = (this)->GetValue(ELECTROSTATIC_VOLUME_CHARGE);

    // main loop
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));

    /*		for(unsigned int k = 1; k<integration_points.size(); k++)	//integration points
    		{
    			double w_detj = integration_points[k].Weight()*mDetJo[k];
    */
    rLeftHandSideMatrix *= Area;
//		}


    noalias(rRightHandSideVector) = volume_sources*Area/4.0;
    //subtracting the dirichlet term
    // RHS -= LHS*ELECTROSTATIC_POTENTIALs
    for(unsigned int iii = 0; iii<number_of_points; iii++)
        ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);

    //multiplying by area, rho and density
    //rRightHandSideVector *= (Area * permittivity);
    //rLeftHandSideMatrix *= (Area * permittivity);
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void Electrostatic3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
}

//************************************************************************************
//************************************************************************************
// this subroutine calculates the nodal contributions for the explicit steps of the
// fractional step procedure
void Electrostatic3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}


//************************************************************************************
//************************************************************************************
void Electrostatic3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(ELECTROSTATIC_POTENTIAL).EquationId();
}

//************************************************************************************
//************************************************************************************
void Electrostatic3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
        ElementalDofList[i] = GetGeometry()[i].pGetDof(ELECTROSTATIC_POTENTIAL);

}

//************************************************************************************
//************************************************************************************
void Electrostatic3D::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
{
    IntegrationMethod mThisIntegrationMethod;

    mThisIntegrationMethod= GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
    Vector vect_tmp(4);
    vect_tmp[0] = GetGeometry()[0].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);
    vect_tmp[1] = GetGeometry()[1].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);
    vect_tmp[2] = GetGeometry()[2].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);
    vect_tmp[3] = GetGeometry()[3].FastGetSolutionStepValue(ELECTROSTATIC_POTENTIAL);

    //reading properties and conditions
    array_1d<double,3> permittivity = GetProperties()[ELECTRICAL_PERMITTIVITY];
    msD(0,0)=permittivity[0];
    msD(1,1)=permittivity[1];
    msD(2,2)=permittivity[2];
    msD(0,1)=0.0;
    msD(0,2)=0.0;
    msD(1,0)=0.0;
    msD(1,2)=0.0;
    msD(2,0)=0.0;
    msD(2,1)=0.0;

    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

    if(Output.size() != GetGeometry().IntegrationPoints(mThisIntegrationMethod).size())
        Output.resize(GetGeometry().IntegrationPoints(mThisIntegrationMethod).size());

    for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
    {
        if(rVariable==ELECTRIC_FIELD)
        {
            //KRATOS_WATCH("GiD Post Electrostatic - calc - elec - field");
            noalias(Output[PointNumber])=-prod(trans(msDN_DX),vect_tmp);
            //Output[PointNumber][2]=0.0;
//				KRATOS_WATCH(Output[PointNumber]);
        }
        else if(rVariable==ELECTRIC_DISPLACEMENT_FIELD)
        {
            //KRATOS_WATCH(rValues[PointNumber]);
            noalias(Output[PointNumber])=-prod(prod(msD,trans(msDN_DX)),vect_tmp);
//				KRATOS_WATCH(Output[PointNumber]);
        }
        else
        {
//				KRATOS_WATCH("GiD Post Electrostatic - calc - else");
            Output[PointNumber][0]=msD(0,0);
            Output[PointNumber][1]=msD(1,1);
            Output[PointNumber][2]=msD(2,2);
//				KRATOS_WATCH(Output[PointNumber]);
        }
    }



}


//************************************************************************************
//************************************************************************************


void Electrostatic3D::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    //KRATOS_WATCH("GiD Post Electrostatic - GetValueOnIntegrationPoints");

    if(rVariable==ELECTRIC_FIELD)
    {
//			KRATOS_WATCH("GiD Post Electrostatic - get - elec-field");
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else if(rVariable==ELECTRIC_DISPLACEMENT_FIELD)
    {
//			KRATOS_WATCH("GiD Post Electrostatic - get - elec-disp");
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
//			KRATOS_WATCH("GiD Post Electrostatic - get - else");
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
}


} // Namespace Kratos



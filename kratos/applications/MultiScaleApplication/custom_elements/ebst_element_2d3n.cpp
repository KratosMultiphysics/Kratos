/*
==============================================================================
KratosStructuralApplication
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

/* **************************************************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2008-10-13 07:00:53 $
 *   Revision:            $Revision: 1.12 $
 *
 * ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/ebst_element_2d3n.h"
#include "includes/constitutive_law.h"
#include "multiscale_application_variables.h"

#define EBST_ELEM_NUM_IP 1
#if EBST_ELEM_NUM_IP==1
#define EBST_SIMULATE_2D  
#endif // EBST_ELEM_NUM_IP==1


namespace Kratos
{

inline void GetRotationMatrixForStrains(double radians, Matrix & T)
{
	double c = std::cos(radians);
	double s = std::sin(radians);
	if(T.size1() != 3 || T.size2() != 3)
		T.resize(3, 3, false);
	noalias( T ) = ZeroMatrix(3, 3);
	T(0, 0) = c * c;			T(0, 1) =   s * s;				T(0, 2) = - s * c;
	T(1, 0) = s * s;			T(1, 1) =   c * c;				T(1, 2) =   s * c;
	T(2, 0) = 2.0 * s * c;		T(2, 1) = - 2.0 * s * c;		T(2, 2) = c * c - s * s;
}
inline void GetRotationMatrixForStresses(double radians, Matrix & T)
{
	double c = std::cos(radians);
	double s = std::sin(radians);
	if(T.size1() != 3 || T.size2() != 3)
		T.resize(3, 3, false);
	noalias( T ) = ZeroMatrix(3, 3);
	T(0, 0) = c * c;		T(0, 1) =   s * s;		T(0, 2) = - 2.0 * s * c;
	T(1, 0) = s * s;		T(1, 1) =   c * c;		T(1, 2) =   2.0 * s * c;
	T(2, 0) = s * c;		T(2, 1) = - s * c;		T(2, 2) = c * c - s * s;
}


// Constructor

EBSTElement2D3N::EBSTElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
	, mInitialized(false)
{
}

// Constructor

EBSTElement2D3N::EBSTElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
	, mInitialized(false)
{
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer EBSTElement2D3N::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new EBSTElement2D3N(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

EBSTElement2D3N::~EBSTElement2D3N()
{
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);

    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);
    unsigned int dim = number_of_nodes * 3;

    if (rResult.size() != dim)
        rResult.resize(dim, false);

    //nodes of the central element
    for (int i = 0; i < 3; i++)
    {
        int index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

    //adding the ids ofthe neighbouring nodes
    int index = 9;
    for (unsigned int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            rResult[index] = neigb[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = neigb[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = neigb[i].GetDof(DISPLACEMENT_Z).EquationId();
            index += 3;
        }
    }


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void EBSTElement2D3N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
    ElementalDofList.resize(0);

    //nodes of the central element
    for (unsigned int i = 0; i < GetGeometry().size(); i++)
    {
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    //adding the dofs ofthe neighbouring nodes
    for (unsigned int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_Z));
        }
    }
    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::GetValuesVector(
    Vector& values,
    int Step)
{
    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize);

    //nodes of the central element
    for (unsigned int i = 0; i < 3; i++)
    {
        const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * 3;
        values[index] = disp[0];
        values[index + 1] = disp[1];
        values[index + 2] = disp[2];
    }

    //neighbour nodes
    int index = 9;
    for (int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            const array_1d<double, 3 > & disp = neigb[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            values[index] = disp[0];
            values[index + 1] = disp[1];
            values[index + 2] = disp[2];
            index += 3;
        }
    }
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::GetFirstDerivativesVector(
    Vector& values,
    int Step)
{
    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize);

    //nodes of the central element
    for (unsigned int i = 0; i < 3; i++)
    {
        const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        unsigned int index = i * 3;
        values[index] = vel[0];
        values[index + 1] = vel[1];
        values[index + 2] = vel[2];
    }

    //neighbour nodes
    int index = 9;
    for (int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            const array_1d<double, 3 > & vel = neigb[i].FastGetSolutionStepValue(VELOCITY, Step);
            values[index] = vel[0];
            values[index + 1] = vel[1];
            values[index + 2] = vel[2];
            index += 3;
        }
    }
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::GetSecondDerivativesVector(
    Vector& values,
    int Step)
{
    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
    unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

    const unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize);

    //nodes of the central element
    for (unsigned int i = 0; i < 3; i++)
    {
        const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        unsigned int index = i * 3;
        values[index] = acc[0];
        values[index + 1] = acc[1];
        values[index + 2] = acc[2];
    }

    //neighbour nodes
    int index = 9;
    for (int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            const array_1d<double, 3 > & acc = neigb[i].FastGetSolutionStepValue(ACCELERATION, Step);
            values[index] = acc[0];
            values[index + 1] = acc[1];
            values[index + 2] = acc[2];
            index += 3;
        }
    }
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

// KRATOS_WATCH(rRightHandSideVector);


}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

//        KRATOS_WATCH(mstrains);
//        KRATOS_WATCH(Id());
//                            KRATOS_WATCH(mstrains);
//        Vector prova1(3);
//        prova1[0] = 10e-3;
//        prova1[1] = 5e-3;
//        prova1[2] = 5e-3;
//
//        Vector stress1(3);
//        Matrix tangent1(3,3);
//        mConstitutiveLawVector[0]->CalculateMaterialResponse(prova1,stress1,tangent1,true,true,false);
//        KRATOS_WATCH(prova1);
//        KRATOS_WATCH(stress1);
//        Vector dp(3);
//        dp[0] = 1e-4;
//        dp[1] = 1e-4;
//        dp[2] = 1e-4;
//
//        Vector prova2(3);
//        Matrix tangent2(3,3);
//        Vector stress2(3);
//        prova2 = prova1 + dp;
//        mConstitutiveLawVector[0]->CalculateMaterialResponse(prova2,stress2,tangent2,true,true,false);
//        KRATOS_WATCH(prova2);
//        KRATOS_WATCH(stress2);
//
//        Vector by_tangent = stress1;
//        noalias(by_tangent) += prod(tangent1,dp);
//        KRATOS_WATCH(by_tangent);

}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);

    //rMassMatrix.resize(0,0);
    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
    unsigned int MatSize = number_of_nodes * 3;
    if (rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize, MatSize, false);
    noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);

    double Area = GetGeometry().Area();
    double TotalMass = Area * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
    Vector LumpFact;
    GetGeometry().LumpingFactors(LumpFact);
// KRATOS_WATCH(GetProperties()[DENSITY]);
    for (unsigned int i = 0; i < 3; i++)
    {
        double temp = LumpFact[i] * TotalMass;
        for (unsigned int j = 0; j < 3; j++)
        {
            unsigned int index = i * 3 + j;
            rMassMatrix(index, index) = temp;
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);

    //rMassMatrix.resize(0,0);
    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
    unsigned int MatSize = number_of_nodes * 3;
    if (rDampingMatrix.size1() != MatSize)
        rDampingMatrix.resize(MatSize, MatSize, false);

    noalias(rDampingMatrix) = ZeroMatrix(MatSize, MatSize);
//         if (rDampingMatrix.size1() != 0)
//             rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void EBSTElement2D3N::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo)
{
    Vector N(3);
    N[0] = 0.33333333333333333333;
    N[1] = 0.33333333333333333333;
    N[2] = 0.33333333333333333333;

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        //call the finalize solution step function
        mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),GetGeometry(),N, rCurrentProcessInfo);
    }
}






void EBSTElement2D3N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
																	 std::vector<double>& rOutput,
																	 const ProcessInfo& rCurrentProcessInfo)
{
	GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
}

void EBSTElement2D3N::CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
																	 std::vector< Vector >& rOutput,
																	 const ProcessInfo& rCurrentProcessInfo)
{
	GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
}

void EBSTElement2D3N::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
																	 std::vector< Matrix >& rOutput,
																	 const ProcessInfo& rCurrentProcessInfo)
{
	GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
}

void EBSTElement2D3N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                       std::vector<double>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() != 1) rValues.resize(1);
	rValues[0] = 0.0;
	mConstitutiveLawVector[0]->GetValue(rVariable, rValues[0]);
}

void EBSTElement2D3N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                       std::vector<Vector>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() != 1) rValues.resize(1);
	rValues[0] = ZeroVector(1);
	mConstitutiveLawVector[0]->GetValue(rVariable, rValues[0]);
}

void EBSTElement2D3N::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                       std::vector<Matrix>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() != 1) rValues.resize(1);

	if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||
		   rVariable == PK2_STRESS_TENSOR ||
		   rVariable == CAUCHY_STRESS_TENSOR)
	{
		boost::numeric::ublas::bounded_matrix<double, 6, 18 > msL1;
		boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_f;
		boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_m;
		boost::numeric::ublas::bounded_matrix<double, 18, 18 > msK;
		boost::numeric::ublas::bounded_matrix<double, 6, 3 > ms_coord;

		WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
	
		//fill the aux matrix of coordinates
		for (unsigned int i = 0; i < 3; i++)
		{
			ms_coord(i, 0) = GetGeometry()[i].X();
			ms_coord(i, 1) = GetGeometry()[i].Y();
			ms_coord(i, 2) = GetGeometry()[i].Z();
		}
		for (unsigned int i = 0; i < 3; i++)
		{
			ms_coord(i + 3, 0) = neigb[i].X();
			ms_coord(i + 3, 1) = neigb[i].Y();
			ms_coord(i + 3, 2) = neigb[i].Z();
		}
	
		//compute phis on center and all gauss points
		boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiM;
		boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG1;
		boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG2;
		boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG3;
		CalculatePhiM(phiM, mdcgM, ms_coord);
		CalculatePhiG(phiG1, mdcg1, ms_coord);
		CalculatePhiG(phiG2, mdcg2, ms_coord);
		CalculatePhiG(phiG3, mdcg3, ms_coord);

		//calculating the normal to the element
		array_1d<double, 3 > phi1, phi2, t3e;
		phi1[0] = phiM(0, 0);
		phi1[1] = phiM(0, 1);
		phi1[2] = phiM(0, 2);
		phi2[0] = phiM(1, 0);
		phi2[1] = phiM(1, 1);
		phi2[2] = phiM(1, 2);
		MathUtils<double>::CrossProduct(t3e, phi1, phi2);
		double nze = norm_2(t3e);
		t3e /= nze;

		//calculating derivative of h --> see Eqn 62
		unsigned int ind = 0;
		for (unsigned int i = 0; i < 6; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				msL1(0, ind) = mdcg1(0, i) * t3e[j];
				msL1(1, ind) = mdcg1(1, i) * t3e[j];

				msL1(2, ind) = mdcg2(0, i) * t3e[j];
				msL1(3, ind) = mdcg2(1, i) * t3e[j];

				msL1(4, ind) = mdcg3(0, i) * t3e[j];
				msL1(5, ind) = mdcg3(1, i) * t3e[j];
				ind++;
			}
		}

		boost::numeric::ublas::bounded_matrix<double, 3, 6 > DN = ZeroMatrix(3, 6);
		DN(0, 0) = mdcgM(0, 0);
		DN(1, 1) = mdcgM(1, 0);
		DN(2, 0) = mdcgM(1, 0);
		DN(2, 1) = mdcgM(0, 0);
		DN(0, 2) = mdcgM(0, 1);
		DN(1, 3) = mdcgM(1, 1);
		DN(2, 2) = mdcgM(1, 1);
		DN(2, 3) = mdcgM(0, 1);
		DN(0, 4) = mdcgM(0, 2);
		DN(1, 5) = mdcgM(1, 2);
		DN(2, 4) = mdcgM(1, 2);
		DN(2, 5) = mdcgM(0, 2);

		noalias(msB_f) = 2.0 * prod(DN, msL1);
	
		//take care ... in msB_f we miss the second part of Eqn 62!
		array_1d<double, 3 > h00, h11, h01;
		Calculate_h_ab(h00, 0, 0, phiG1, phiG2, phiG3, mdcgM);
		Calculate_h_ab(h11, 1, 1, phiG1, phiG2, phiG3, mdcgM);
		Calculate_h_ab(h01, 0, 1, phiG1, phiG2, phiG3, mdcgM);
		array_1d<double, 3 > fhicm1, fhicm2;
		fhicm1[0] = phi2[1] * t3e[2] - phi2[2] * t3e[1];
		fhicm1[1] = phi2[2] * t3e[0] - phi2[0] * t3e[2];
		fhicm1[2] = phi2[0] * t3e[1] - phi2[1] * t3e[0];

		fhicm2[0] = -phi1[1] * t3e[2] + phi1[2] * t3e[1];
		fhicm2[1] = -phi1[2] * t3e[0] + phi1[0] * t3e[2];
		fhicm2[2] = -phi1[0] * t3e[1] + phi1[1] * t3e[0];

		double ro111 = h00[0] * fhicm1[0] + h00[1] * fhicm1[1] + h00[2] * fhicm1[2];
		double ro122 = h11[0] * fhicm1[0] + h11[1] * fhicm1[1] + h11[2] * fhicm1[2];
		double ro211 = h00[0] * fhicm2[0] + h00[1] * fhicm2[1] + h00[2] * fhicm2[2];
		double ro222 = h11[0] * fhicm2[0] + h11[1] * fhicm2[1] + h11[2] * fhicm2[2];
		double ro112 = h01[0] * fhicm1[0] + h01[1] * fhicm1[1] + h01[2] * fhicm1[2];
		double ro212 = h01[0] * fhicm2[0] + h01[1] * fhicm2[1] + h01[2] * fhicm2[2];

		const double& N11 = mdcgM(0, 0);
		const double& N12 = mdcgM(1, 0);
		const double& N21 = mdcgM(0, 1);
		const double& N22 = mdcgM(1, 1);
		const double& N31 = mdcgM(0, 2);
		const double& N32 = mdcgM(1, 2);
	
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > m;
		m(0, 0) = 2 * (N11 * ro111 + N12 * ro211);
		m(1, 0) = 2 * (N11 * ro122 + N12 * ro222);
		m(2, 0) = 2 * (N11 * ro112 + N12 * ro212);

		m(0, 1) = 2 * (N21 * ro111 + N22 * ro211);
		m(1, 1) = 2 * (N21 * ro122 + N22 * ro222);
		m(2, 1) = 2 * (N21 * ro112 + N22 * ro212);

		m(0, 2) = 2 * (N31 * ro111 + N32 * ro211);
		m(1, 2) = 2 * (N31 * ro122 + N32 * ro222);
		m(2, 2) = 2 * (N31 * ro112 + N32 * ro212);

		boost::numeric::ublas::bounded_matrix<double, 3, 18 > n;
		noalias(n) = ZeroMatrix(3, 18);
		n(0, 0) = t3e[0];
		n(0, 1) = t3e[1];
		n(0, 2) = t3e[2];
		n(1, 3) = t3e[0];
		n(1, 4) = t3e[1];
		n(1, 5) = t3e[2];
		n(2, 6) = t3e[0];
		n(2, 7) = t3e[1];
		n(2, 8) = t3e[2];

		//calculate bending strain
		array_1d<double, 3 > bending_strain;
		bending_strain[0] = inner_prod(h00, t3e);
		bending_strain[1] = inner_prod(h11, t3e);
		bending_strain[2] = inner_prod(h01, t3e);

		//if the original curvature was not yet initialized
		if (mK0[0] == 1234567.89)
		{
			noalias(mK0) = bending_strain;
		}
	
		//subtract to the curvature the original curvature
		noalias(bending_strain) -= mK0;

		//calculating here the membrane stiffness
		//strain is computed on all of the gauss points and then averaged following
		//the assumed strain approach
		array_1d<double, 3 > membrane_strain;
		noalias(membrane_strain) = ZeroVector(3);
		CalculateAndAdd_MembraneStrain(membrane_strain, phiG1);
		CalculateAndAdd_MembraneStrain(membrane_strain, phiG2);
		CalculateAndAdd_MembraneStrain(membrane_strain, phiG3);
		membrane_strain *= 0.333333333333333333333333333333333;

		if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
		{

			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// calculate the angle between material and local element axis

			// calculate the angle between the element x direction and the material x direction.
			const array_1d<double, 3 >& a = phi1;
			array_1d<double, 3 > b; // material direction. here in 2d [1,0,0]
			b[0]=1.0; b[1]=0.0; b[2]=0.0;
			double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
			if(a_dot_b < -1.0) a_dot_b = -1.0;
			if(a_dot_b >  1.0) a_dot_b =  1.0;
			double angle = std::acos( a_dot_b );
			// if they are not counter-clock-wise, let's change the sign of the angle
			if(angle != 0.0) {
				array_1d<double, 3 > elem_vz; elem_vz[0]=0.0; elem_vz[1]=0.0; elem_vz[2]=1.0;
				array_1d<double, 3> elem_vy;
				MathUtils<double>::CrossProduct(elem_vy, elem_vz, phi1);
				if( b(0)*elem_vy(0) + b(1)*elem_vy(1) + b(2)*elem_vy(2) < 0.0 )
					angle = -angle;
			}
			Matrix RE;
			GetRotationMatrixForStrains(-angle,RE);

			// rotate the elemental strains to the material coordinate system
			membrane_strain = prod(RE, membrane_strain);

			// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			array_1d<double,6> temp;
			array_1d<double,6> global_strain;
			array_1d<double,3>& v1 = phi1;
			array_1d<double,3>& v2 = phi2;

			//adding the component S11
			noalias(temp)  = VoigtTensorComponents(v1,v1);
			temp *= membrane_strain[0];
			noalias(global_strain) = temp;
			//adding the component S22
			noalias(temp)  = VoigtTensorComponents(v2,v2);
			temp *= membrane_strain[1];
			noalias(global_strain) += temp;
			//adding the component S12 (& S21)
			noalias(temp)  = VoigtTensorComponents(v1,v2);
			noalias(temp) += VoigtTensorComponents(v2,v1);
			temp *= membrane_strain[2];
			noalias(global_strain) += temp;

			global_strain.clear();
			global_strain(0) = membrane_strain(0);
			global_strain(1) = membrane_strain(1);
			global_strain(3) = membrane_strain(2);

			rValues[0] = Matrix(3,3);
			noalias(rValues[0]) = MathUtils<double>::StrainVectorToTensor(global_strain);
		}
		else
		{
			array_1d<double, 3 >& membrane_stress = m_membrane_stress;
			array_1d<double, 3 >& bending_stress = m_bending_stress;
			boost::numeric::ublas::bounded_matrix<double, 3, 3 > Dmat_m;
			boost::numeric::ublas::bounded_matrix<double, 3, 3 > Dmat_f;
			double h_on_h0;
			CalculateEquivalentStresses(
				membrane_strain, bending_strain, 
				Dmat_m, Dmat_f, 
				membrane_stress, bending_stress, 
				h_on_h0,rCurrentProcessInfo,phi1);
			
			array_1d<double,6> temp;
			array_1d<double,6> global_stress;
			array_1d<double,3>& v1 = phi1;
			array_1d<double,3>& v2 = phi2;

			//adding the component S11
			noalias(temp)  = VoigtTensorComponents(v1,v1);
			temp *= membrane_stress[0];
			noalias(global_stress) = temp;
			//adding the component S22
			noalias(temp)  = VoigtTensorComponents(v2,v2);
			temp *= membrane_stress[1];
			noalias(global_stress) += temp;
			//adding the component S12 (& S21)
			noalias(temp)  = VoigtTensorComponents(v1,v2);
			noalias(temp) += VoigtTensorComponents(v2,v1);
			temp *= membrane_stress[2];
			noalias(global_stress) += temp;
			global_stress /=  GetGeometry().Area() * GetProperties()[THICKNESS];

			rValues[0] = Matrix(3,3);
			noalias(rValues[0]) = MathUtils<double>::StressVectorToTensor(global_stress);
		}
	}
	else
	{
		rValues[0] = ZeroMatrix(1,1);
		mConstitutiveLawVector[0]->GetValue(rVariable, rValues[0]);
	}
}

void EBSTElement2D3N::SetValueOnIntegrationPoints(const Variable<double>& rVariable, 
												   std::vector<double>& rValues, 
												   const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() == mConstitutiveLawVector.size())
		for(unsigned int i = 0; i < rValues.size(); i++)
			mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void EBSTElement2D3N::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, 
												   std::vector<Vector>& rValues, 
												   const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() == mConstitutiveLawVector.size())
		for(unsigned int i = 0; i < rValues.size(); i++)
			mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void EBSTElement2D3N::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, 
												   std::vector<Matrix>& rValues, 
												   const ProcessInfo& rCurrentProcessInfo)
{
	if(rValues.size() == mConstitutiveLawVector.size())
		for(unsigned int i = 0; i < rValues.size(); i++)
			mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
}

void EBSTElement2D3N::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
												 std::vector<ConstitutiveLaw::Pointer>& rValues,
												 const ProcessInfo& rCurrentProcessInfo )
{
	if ( mConstitutiveLawVector.size() != rValues.size() )
	{
		mConstitutiveLawVector.resize(rValues.size());
		if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( ) )
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
	}
	for(unsigned int i=0; i<rValues.size(); i++)
		mConstitutiveLawVector[i] = rValues[i];
}






void EBSTElement2D3N::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag)
{
    boost::numeric::ublas::bounded_matrix<double, 6, 18 > msL1;
    boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_f;
    boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_m;
    boost::numeric::ublas::bounded_matrix<double, 18, 18 > msK;
    boost::numeric::ublas::bounded_matrix<double, 6, 3 > ms_coord;

    WeakPointerVector< Node < 3 > >& neigb = this->GetValue(NEIGHBOUR_NODES);
	
    //fill the aux matrix of coordinates
    for (unsigned int i = 0; i < 3; i++)
    {
        ms_coord(i, 0) = GetGeometry()[i].X();
        ms_coord(i, 1) = GetGeometry()[i].Y();
        ms_coord(i, 2) = GetGeometry()[i].Z();
    }
    for (unsigned int i = 0; i < 3; i++)
    {
        ms_coord(i + 3, 0) = neigb[i].X();
        ms_coord(i + 3, 1) = neigb[i].Y();
        ms_coord(i + 3, 2) = neigb[i].Z();
    }
	
    //compute phis on center and all gauss points
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiM;
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG1;
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG2;
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > phiG3;
    CalculatePhiM(phiM, mdcgM, ms_coord);
    CalculatePhiG(phiG1, mdcg1, ms_coord);
    CalculatePhiG(phiG2, mdcg2, ms_coord);
    CalculatePhiG(phiG3, mdcg3, ms_coord);

    //calculating the normal to the element
    array_1d<double, 3 > phi1, phi2, t3e;
    phi1[0] = phiM(0, 0);
    phi1[1] = phiM(0, 1);
    phi1[2] = phiM(0, 2);
    phi2[0] = phiM(1, 0);
    phi2[1] = phiM(1, 1);
    phi2[2] = phiM(1, 2);
    MathUtils<double>::CrossProduct(t3e, phi1, phi2);
    double nze = norm_2(t3e);
    t3e /= nze;

    //calculating derivative of h --> see Eqn 62
    unsigned int ind = 0;
    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            msL1(0, ind) = mdcg1(0, i) * t3e[j];
            msL1(1, ind) = mdcg1(1, i) * t3e[j];

            msL1(2, ind) = mdcg2(0, i) * t3e[j];
            msL1(3, ind) = mdcg2(1, i) * t3e[j];

            msL1(4, ind) = mdcg3(0, i) * t3e[j];
            msL1(5, ind) = mdcg3(1, i) * t3e[j];
            ind++;
        }
    }

    boost::numeric::ublas::bounded_matrix<double, 3, 6 > DN = ZeroMatrix(3, 6);
    DN(0, 0) = mdcgM(0, 0);
    DN(1, 1) = mdcgM(1, 0);
    DN(2, 0) = mdcgM(1, 0);
    DN(2, 1) = mdcgM(0, 0);
    DN(0, 2) = mdcgM(0, 1);
    DN(1, 3) = mdcgM(1, 1);
    DN(2, 2) = mdcgM(1, 1);
    DN(2, 3) = mdcgM(0, 1);
    DN(0, 4) = mdcgM(0, 2);
    DN(1, 5) = mdcgM(1, 2);
    DN(2, 4) = mdcgM(1, 2);
    DN(2, 5) = mdcgM(0, 2);

    noalias(msB_f) = 2.0 * prod(DN, msL1);
	
    //take care ... in msB_f we miss the second part of Eqn 62!
    array_1d<double, 3 > h00, h11, h01;
    Calculate_h_ab(h00, 0, 0, phiG1, phiG2, phiG3, mdcgM);
    Calculate_h_ab(h11, 1, 1, phiG1, phiG2, phiG3, mdcgM);
    Calculate_h_ab(h01, 0, 1, phiG1, phiG2, phiG3, mdcgM);
    array_1d<double, 3 > fhicm1, fhicm2;
    fhicm1[0] = phi2[1] * t3e[2] - phi2[2] * t3e[1];
    fhicm1[1] = phi2[2] * t3e[0] - phi2[0] * t3e[2];
    fhicm1[2] = phi2[0] * t3e[1] - phi2[1] * t3e[0];

    fhicm2[0] = -phi1[1] * t3e[2] + phi1[2] * t3e[1];
    fhicm2[1] = -phi1[2] * t3e[0] + phi1[0] * t3e[2];
    fhicm2[2] = -phi1[0] * t3e[1] + phi1[1] * t3e[0];

    double ro111 = h00[0] * fhicm1[0] + h00[1] * fhicm1[1] + h00[2] * fhicm1[2];
    double ro122 = h11[0] * fhicm1[0] + h11[1] * fhicm1[1] + h11[2] * fhicm1[2];
    double ro211 = h00[0] * fhicm2[0] + h00[1] * fhicm2[1] + h00[2] * fhicm2[2];
    double ro222 = h11[0] * fhicm2[0] + h11[1] * fhicm2[1] + h11[2] * fhicm2[2];
    double ro112 = h01[0] * fhicm1[0] + h01[1] * fhicm1[1] + h01[2] * fhicm1[2];
    double ro212 = h01[0] * fhicm2[0] + h01[1] * fhicm2[1] + h01[2] * fhicm2[2];

    const double& N11 = mdcgM(0, 0);
    const double& N12 = mdcgM(1, 0);
    const double& N21 = mdcgM(0, 1);
    const double& N22 = mdcgM(1, 1);
    const double& N31 = mdcgM(0, 2);
    const double& N32 = mdcgM(1, 2);
	
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > m;
    m(0, 0) = 2 * (N11 * ro111 + N12 * ro211);
    m(1, 0) = 2 * (N11 * ro122 + N12 * ro222);
    m(2, 0) = 2 * (N11 * ro112 + N12 * ro212);

    m(0, 1) = 2 * (N21 * ro111 + N22 * ro211);
    m(1, 1) = 2 * (N21 * ro122 + N22 * ro222);
    m(2, 1) = 2 * (N21 * ro112 + N22 * ro212);

    m(0, 2) = 2 * (N31 * ro111 + N32 * ro211);
    m(1, 2) = 2 * (N31 * ro122 + N32 * ro222);
    m(2, 2) = 2 * (N31 * ro112 + N32 * ro212);

    boost::numeric::ublas::bounded_matrix<double, 3, 18 > n;
    noalias(n) = ZeroMatrix(3, 18);
    n(0, 0) = t3e[0];
    n(0, 1) = t3e[1];
    n(0, 2) = t3e[2];
    n(1, 3) = t3e[0];
    n(1, 4) = t3e[1];
    n(1, 5) = t3e[2];
    n(2, 6) = t3e[0];
    n(2, 7) = t3e[1];
    n(2, 8) = t3e[2];

    //calculate bending strain
    array_1d<double, 3 > bending_strain;
    bending_strain[0] = inner_prod(h00, t3e);
    bending_strain[1] = inner_prod(h11, t3e);
    bending_strain[2] = inner_prod(h01, t3e);

    //if the original curvature was not yet initialized
    if (mK0[0] == 1234567.89)
    {
        noalias(mK0) = bending_strain;
    }
	
    //subtract to the curvature the original curvature
    noalias(bending_strain) -= mK0;

    //calculating here the membrane stiffness
    //strain is computed on all of the gauss points and then averaged following
    //the assumed strain approach
    array_1d<double, 3 > membrane_strain;
    noalias(membrane_strain) = ZeroVector(3);
    CalculateAndAdd_MembraneStrain(membrane_strain, phiG1);
    CalculateAndAdd_MembraneStrain(membrane_strain, phiG2);
    CalculateAndAdd_MembraneStrain(membrane_strain, phiG3);
    membrane_strain *= 0.333333333333333333333333333333333;

    //calculating the membrane strain-displacement matrix
    noalias(msB_m) = ZeroMatrix(3, 18);
    CalculateAndAdd_MembraneB(msB_m, mdcg1, phiG1);
    CalculateAndAdd_MembraneB(msB_m, mdcg2, phiG2);
    CalculateAndAdd_MembraneB(msB_m, mdcg3, phiG3);
    msB_m *= 0.333333333333333333333333333333333;

    array_1d<double, 3 >& membrane_stress = m_membrane_stress;
    array_1d<double, 3 >& bending_stress = m_bending_stress;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Dmat_m;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Dmat_f;
    double h_on_h0;
    CalculateEquivalentStresses(
		membrane_strain, bending_strain, 
		Dmat_m, Dmat_f, 
		membrane_stress, bending_stress, 
		h_on_h0,rCurrentProcessInfo,phi1);
	
    //bending contribution to the elemental stiffness
    boost::numeric::ublas::bounded_matrix<double, 3, 18 > aux;
    noalias(aux) = prod(Dmat_f, msB_f);
    noalias(msK) = prod(trans(msB_f), aux);

    //adding the membrane contribution to the elemental stiffness
    noalias(aux) = prod(Dmat_m, msB_m);
    noalias(msK) += prod(trans(msB_m), aux);

    //adding the geometric membrane stiffness
    CalculateAndAdd_Membrane_Kg(msK, mdcg1, mdcg2, mdcg3, membrane_stress);
	
    //calculate rhs
    array_1d<double, 18 > rhs_full = ZeroVector(18);
    /*const array_1d<double, 3 > & aaa = GetProperties()[BODY_FORCE];
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            rhs_full[i * 3 + j] = aaa[j] * Area0 * 0.3333333333333333333333;*/
	// body force with density and volume acceleration
	array_1d<double,3> body_force_vec;
	body_force_vec.clear();
	for(unsigned int i=0; i<3; i++) {
		if(GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
			body_force_vec += GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
		}
	}
	body_force_vec /= 3.0; // constant value at the centroid
	body_force_vec *= GetProperties()[DENSITY] * Area0 * GetProperties()[THICKNESS];
	for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            rhs_full[i * 3 + j] = body_force_vec(j)/3.0; // one third on each node

    noalias(rhs_full) -= prod(trans(msB_f), bending_stress);

    noalias(rhs_full) -= prod(trans(msB_m), membrane_stress);

    //assembling the elemental contribution to the correct size.
    unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
    unsigned int MatSize = number_of_nodes * 3;

    //resizing as needed the LHS
    if (rLeftHandSideMatrix.size1() != MatSize)
        rLeftHandSideMatrix.resize(MatSize, MatSize);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS

    if (rRightHandSideVector.size() != MatSize)
        rRightHandSideVector.resize(MatSize, false);

    array_1d<unsigned int, 18 > id_vec;
    for (unsigned int i = 0; i < 9; i++) id_vec[i] = i;
    unsigned int index = 9;
    for (unsigned int i = 0; i < 3; i++)
    {
        if (HasNeighbour(i, neigb[i]))
        {
            for (unsigned int j = 0; j < 3; j++)
                id_vec[9 + i * 3 + j] = index + j;
            index += 3;
        }
        else
        {
            for (unsigned int j = 0; j < 3; j++)
                id_vec[9 + i * 3 + j] = 1000;
        }
    }

	

    //add the first 9*9 block
    for (unsigned int i = 0; i < 18; i++)
    {
        if (id_vec[i] < 1000)
        {
            rRightHandSideVector(id_vec[i]) = rhs_full[i];
            for (unsigned int j = 0; j < 18; j++)
            {
                if (id_vec[j] < 1000) rLeftHandSideMatrix(id_vec[i], id_vec[j]) = msK(i, j);
            }
        }
    }
}

bool EBSTElement2D3N::HasNeighbour(unsigned int index, const Node < 3 > & neighb)
{
    if (neighb.Id() == GetGeometry()[index].Id())
        return false;
    else
        return true;
}

unsigned int EBSTElement2D3N::NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs)
{
    unsigned int active_neighbours = 0;
    for (unsigned int i = 0; i < neighbs.size(); i++)
        if (HasNeighbour(i, neighbs[i])) active_neighbours++;
    return active_neighbours;
}

void EBSTElement2D3N::Initialize()
{
	if(mInitialized) return;
	mInitialized = true;

    KRATOS_TRY

//         boost::numeric::ublas::bounded_matrix<double, 6, 18 > msL1;
//         boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_f;
//         boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_m;
//         boost::numeric::ublas::bounded_matrix<double, 18, 18 > msK;
    boost::numeric::ublas::bounded_matrix<double, 6, 3 > ms_coord;
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > ms_loc_der_central;
    boost::numeric::ublas::bounded_matrix<double, 6, 2 > ms_loc_der_patch;

    //find the "nodal neighbours" given the elemental neighbours
    WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
    if (elem_neigb.size() == 0) KRATOS_THROW_ERROR(std::logic_error, "the neighbour elements are not calculated", "")
        WeakPointerVector< Node < 3 > >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
    nodal_neigb.resize(3);
    Geometry< Node < 3 > >& center_geom = GetGeometry();
    /*KRATOS_WATCH(this->GetValue(NEIGHBOUR_ELEMENTS).size());
        std::cout << "I am elem" << Id() << std::endl;
        std::cout << "neighbours =" << elem_neigb[0].Id() << " " << elem_neigb[1].Id() << " " << elem_neigb[2].Id() << " " << std::endl;*/
    for (unsigned int i = 0; i < center_geom.size(); i++)
    {
        if (elem_neigb[i].Id() != Id()) //if the elemental neighbour exists
        {
            Geometry< Node < 3 > >& geom = elem_neigb[i].GetGeometry();
            for (unsigned int j = 0; j < geom.size(); j++)
            {
                bool aux = false;
                for (unsigned int k = 0; k < center_geom.size(); k++)
                {
                    if (geom[j].Id() == center_geom[k].Id())
                        aux = true;
                }

                if (aux == false) nodal_neigb(i) = Node < 3 > ::WeakPointer(geom(j));
            }
        }
        else   //the elemenetal neighbour does not exist
            nodal_neigb(i) = Node < 3 > ::WeakPointer(center_geom(i));
    }

    //    std::cout << "node1" << GetGeometry()[0].Id() << "opposite node =" << nodal_neigb[0].Id() << std::endl;
    //    std::cout << "node2" << GetGeometry()[1].Id() << "opposite node =" << nodal_neigb[1].Id() << std::endl;
    //    std::cout << "node3" << GetGeometry()[2].Id() << "opposite node =" << nodal_neigb[2].Id() << std::endl;

    boost::numeric::ublas::bounded_matrix<double, 2, 2 > ijac;


    //**********************************************************************************
    //fill the aux matrix of coordinates
    for (unsigned int i = 0; i < 3; i++)
    {
        ms_coord(i, 0) = GetGeometry()[i].X();
        ms_coord(i, 1) = GetGeometry()[i].Y();
        ms_coord(i, 2) = GetGeometry()[i].Z();
    }
    for (unsigned int i = 0; i < 3; i++)
    {
        ms_coord(i + 3, 0) = nodal_neigb[i].X();
        ms_coord(i + 3, 1) = nodal_neigb[i].Y();
        ms_coord(i + 3, 2) = nodal_neigb[i].Z();
    }

    //**********************************************************************************
    //calculate local system of coordinates of the elem ->SysCartE
    array_1d<double, 3 > v12, v13, vze;
    v12[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
    v12[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
    v12[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();
    double n12 = norm_2(v12);
    v12 /= n12;

    v13[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
    v13[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
    v13[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();
    double n13 = norm_2(v13);
    v13 /= n13;

    //vze = cross prod v12 v13 --> then normalize
    MathUtils<double>::CrossProduct(vze, v12, v13);
    double nze = norm_2(vze);
    vze /= nze;

    //version by Riccardo
    noalias(mvxe) = v12;

    MathUtils<double>::CrossProduct(mvye, vze, mvxe);

    //    //to compare with francisco
    //    mvxe[0] = 1.0;
    //    mvxe[1] = 0.0;
    //    mvxe[2] = 0.0;
    //
    //    mvye[0] = 0.0;
    //    mvye[1] = 1.0;
    //    mvye[2] = 0.0;
    //
    //    vze[0] = 0.0;
    //    vze[1] = 0.0;
    //    vze[2] = 1.0;

    //*****************************************************************************
    //calculate cartesian derivatives for the central element
    ms_loc_der_central(0, 0) = -1.0;
    ms_loc_der_central(0, 1) = -1.0;
    ms_loc_der_central(1, 0) = 1.0;
    ms_loc_der_central(1, 1) = 0.0;
    ms_loc_der_central(2, 0) = 0.0;
    ms_loc_der_central(2, 1) = 1.0;
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > phiM = ZeroMatrix(3, 2);
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            phiM(j, 0) += ms_loc_der_central(i, 0) * ms_coord(i, j);
            phiM(j, 1) += ms_loc_der_central(i, 1) * ms_coord(i, j);
        }
    }

    boost::numeric::ublas::bounded_matrix<double, 2, 2 > jac;
    jac(0, 0) = phiM(0, 0) * mvxe[0] + phiM(1, 0) * mvxe[1] + phiM(2, 0) * mvxe[2];
    jac(1, 0) = phiM(0, 1) * mvxe[0] + phiM(1, 1) * mvxe[1] + phiM(2, 1) * mvxe[2];
    jac(0, 1) = phiM(0, 0) * mvye[0] + phiM(1, 0) * mvye[1] + phiM(2, 0) * mvye[2];
    jac(1, 1) = phiM(0, 1) * mvye[0] + phiM(1, 1) * mvye[1] + phiM(2, 1) * mvye[2];

    double detJ = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
    ijac(0, 0) = jac(1, 1) / detJ;
    ijac(0, 1) = -jac(0, 1) / detJ;
    ijac(1, 0) = -jac(1, 0) / detJ;
    ijac(1, 1) = jac(0, 0) / detJ;
    Area0 = 0.5 * detJ;

    noalias(mdcgM) = prod(ijac, trans(ms_loc_der_central));

    //*****************************************************************************
    double eta1 = 0.5;
    double eta2 = 0.5;
    if (HasNeighbour(0, nodal_neigb[0]))
        CalculateCartesianDerOnGauss(eta1, eta2, ms_coord, mvxe, mvye, ijac, mdcg1);
    else
    {
        noalias(mdcg1) = ZeroMatrix(6, 2);
        for (unsigned int i = 0; i < 3; i++)
        {
            mdcg1(0, i) = mdcgM(0, i);
            mdcg1(1, i) = mdcgM(1, i);
        }
    }

    eta1 = 0.0;
    eta2 = 0.5;
    if (HasNeighbour(1, nodal_neigb[1]))
        CalculateCartesianDerOnGauss(eta1, eta2, ms_coord, mvxe, mvye, ijac, mdcg2);
    else
    {
        noalias(mdcg2) = ZeroMatrix(6, 2);
        for (unsigned int i = 0; i < 3; i++)
        {
            mdcg2(0, i) = mdcgM(0, i);
            mdcg2(1, i) = mdcgM(1, i);
        }
    }

    eta1 = 0.5;
    eta2 = 0.0;
    if (HasNeighbour(2, nodal_neigb[2]))
        CalculateCartesianDerOnGauss(eta1, eta2, ms_coord, mvxe, mvye, ijac, mdcg3);
    else
    {
        noalias(mdcg3) = ZeroMatrix(6, 2);
        for (unsigned int i = 0; i < 3; i++)
        {
            mdcg3(0, i) = mdcgM(0, i);
            mdcg3(1, i) = mdcgM(1, i);
        }
    }

    //initialize the initial curvature to a well defined value
    mK0[0] = 1234567.89;
    mK0[1] = 1234567.89;
    mK0[2] = 1234567.89;

	//Constitutive Law initialization

	Vector N(3);
    N[0] = 0.333333333333333333333333333;
    N[1] = 0.333333333333333333333333333;
    N[2] = 0.333333333333333333333333333;

	if ( mConstitutiveLawVector.size() != EBST_ELEM_NUM_IP )
	{
		mConstitutiveLawVector.resize( EBST_ELEM_NUM_IP );
	}
	else
	{
		// check whether the constitutive law pointers have been already set up
		bool already_set_up = true;
		for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
		{
			if(mConstitutiveLawVector[i] == NULL)
				already_set_up = false;
		}
		if(already_set_up)
		{
			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), N );
			}

			return; // if so, we are done here!
		}
	}

	// NOTE:
	// This is the standard (previous) implementation:
	// If we are here, it means that no one already set up the constitutive law vector
	// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

	if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
	{
		for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
		{
			mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
			mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(), N);
		}
	}
	else
	{
		KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
	}

	KRATOS_CATCH("");
}

void EBSTElement2D3N::CalculateCartesianDerOnGauss(
    const double eta1,
    const double eta2,
    const boost::numeric::ublas::bounded_matrix<double, 6, 3 > & ms_coord,
    const array_1d<double, 3 > & vxe,
    const array_1d<double, 3 > & vye,
    boost::numeric::ublas::bounded_matrix<double, 2, 2 > & Jinv,
    boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcg
)
{
//               boost::numeric::ublas::bounded_matrix<double, 6, 18 > msL1;
//         boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_f;
//         boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_m;
//         boost::numeric::ublas::bounded_matrix<double, 18, 18 > msK;
//         boost::numeric::ublas::bounded_matrix<double, 6, 3 > ms_coord;
    boost::numeric::ublas::bounded_matrix<double, 3, 2 > ms_loc_der_central;
    boost::numeric::ublas::bounded_matrix<double, 6, 2 > ms_loc_der_patch;

    double eta3 = 1.0 - eta2 - eta1;

    array_1d<double, 6 > N;
    N[0] = eta3 + eta1*eta2;
    N[1] = eta1 + eta2*eta3;
    N[2] = eta2 + eta3*eta1;
    N[3] = eta3 * (eta3 - 1.0)*0.5;
    N[4] = eta1 * (eta1 - 1.0)*0.5;
    N[5] = eta2 * (eta2 - 1.0)*0.5;

    ms_loc_der_patch(0, 0) = -1.0 + eta2;
    ms_loc_der_patch(0, 1) = -1.0 + eta1;
    ms_loc_der_patch(1, 0) = 1.0 - eta2;
    ms_loc_der_patch(1, 1) = 1.0 - eta1 - 2.0 * eta2;
    ms_loc_der_patch(2, 0) = -2.0 * eta1 + 1.0 - eta2;
    ms_loc_der_patch(2, 1) = 1.0 - eta1;
    ms_loc_der_patch(3, 0) = eta1 + eta2 - 0.5;
    ms_loc_der_patch(3, 1) = eta1 + eta2 - 0.5;
    ms_loc_der_patch(4, 0) = eta1 - 0.5;
    ms_loc_der_patch(4, 1) = 0.0;
    ms_loc_der_patch(5, 0) = 0.0;
    ms_loc_der_patch(5, 1) = eta2 - 0.5;

    boost::numeric::ublas::bounded_matrix<double, 3, 2 > phi;
    noalias(phi) = prod(trans(ms_coord), ms_loc_der_patch);

    array_1d<double, 3 > phi1, phi2;
    phi1[0] = phi(0, 0);
    phi1[1] = phi(1, 0);
    phi1[2] = phi(2, 0);
    phi2[0] = phi(0, 1);
    phi2[1] = phi(1, 1);
    phi2[2] = phi(2, 1);

    array_1d<double, 3 > t3g;
    MathUtils<double>::CrossProduct(t3g, phi1, phi2);
    double n3 = norm_2(t3g);
    t3g /= n3;
    // KRATOS_WATCH(t3g);

    array_1d<double, 3 > t2g;
    MathUtils<double>::CrossProduct(t2g, t3g, vxe);
    double n2 = norm_2(t2g);
    t2g /= n2;
    // KRATOS_WATCH(t2g);

    array_1d<double, 3 > t1g;
    MathUtils<double>::CrossProduct(t1g, t2g, t3g);
    double n1 = norm_2(t1g);
    t1g /= n1;
    // KRATOS_WATCH(t1g);

    boost::numeric::ublas::bounded_matrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(phi1, t1g);
    jac(0, 1) = inner_prod(phi1, t2g);
    jac(1, 0) = inner_prod(phi2, t1g);
    jac(1, 1) = inner_prod(phi2, t2g);
    // KRATOS_WATCH(jac);

    double detJ = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
    Jinv(0, 0) = jac(1, 1) / detJ;
    Jinv(0, 1) = -jac(0, 1) / detJ;
    Jinv(1, 0) = -jac(1, 0) / detJ;
    Jinv(1, 1) = jac(0, 0) / detJ;

    noalias(dcg) = prod(Jinv, trans(ms_loc_der_patch));

}

void EBSTElement2D3N::CalculatePhiM(
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiM,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & dcM,
    const boost::numeric::ublas::bounded_matrix<double, 6, 3 > & ms_coord
)
{
    noalias(phiM) = ZeroMatrix(2, 3);
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            phiM(0, j) += dcM(0, i) * ms_coord(i, j);
            phiM(1, j) += dcM(1, i) * ms_coord(i, j);
            //	      phiM(j,0) += dcM(i,0)*ms_coord(i,j);
            //	      phiM(j,1) += dcM(i,1)*ms_coord(i,j);
        }
    }
}

void EBSTElement2D3N::CalculatePhiG(
    boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG,
    const boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcG,
    const boost::numeric::ublas::bounded_matrix<double, 6, 3 > & ms_coord
)
{
    noalias(phiG) = ZeroMatrix(2, 3);
    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            phiG(0, j) += dcG(0, i) * ms_coord(i, j);
            phiG(1, j) += dcG(1, i) * ms_coord(i, j);
            //	      phiG(j,0) += dcG(i,0)*ms_coord(i,j);
            //	      phiG(j,1) += dcG(i,1)*ms_coord(i,j);
        }
    }
}

void EBSTElement2D3N::CalculateAndAdd_MembraneStrain(
    array_1d<double, 3 > & strain,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG
)
{
    strain[0] += 0.5 * (phiG(0, 0) * phiG(0, 0) + phiG(0, 1) * phiG(0, 1) + phiG(0, 2) * phiG(0, 2) - 1.0);
    strain[1] += 0.5 * (phiG(1, 0) * phiG(1, 0) + phiG(1, 1) * phiG(1, 1) + phiG(1, 2) * phiG(1, 2) - 1.0);
    strain[2] += phiG(0, 0) * phiG(1, 0) + phiG(0, 1) * phiG(1, 1) + phiG(0, 2) * phiG(1, 2);
}

void EBSTElement2D3N::CalculateAndAdd_MembraneB(
    boost::numeric::ublas::bounded_matrix<double, 3, 18 > & B_m,
    const boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcgG,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG
)
{
    //implementing as in Eqn 75
    boost::numeric::ublas::bounded_matrix<double, 6, 18 > aux;
    noalias(aux) = ZeroMatrix(6, 18);
    for (unsigned int i = 0; i < 6; i++)
    {
        unsigned int base = i * 3;
        aux(0, base) = dcgG(0, i);
        aux(1, base + 1) = dcgG(0, i);
        aux(2, base + 2) = dcgG(0, i);
        aux(3, base) = dcgG(1, i);
        aux(4, base + 1) = dcgG(1, i);
        aux(5, base + 2) = dcgG(1, i);
    }

    boost::numeric::ublas::bounded_matrix<double, 3, 6 > aux_phi;
    noalias(aux_phi) = ZeroMatrix(3, 6);
    aux_phi(0, 0) = phiG(0, 0);
    aux_phi(0, 1) = phiG(0, 1);
    aux_phi(0, 2) = phiG(0, 2);

    aux_phi(1, 3) = phiG(1, 0);
    aux_phi(1, 4) = phiG(1, 1);
    aux_phi(1, 5) = phiG(1, 2);

    aux_phi(2, 0) = phiG(1, 0);
    aux_phi(2, 1) = phiG(1, 1);
    aux_phi(2, 2) = phiG(1, 2);
    aux_phi(2, 3) = phiG(0, 0);
    aux_phi(2, 4) = phiG(0, 1);
    aux_phi(2, 5) = phiG(0, 2);

    noalias(B_m) += prod(aux_phi, aux);
}

void EBSTElement2D3N::Calculate_h_ab(
    array_1d<double, 3 > & h_ab,
    const unsigned int alpha,
    const unsigned int beta,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG1,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG2,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & phiG3,
    const boost::numeric::ublas::bounded_matrix<double, 2, 3 > & dcgM
)
{
    //GAUSS 1
    double Nalpha = dcgM(alpha, 0);
    double Nbeta = dcgM(beta, 0);
    h_ab[0] = Nalpha * phiG1(beta, 0) + Nbeta * phiG1(alpha, 0);
    h_ab[1] = Nalpha * phiG1(beta, 1) + Nbeta * phiG1(alpha, 1);
    h_ab[2] = Nalpha * phiG1(beta, 2) + Nbeta * phiG1(alpha, 2);

    //GAUSS 2
    Nalpha = dcgM(alpha, 1);
    Nbeta = dcgM(beta, 1);
    h_ab[0] += Nalpha * phiG2(beta, 0) + Nbeta * phiG2(alpha, 0);
    h_ab[1] += Nalpha * phiG2(beta, 1) + Nbeta * phiG2(alpha, 1);
    h_ab[2] += Nalpha * phiG2(beta, 2) + Nbeta * phiG2(alpha, 2);

    //GAUSS 3
    Nalpha = dcgM(alpha, 2);
    Nbeta = dcgM(beta, 2);
    h_ab[0] += Nalpha * phiG3(beta, 0) + Nbeta * phiG3(alpha, 0);
    h_ab[1] += Nalpha * phiG3(beta, 1) + Nbeta * phiG3(alpha, 1);
    h_ab[2] += Nalpha * phiG3(beta, 2) + Nbeta * phiG3(alpha, 2);

}

void EBSTElement2D3N::CalculateAndAdd_Membrane_Kg(
    boost::numeric::ublas::bounded_matrix<double, 18, 18 > & K,
    const boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcgG1,
    const boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcgG2,
    const boost::numeric::ublas::bounded_matrix<double, 2, 6 > & dcgG3,
    const array_1d<double, 3 > & membrane_stress
)
{
    boost::numeric::ublas::bounded_matrix<double, 6, 6 > A;


    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 6; j++)
        {
            //gauss 1
            A(i, j) = membrane_stress[0] * dcgG1(0, i) * dcgG1(0, j)
                      + membrane_stress[1] * dcgG1(1, i) * dcgG1(1, j)
                      + membrane_stress[2]*(dcgG1(0, i) * dcgG1(1, j) + dcgG1(1, i) * dcgG1(0, j));

            //gauss 2
            A(i, j) += membrane_stress[0] * dcgG2(0, i) * dcgG2(0, j)
                       + membrane_stress[1] * dcgG2(1, i) * dcgG2(1, j)
                       + membrane_stress[2]*(dcgG2(0, i) * dcgG2(1, j) + dcgG2(1, i) * dcgG2(0, j));

            //gauss 3
            A(i, j) += membrane_stress[0] * dcgG3(0, i) * dcgG3(0, j)
                       + membrane_stress[1] * dcgG3(1, i) * dcgG3(1, j)
                       + membrane_stress[2]*(dcgG3(0, i) * dcgG3(1, j) + dcgG3(1, i) * dcgG3(0, j));
        }
    }

    A *= 0.333333333333333333333333333;

    //now assemble in K
    for (unsigned int i = 0; i < 6; i++)
    {
        unsigned int base_i = i * 3;
        for (unsigned int j = 0; j < 6; j++)
        {
            unsigned int base_j = j * 3;
            K(base_i, base_j) += A(i, j);
            K(base_i + 1, base_j + 1) += A(i, j);
            K(base_i + 2, base_j + 2) += A(i, j);
        }
    }


}

void EBSTElement2D3N::CalculateEquivalentStresses(
    const array_1d<double, 3 > & membrane_strain_el,
    const array_1d<double, 3 > & bending_strain_el,
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Dmat_m,
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Dmat_f,
    array_1d<double, 3 > & membrane_stress,
    array_1d<double, 3 > & bending_stress,
    double& h_on_h0, //  ratio between current thickness and original thickness h/h0
    const ProcessInfo& rCurrentProcessInfo,
	const array_1d<double, 3 > & elem_vx
)
{
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculate the angle between material and local element axis

	// calculate the angle between the element x direction and the material x direction.
	const array_1d<double, 3 >& a = elem_vx;
	array_1d<double, 3 > b; // material direction. here in 2d [1,0,0]
	b[0]=1.0; b[1]=0.0; b[2]=0.0;
	double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
	if(a_dot_b < -1.0) a_dot_b = -1.0;
	if(a_dot_b >  1.0) a_dot_b =  1.0;
	double angle = std::acos( a_dot_b );
	// if they are not counter-clock-wise, let's change the sign of the angle
	if(angle != 0.0) {
		array_1d<double, 3 > elem_vz; elem_vz[0]=0.0; elem_vz[1]=0.0; elem_vz[2]=1.0;
		array_1d<double, 3> elem_vy;
		MathUtils<double>::CrossProduct(elem_vy, elem_vz, elem_vx);
		if( b(0)*elem_vy(0) + b(1)*elem_vy(1) + b(2)*elem_vy(2) < 0.0 )
			angle = -angle;
	}
	Matrix RE;
	GetRotationMatrixForStrains(-angle,RE);
	Matrix RS;
	GetRotationMatrixForStresses(angle,RS);

	// rotate the elemental strains to the material coordinate system
	array_1d<double,3> membrane_strain;
	array_1d<double,3> bending_strain;
	noalias(membrane_strain) = prod(RE, membrane_strain_el);
	noalias(bending_strain) = prod(RE, bending_strain_el);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    //compute current area
    double current_area = GetGeometry().Area();

    //compute h/h0 in the hipothesis of incompressible behaviour (this is indeed an approx)
    double lambda = Area0 / current_area;
    h_on_h0 = lambda;

    //define aux variables
    double t = GetProperties()[THICKNESS];
    double t_half = 0.5*t;
    double eta,membrane_weight,stress_bending_weight,Df_weight;
    Vector strain(3);
    Vector stress(3);
    Matrix D(3,3);
    Vector N(3);
    N[0] = 0.333333333333333333333333333;
    N[1] = 0.333333333333333333333333333;
    N[2] = 0.333333333333333333333333333;

	// constitutive law parameters
	ConstitutiveLaw::Parameters parameters(GetGeometry(), GetProperties(), rCurrentProcessInfo);
	Flags& options = parameters.GetOptions();
	options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
	options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
	double detF = 1.0;
	Matrix F(IdentityMatrix(2,2));
	Matrix DNDX(ZeroMatrix(3,2));
	parameters.SetDeterminantF(detF);
	parameters.SetDeformationGradientF(F);
	parameters.SetShapeFunctionsValues(N);
	parameters.SetShapeFunctionsDerivatives(DNDX);
	parameters.SetStrainVector( strain );
	parameters.SetStressVector( stress );
	parameters.SetConstitutiveMatrix( D );

    //gauss point 1
    eta = -0.8611363116   ;
    membrane_weight =  0.3478548451*t_half*Area0;
    stress_bending_weight = membrane_weight*(eta)*t_half;
    Df_weight = membrane_weight*(eta*eta)*t_half*t_half;
    strain[0] = membrane_strain[0] + bending_strain[0]*eta*lambda*t_half;
    strain[1] = membrane_strain[1] + bending_strain[1]*eta*lambda*t_half;
    strain[2] = membrane_strain[2] + bending_strain[2]*eta*lambda*t_half; //note that the 2 is already included in the membrane and bending strain;
    mstrains(0,0) = strain[0];
    mstrains(0,1) = strain[1];
    mstrains(0,2) = strain[2];
    //mConstitutiveLawVector[0]->CalculateMaterialResponse(strain,ZeroMatrix(1,1),stress,D,rCurrentProcessInfo,GetProperties(),GetGeometry(),N,true,1,false);
    mConstitutiveLawVector[0]->CalculateMaterialResponseCauchy(parameters);
	noalias(Dmat_m)             = membrane_weight * D;
    noalias(membrane_stress)    = membrane_weight * stress;
    noalias(Dmat_f)             = Df_weight * D;
    noalias(bending_stress)     = stress_bending_weight * stress;

    //gauss point 2
    eta = -0.3399810436   ;
    membrane_weight =  0.6521451549*t_half*Area0;
    stress_bending_weight = membrane_weight*(eta)*t_half;
    Df_weight = membrane_weight*(eta*eta)*t_half*t_half;
    strain[0] = membrane_strain[0] + bending_strain[0]*eta*lambda*t_half;
    strain[1] = membrane_strain[1] + bending_strain[1]*eta*lambda*t_half;
    strain[2] = membrane_strain[2] + bending_strain[2]*eta*lambda*t_half; //note that the 2 is already included in the membrane and bending strain;
    mstrains(1,0) = strain[0];
    mstrains(1,1) = strain[1];
    mstrains(1,2) = strain[2];
#ifndef EBST_SIMULATE_2D
	mConstitutiveLawVector[1]->CalculateMaterialResponseCauchy(parameters);  
#endif // !EBST_SIMULATE_2D
	noalias(Dmat_m)             += membrane_weight * D;
    noalias(membrane_stress)    += membrane_weight * stress;
    noalias(Dmat_f)             += Df_weight * D;
    noalias(bending_stress)     += stress_bending_weight * stress;

    //gauss point 3
    eta =  0.3399810436   ;
    membrane_weight =  0.6521451549*t_half*Area0;
    stress_bending_weight = membrane_weight*(eta)*t_half;
    Df_weight = membrane_weight*eta*eta*t_half*t_half;
    strain[0] = membrane_strain[0] + bending_strain[0]*eta*lambda*t_half;
    strain[1] = membrane_strain[1] + bending_strain[1]*eta*lambda*t_half;
    strain[2] = membrane_strain[2] + bending_strain[2]*eta*lambda*t_half; //note that the 2 is already included in the membrane and bending strain;
    mstrains(2,0) = strain[0];
    mstrains(2,1) = strain[1];
    mstrains(2,2) = strain[2];
#ifndef EBST_SIMULATE_2D
	mConstitutiveLawVector[2]->CalculateMaterialResponseCauchy(parameters);  
#endif // !EBST_SIMULATE_2D
	noalias(Dmat_m)             += membrane_weight * D;
    noalias(membrane_stress)    += membrane_weight * stress;
    noalias(Dmat_f)             += Df_weight * D;
    noalias(bending_stress)     += stress_bending_weight * stress;

    //gauss point 4
    eta = 0.8611363116   ;
    membrane_weight =  0.3478548451*t_half*Area0;
    stress_bending_weight = membrane_weight*(eta)*t_half;
    Df_weight = membrane_weight*eta*eta*t_half*t_half;
    strain[0] = membrane_strain[0] + bending_strain[0]*eta*lambda*t_half;
    strain[1] = membrane_strain[1] + bending_strain[1]*eta*lambda*t_half;
    strain[2] = membrane_strain[2] + bending_strain[2]*eta*lambda*t_half; //note that the 2 is already included in the membrane and bending strain;
    mstrains(3,0) = strain[0];
    mstrains(3,1) = strain[1];
    mstrains(3,2) = strain[2];
#ifndef EBST_SIMULATE_2D
	mConstitutiveLawVector[3]->CalculateMaterialResponseCauchy(parameters);  
#endif // !EBST_SIMULATE_2D
	noalias(Dmat_m)             += membrane_weight * D;
    noalias(membrane_stress)    += membrane_weight * stress;
    noalias(Dmat_f)             += Df_weight * D;
    noalias(bending_stress)     += stress_bending_weight * stress;


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// back to elemental coordinate system
	Matrix aux_1(3,3);
	membrane_stress = prod(RS, membrane_stress);
	bending_stress = prod(RS, bending_stress);

	noalias(aux_1) = prod(Dmat_m, trans(RS));
	noalias(Dmat_m) = prod(RS, aux_1);
	noalias(aux_1) = prod(Dmat_f, trans(RS));
	noalias(Dmat_f) = prod(RS, aux_1);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}

//auxiliary function needed in the calculation of output stresses
inline array_1d<double,6> EBSTElement2D3N::VoigtTensorComponents(
    array_1d<double,3>& a,
    array_1d<double,3>& b)

{
    array_1d<double,6> v;

    v[0] = a[0]*b[0];
    v[1] = a[1]*b[1];
    v[2] = a[2]*b[2];
    v[3] = a[0]*b[1];
    v[4] = a[1]*b[2];
    v[5] = a[0]*b[2];

    return v;
}

void EBSTElement2D3N::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    Vector N(3);
    N[0] = 0.33333333333333333333;
    N[1] = 0.33333333333333333333;
    N[2] = 0.33333333333333333333;
    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),GetGeometry(),N, CurrentProcessInfo);
    }
}

} // Namespace Kratos.

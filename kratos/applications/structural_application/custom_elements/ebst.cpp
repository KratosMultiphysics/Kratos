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
#include "custom_elements/ebst.h"
#include "includes/constitutive_law.h"
#include "structural_application.h"


namespace Kratos {

    namespace EbstAuxiliaries {
        boost::numeric::ublas::bounded_matrix<double, 6, 18 > msL1;
        boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_f;
        boost::numeric::ublas::bounded_matrix<double, 3, 18 > msB_m;
        boost::numeric::ublas::bounded_matrix<double, 18, 18 > msK;
        
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > msDmat_m;	
        boost::numeric::ublas::bounded_matrix<double, 6, 3 > ms_coord;

        boost::numeric::ublas::bounded_matrix<double, 3, 2 > ms_loc_der_central;
        boost::numeric::ublas::bounded_matrix<double, 6, 2 > ms_loc_der_patch;

    }

    using namespace EbstAuxiliaries;


    //***********************************************************************************
    //***********************************************************************************
    // -------- //
    //  PUBLIC  //
    // -------- //

    // Constructor

    Ebst::Ebst() {
    }

    // Constructor

    Ebst::Ebst(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
    }

    // Constructor

    Ebst::Ebst(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
    }

    //***********************************************************************************
    //***********************************************************************************

    Element::Pointer Ebst::Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new Ebst(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    //***********************************************************************************
    //***********************************************************************************
    // Destructor

    Ebst::~Ebst() {
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);

        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);
        unsigned int dim = number_of_nodes * 3;

        if (rResult.size() != dim)
            rResult.resize(dim, false);

        //nodes of the central element
        for (int i = 0; i < 3; i++) {
            int index = i*3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

//        KRATOS_WATCH(Id());
//
//        for(unsigned int k = 0; k<rResult.size(); k++)
//            std::cout << rResult[k] << " ";
//        std::cout << std::endl;

        //adding the ids ofthe neighbouring nodes
        int index = 9;
        for (unsigned int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
                rResult[index] = neigb[i].GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = neigb[i].GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = neigb[i].GetDof(DISPLACEMENT_Z).EquationId();
                index += 3;
            }
        }
        
//        for(unsigned int k = 0; k<rResult.size(); k++)
//            std::cout << rResult[k] << " ";
//        std::cout << std::endl;

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void Ebst::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
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
            if (HasNeighbour(i,neigb[i]))
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

    void Ebst::GetValuesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < 3; i++) {
            const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * 3;
            values[index] = disp[0];
            values[index + 1] = disp[1];
            values[index + 2] = disp[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
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

    void Ebst::GetFirstDerivativesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < 3; i++) {
            const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            unsigned int index = i * 3;
            values[index] = vel[0];
            values[index + 1] = vel[1];
            values[index + 2] = vel[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
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

    void Ebst::GetSecondDerivativesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < 3; i++) {
            const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            unsigned int index = i * 3;
            values[index] = acc[0];
            values[index + 1] = acc[1];
            values[index + 2] = acc[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
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

    void Ebst::CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateOnIntegrationPoints(
            const Variable<Matrix>& rVariable,
            std::vector<Matrix>& Output,
            const ProcessInfo& rCurrentProcessInfo) {

    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::MassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);

        //rMassMatrix.resize(0,0);
        // LUMPED MASS MATRIX
        unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
        unsigned int MatSize = number_of_nodes * 3;
        if (rMassMatrix.size1() != MatSize)
            rMassMatrix.resize(MatSize, MatSize, false);
        rMassMatrix = ZeroMatrix(MatSize, MatSize);

        double Area = GetGeometry().Area();
        double TotalMass = Area * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
        Vector LumpFact;
        LumpFact = GetGeometry().LumpingFactors(LumpFact);

        for (unsigned int i = 0; i < number_of_nodes; i++) {
            double temp = LumpFact[i] * TotalMass;
            for (unsigned int j = 0; j < 3; j++) {
                unsigned int index = i * 3 + j;
                rMassMatrix(index, index) = temp;
            }
        }

        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::DampMatrix(
            MatrixType& rDampMatrix,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

        if (rDampMatrix.size1() != 0)
            rDampMatrix.resize(0, 0, false);

        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::FinalizeSolutionStep(
            ProcessInfo& rCurrentProcessInfo) {
    }





    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag,
            bool CalculateResidualVectorFlag) {
        KRATOS_TRY

	WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
	
// KRATOS_WATCH(Id());
	
	array_1d<double,3> v12,v13,vze,vxe,vye;
	
	//fill the aux matrix of coordinates
	for(unsigned int i=0;i<3; i++) { 
	    ms_coord(i,0)=GetGeometry()[i].X(); 
	    ms_coord(i,1)=GetGeometry()[i].Y(); 
	    ms_coord(i,2)=GetGeometry()[i].Z();  
	}
	for(unsigned int i=0;i<3; i++) { 
	    ms_coord(i+3,0)=neigb[i].X(); 
	    ms_coord(i+3,1)=neigb[i].Y(); 
	    ms_coord(i+3,2)=neigb[i].Z();  
	}
// KRATOS_WATCH(coord);
	
	//calculate local system of coordinates of the elem ->SysCartE
	v12[0] = GetGeometry()[1].X() -  GetGeometry()[0].X();
	v12[1] = GetGeometry()[1].Y() -  GetGeometry()[0].Y();
	v12[2] = GetGeometry()[1].Z() -  GetGeometry()[0].Z();
	double n12 = norm_2(v12);
	v12/=n12;
	
	v13[0] = GetGeometry()[2].X() -  GetGeometry()[0].X();
	v13[1] = GetGeometry()[2].Y() -  GetGeometry()[0].Y();
	v13[2] = GetGeometry()[2].Z() -  GetGeometry()[0].Z();
	double n13 = norm_2(v13);
	v13/=n13;
	
	//vze = cross prod v12 v13 --> then normalize
	MathUtils<double>::CrossProduct(vze,v12,v13);
	double nze = norm_2(vze);
	vze/=nze;
	
	//version by Riccardo
	noalias(vxe) = v12;
	
	MathUtils<double>::CrossProduct(vye,vze,vxe);
	
// am
	
/*	KRATOS_WATCH(v12);
	KRATOS_WATCH(v13);
	KRATOS_WATCH(vze);
	KRATOS_WATCH(vxe);
	KRATOS_WATCH(vye);*/
	
	//*****************************************************************************
	//calculate cartesian derivatives for the central element
	ms_loc_der_central(0,0) = -1.0  ; ms_loc_der_central(0,1) = -1.0  ;
	ms_loc_der_central(1,0) =  1.0  ; ms_loc_der_central(1,1) =  0.0  ;
	ms_loc_der_central(2,0) =  0.0  ; ms_loc_der_central(2,1) =  1.0  ;
	boost::numeric::ublas::bounded_matrix<double, 3,2 > dd = ZeroMatrix(3,2);
	for(unsigned int i=0; i<3; i++) {
	  for(unsigned int j=0; j<3; j++){
	      dd(j,0) += ms_loc_der_central(i,0)*ms_coord(i,j);
	      dd(j,1) += ms_loc_der_central(i,1)*ms_coord(i,j);
	    }
	 }

	boost::numeric::ublas::bounded_matrix<double, 2,2 > jac;
	jac(0,0) = dd(0,0)*vxe[0] + dd(1,0)*vxe[1] + dd(2,0)*vxe[2];
	jac(1,0) = dd(0,1)*vxe[0] + dd(1,1)*vxe[1] + dd(2,1)*vxe[2];
	jac(0,1) = dd(0,0)*vye[0] + dd(1,0)*vye[1] + dd(2,0)*vye[2];
	jac(1,1) = dd(0,1)*vye[0] + dd(1,1)*vye[1] + dd(2,1)*vye[2];
	
	double detJ = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
	boost::numeric::ublas::bounded_matrix<double, 2,2 > ijac;
	ijac(0,0) =  jac(1,1) / detJ;
	ijac(0,1) = -jac(0,1) / detJ;
	ijac(1,0) = -jac(1,0) / detJ;
	ijac(1,1) =  jac(0,0) / detJ;
	double Area = 0.5 * detJ;
// KRATOS_WATCH(Area);
	boost::numeric::ublas::bounded_matrix<double, 3,2 > dcgM = prod(ijac,trans(ms_loc_der_central));

// KRATOS_WATCH(dcgM);
	//*****************************************************************************
	boost::numeric::ublas::bounded_matrix<double, 3,2 > phi;	
	
	double eta1 = 0.5   ; double eta2 = 0.5;
	boost::numeric::ublas::bounded_matrix<double, 2,6 > dcg1; //cartesian derivatives on gauss 1
	if (HasNeighbour(0,neigb[0])) 
	  CalculateCartesianDerOnGauss(eta1,eta2, ms_coord,vxe,vye,phi,dcg1 );
	else{
	  noalias(dcg1) = ZeroMatrix(2,6);
	  for(unsigned int i=0;i<3;i++){
	    dcg1(0,i) = dcgM(0,i); dcg1(1,i) = dcgM(1,i);
	  }
	}
// KRATOS_WATCH(phi);
// KRATOS_WATCH(dcg1);

	eta1 = 0.0   ; eta2 = 0.5;
	boost::numeric::ublas::bounded_matrix<double, 2,6 > dcg2; //cartesian derivatives on gauss 2
	if (HasNeighbour(1,neigb[1])) 
	  CalculateCartesianDerOnGauss(eta1,eta2, ms_coord,vxe,vye,phi,dcg2 );
	else{
	  noalias(dcg2) = ZeroMatrix(2,6);
	  for(unsigned int i=0;i<3;i++){
	    dcg2(0,i) = dcgM(0,i); dcg2(1,i) = dcgM(1,i);
	  }
	}
// KRATOS_WATCH(phi);
// KRATOS_WATCH(dcg2);

	eta1 = 0.5   ; eta2 = 0.0;
	boost::numeric::ublas::bounded_matrix<double, 2,6 > dcg3; //cartesian derivatives on gauss 3
	if (HasNeighbour(2,neigb[2])) 
	  CalculateCartesianDerOnGauss(eta1,eta2, ms_coord,vxe,vye,phi,dcg3 );
	else{
	  noalias(dcg3) = ZeroMatrix(2,6);
	  for(unsigned int i=0;i<3;i++){
	    dcg3(0,i) = dcgM(0,i); dcg3(1,i) = dcgM(1,i);
	  }
	}
// KRATOS_WATCH(phi);
// KRATOS_WATCH(dcg3);

	//*****************************************************************************
	//calculating derivative of h --> see Eqn 62
	unsigned int ind = 0;
	for(unsigned int i = 0; i < 6; i++)
	{
	  for(unsigned int j = 0; j < 3; j++)
	  {
	    msL1(0,ind) = dcg1(0,i)*vze[j];
	    msL1(1,ind) = dcg1(1,i)*vze[j];

	    msL1(2,ind) = dcg2(0,i)*vze[j];
	    msL1(3,ind) = dcg2(1,i)*vze[j];

	    msL1(4,ind) = dcg3(0,i)*vze[j];
	    msL1(5,ind) = dcg3(1,i)*vze[j];
	    ind++;
	  }
	}
// KRATOS_WATCH(msL1);
	
	boost::numeric::ublas::bounded_matrix<double, 3,6 > DN = ZeroMatrix(3,6);
	DN(0,0) = dcgM(0,0);
	DN(1,1) = dcgM(1,0);
	DN(2,0) = dcgM(1,0);
	DN(2,1) = dcgM(0,0);
	DN(0,2) = dcgM(0,1);
	DN(1,3) = dcgM(1,1);
	DN(2,2) = dcgM(1,1);
	DN(2,3) = dcgM(0,1);
	DN(0,4) = dcgM(0,2);
	DN(1,5) = dcgM(1,2);
	DN(2,4) = dcgM(1,2);
	DN(2,5) = dcgM(0,2);
// KRATOS_WATCH(DN);

	noalias(msB_f) = 2.0*prod(DN,msL1);
// KRATOS_WATCH(msB_f);


	//flexural contribution
	double E = GetProperties()[YOUNG_MODULUS];
	double nu = GetProperties()[POISSON_RATIO];
	double t = GetProperties()[THICKNESS];
	double aux0 = t*t*t/12.0;
	double aux1 = aux0 * E / (1.0 - nu*nu);
	double aux2 = aux1 * nu;
	double aux3 = aux0 * E * 0.5 / (1+nu);
	boost::numeric::ublas::bounded_matrix<double, 3, 3 > msDmat_f;
	msDmat_f(0,0) = aux1  ; msDmat_f(0,1) = aux2  ; msDmat_f(0,2) = 0.0  ;
	msDmat_f(1,0) = aux2  ; msDmat_f(1,1) = aux1  ; msDmat_f(1,2) = 0.0  ;
	msDmat_f(2,0) = 0.0   ; msDmat_f(2,1) = 0.0   ; msDmat_f(2,2) = aux3 ;
// KRATOS_WATCH(msDmat_f);	
	msDmat_f *= Area;
	
	boost::numeric::ublas::bounded_matrix<double, 3,18 > aux;
	noalias(aux) = prod(msDmat_f,msB_f);
	noalias(msK) = prod(trans(msB_f),aux);
	
	unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
        unsigned int MatSize = number_of_nodes * 3;

// KRATOS_WATCH(msK);
// KRATOS_WATCH(MatSize);
        //resizing as needed the LHS

            if (rLeftHandSideMatrix.size1() != MatSize)
                rLeftHandSideMatrix.resize(MatSize, MatSize);
             noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
	    
	    array_1d<unsigned int,18> id_vec;
	    for(unsigned int i=0; i<9;i++) id_vec[i] = i;
	    unsigned int index = 9;
	    for(unsigned int i=0; i<3;i++)
	    {
	      if (HasNeighbour(i,neigb[i]))
	      {
		for(unsigned int j=0; j<3;j++)
		  id_vec[9+i*3+j] = index+j;
		index += 3;
	      }
	      else
	      {
		for(unsigned int j=0; j<3;j++)
		  id_vec[9+i*3+j] = 1000;
	      }
	    }
	      
	    
	    
	    //add the first 9*9 block
	    for(unsigned int i=0; i<18;i++)
	    {
	      if(id_vec[i] < 1000)
	      {
	      for(unsigned int j=0; j<18; j++)
		{
		  if(id_vec[j] < 1000)  rLeftHandSideMatrix(id_vec[i],id_vec[j]) = msK(i,j);
		}
	      }
	    }
	    
// KRATOS_WATCH(id_vec);
	    
	    
        

            if (rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize);
            rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
	    
	    //apply body force
	    array_1d<double,3> aaa;
	    noalias(aaa) = GetProperties()[BODY_FORCE];
	    aaa *= Area*0.3333333333333333333333;
	    for(unsigned int i = 0; i<3; i++)
	      for(unsigned int j=0; j<3; j++)
		rRightHandSideVector[i*3+j] = aaa[j];
// KRATOS_WATCH(rRightHandSideVector);
	    
	    //RHS -= K*disp;
	    Vector values(MatSize);
	    GetValuesVector(values,0);
// KRATOS_WATCH(values);
	    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,values);
/*KRATOS_WATCH(Id());
KRATOS_WATCH(rLeftHandSideMatrix);
KRATOS_WATCH(rRightHandSideVector);*/
	    
 

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true) {
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
        }

        KRATOS_CATCH("");
    }


    //***********************************************************************************
    //***********************************************************************************
    bool Ebst::HasNeighbour(unsigned int index, const Node < 3 > & neighb)
    {
        if (neighb.Id() == GetGeometry()[index].Id())
            return false;
        else
            return true;
    }

    //***********************************************************************************
    //***********************************************************************************

    unsigned int Ebst::NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs)
    {
        unsigned int active_neighbours = 0;
        for (unsigned int i = 0; i < neighbs.size(); i++)
            if (HasNeighbour(i,neighbs[i])) active_neighbours++;
        return active_neighbours;
    }

    //***********************************************************************************
    //***********************************************************************************
    void Ebst::Initialize()
    {
        KRATOS_TRY
        //find the "nodal neighbours" given the elemental neighbours
        WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
        if(elem_neigb.size() == 0) KRATOS_ERROR(std::logic_error,"the neighbour elements are not calculated","")
        WeakPointerVector< Node<3> >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
        nodal_neigb.resize(3);
        Geometry< Node<3> >& center_geom = GetGeometry();

        std::cout << "I am elem" << Id() <<std::endl;
        std::cout << "neighbours =" << elem_neigb[0].Id() << " " << elem_neigb[1].Id() << " " << elem_neigb[2].Id() << " " << std::endl;
        for (unsigned int i = 0; i < center_geom.size(); i++)
        {
            if(elem_neigb[i].Id() != Id() ) //if the elemental neighbour exists
            {
                Geometry< Node<3> >& geom = elem_neigb[i].GetGeometry();
                for (unsigned int j = 0; j < geom.size(); j++)
                {
                    bool aux = false;
                    for (unsigned int k = 0; k < center_geom.size(); k++)
                    {
                        if(geom[j].Id()==center_geom[k].Id())
                            aux = true;
                    }

                    if(aux == false) nodal_neigb(i) = Node<3>::WeakPointer( geom(j) );
                }
            }
            else //the elemenetal neighbour does not exist
	      nodal_neigb(i) = Node<3>::WeakPointer( center_geom(i) );
        }

        std::cout << "node1" << GetGeometry()[0].Id() << "opposite node =" << nodal_neigb[0].Id() << std::endl;
        std::cout << "node2" << GetGeometry()[1].Id() << "opposite node =" << nodal_neigb[1].Id() << std::endl;
        std::cout << "node3" << GetGeometry()[2].Id() << "opposite node =" << nodal_neigb[2].Id() << std::endl;

        KRATOS_CATCH("");
    }
    
    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateCartesianDerOnGauss(
	  const double eta1, 
	  const double eta2, 
	  const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& ms_coord,
	  const array_1d<double,3>& vxe,
	  const array_1d<double,3>& vye,
	  boost::numeric::ublas::bounded_matrix<double, 3,2 >& phi,
	  boost::numeric::ublas::bounded_matrix<double, 2,6 >& dcg
	  )
    {
	double eta3 = 1.0 - eta2 - eta1;
	
	array_1d<double,6> N;
	N[0] = eta3 + eta1*eta2;
	N[1] = eta1 + eta2*eta3;
	N[2] = eta2 + eta3*eta1;
	N[3] = eta3*(eta3-1.0)*0.5;
	N[4] = eta1*(eta1-1.0)*0.5;
	N[5] = eta2*(eta2-1.0)*0.5;
	
	ms_loc_der_patch(0,0) = -1.0 + eta2             ; ms_loc_der_patch(0,1) = -1.0 + eta1 ;
	ms_loc_der_patch(1,0) =  1.0 - eta2             ; ms_loc_der_patch(1,1) =  1.0 - eta1 -2.0*eta2 ;
	ms_loc_der_patch(2,0) =  -2.0*eta1 + 1.0 -eta2  ; ms_loc_der_patch(2,1) =  1.0 - eta1 ;
	ms_loc_der_patch(3,0) =  eta1+eta2-0.5  	  ; ms_loc_der_patch(3,1) =  eta1 + eta2 - 0.5  ;
	ms_loc_der_patch(4,0) =  eta1-0.5  		  ; ms_loc_der_patch(4,1) =  0.0  ;
	ms_loc_der_patch(5,0) =  0.0  		  ; ms_loc_der_patch(5,1) =  eta2 - 0.5  ;
	
// KRATOS_WATCH(ms_loc_der_patch);
	noalias(phi) = prod(trans(ms_coord),ms_loc_der_patch);

	array_1d<double,3> phi1,phi2;
	phi1[0] = phi(0,0); phi1[1] = phi(1,0); phi1[2] = phi(2,0);
	phi2[0] = phi(0,1); phi2[1] = phi(1,1); phi2[2] = phi(2,1);
	
	array_1d<double,3> t3g;
	MathUtils<double>::CrossProduct(t3g,phi1,phi2);
	double n3 = norm_2(t3g);
	t3g /= n3;
// KRATOS_WATCH(t3g);
	
	array_1d<double,3> t2g;
	MathUtils<double>::CrossProduct(t2g,t3g,vxe);
	double n2 = norm_2(t2g);
	t2g /= n2;	
// KRATOS_WATCH(t2g);
	
	array_1d<double,3> t1g;
	MathUtils<double>::CrossProduct(t1g,t2g,t3g);
	double n1 = norm_2(t1g);
	t1g /= n1;	
// KRATOS_WATCH(t1g);
		
	boost::numeric::ublas::bounded_matrix<double, 2,2 > jac;
	jac(0,0) = inner_prod(phi1,t1g); jac(0,1) = inner_prod(phi1,t2g);
	jac(1,0) = inner_prod(phi2,t1g); jac(1,1) = inner_prod(phi2,t2g);
// KRATOS_WATCH(jac);
	
	double detJ = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
	boost::numeric::ublas::bounded_matrix<double, 2,2 > ijac;
	ijac(0,0) =  jac(1,1) / detJ;
	ijac(0,1) = -jac(0,1) / detJ;
	ijac(1,0) = -jac(1,0) / detJ;
	ijac(1,1) =  jac(0,0) / detJ;
	
	noalias(dcg) = prod(ijac,trans(ms_loc_der_patch));
	
    }
    


} // Namespace Kratos.

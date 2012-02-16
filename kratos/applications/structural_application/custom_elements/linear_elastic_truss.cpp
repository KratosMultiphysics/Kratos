/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
1
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
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pavel $
//   Date:                $Date: 2011-12-20 17:14:42 $
//   Revision:            $Revision: 1.0 $
//

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/linear_elastic_truss.h"
#include "utilities/math_utils.h"
#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"

namespace Kratos
{
    
    /**
    * Constructor.
    * This deals without DOFs
    * @param NewId element ID
    * @param pGeometry geometry pointer
    * @see LinearElasticTruss(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    */
    LinearElasticTruss::LinearElasticTruss(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {       
        //DO NOT ADD DOFS HERE!!!
    }

    /**
    * Constructor.
    * This deals with DOFs
    * @param NewId element ID
    * @param pGeometry geometry pointer
    * @param pProperties properties pointer
    * @see LinearElasticTruss(IndexType NewId, GeometryType::Pointer pGeometry)
    */
    LinearElasticTruss::LinearElasticTruss(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
        dimension = GetGeometry().WorkingSpaceDimension();
        number_of_nodes = GetGeometry().size();

        //setting up the nodal degrees of freedom
        for(unsigned int i = 0 ; i != number_of_nodes ; ++i) 
        {
            (GetGeometry()[i].pAddDof(DISPLACEMENT_X,REACTION_X));
            (GetGeometry()[i].pAddDof(DISPLACEMENT_Y,REACTION_Y));
            if(dimension == 3)
            {
                (GetGeometry()[i].pAddDof(DISPLACEMENT_Z,REACTION_Z));
            }
        }

        //initializing static variables
        unsigned int dof = number_of_nodes * dimension;
        msA.resize(dof, dof, false);
        msX.resize(dof, false);
        msU.resize(dof, false);
    }
    
    /**
    * Destructor.
    */
    LinearElasticTruss::~LinearElasticTruss()
    {
    }
        
    /**
    * Create a new Crisfield truss element.
    * @return crisfield truss elment
    * @param NewId element ID
    * @param ThisNodes array of nodes
    * @param pProperties properties pointer
    */
    Element::Pointer LinearElasticTruss::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new LinearElasticTruss(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /**
    * Initialization of the Crisfield truss element. 
    * This initializes the cross-section, length, position vector and matrix A for the element
    */
    void LinearElasticTruss::Initialize()
    {
        KRATOS_TRY
        unsigned int dof = number_of_nodes * dimension;  
        
        //cross-section mArea
        mArea = GetProperties()[THICKNESS];
        
        //length mLength (in ref. config)
        mLength = 0.0;
        mLength += (GetGeometry()[1].X0()-GetGeometry()[0].X0())*(GetGeometry()[1].X0()-GetGeometry()[0].X0());
        mLength += (GetGeometry()[1].Y0()-GetGeometry()[0].Y0())*(GetGeometry()[1].Y0()-GetGeometry()[0].Y0());
        if(dimension == 3)
            mLength += (GetGeometry()[1].Z0()-GetGeometry()[0].Z0())*(GetGeometry()[1].Z0()-GetGeometry()[0].Z0());
        mLength = sqrt( mLength );
//      KRATOS_WATCH(mLength);
                
        //position vector msX (in ref. config)
        for (unsigned int i = 0; i< number_of_nodes; i++) {
            int index= i*dimension;
            msX[index] = GetGeometry()[i].X0();
            msX[index+1] = GetGeometry()[i].Y0();
            if(dimension==3)
                msX[index+2] = GetGeometry()[i].Z0();
        }
//      KRATOS_WATCH( msX );    
        
        //matrix A
        noalias(msA) = ZeroMatrix(dof,dof);
        for (unsigned int i = 0 ; i < dof ; i++) {
            for(unsigned int j = 0 ; j < dof ; j++) {
                if(i==j)
                    msA(i,j) = 1;
                if(abs(static_cast<int>(i-j))==dimension)
                    msA(i,j) = -1;
            }
        }
        
        KRATOS_CATCH("")
    }

    /**
    * Calculation of the local system.
    * This calculates both the elemental stiffness matrix and the elemental residual vector
    * @param rLeftHandSideMatrix elemental stiffness matrix
    * @param rRightHandSideVector elemental residual vetor
    * @param rCurrentProcessInfo process info
    * @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
    * @see CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    */
    void LinearElasticTruss::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }
    
    /**
    * Calculation of the left hand side.
    * This calculates only the elemental stiffness matrix
    * @param rLeftHandSideMatrix elemental stiffness matrix
    * @param rCurrentProcessInfo process info
    * @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    * @see CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    */
    void LinearElasticTruss::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = false;
        VectorType temp = Vector();
        
        CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }
    
    /**
    * Calculation of the right hand side.
    * This calculates only the elemental residual vector
    * @param rRightHandSideVector elemental residual vetor
    * @param rCurrentProcessInfo process info
    * @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    * @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
    */
    void LinearElasticTruss::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();
        
        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /**
    * Get the equation ID vector of the element.
    * @param rResult equation ID vector
    * @param rCurrentProcessInfo process info
    */
    void LinearElasticTruss::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int dof = number_of_nodes*dimension;
        if(rResult.size() != dof)
            rResult.resize(dof);
        for (unsigned int i=0; i<number_of_nodes; i++) {
            int index = i*dimension;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            if(dimension == 3)
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    /**
    * Get the DOF list of the element.
    * @param ElementalDofList elemental DOF vector
    * @param rCurrentProcessInfo process info
    */
    void LinearElasticTruss::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        for (unsigned int i=0; i<number_of_nodes; i++) {
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            if(GetGeometry().WorkingSpaceDimension() == 3)
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    /**
    * Get the mass matrix of the element.
    * @param rMassMatrix mass matrix
    * @param rCurrentProcessInfo process info
    * 
    */
    void LinearElasticTruss::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
	const double& density = GetProperties()[DENSITY];
	//lumped
	unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int NumberOfNodes = GetGeometry().size();
	unsigned int MatSize = dimension * NumberOfNodes;
		
	if(rMassMatrix.size1() != MatSize)
		rMassMatrix.resize(MatSize,MatSize,false);

	noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);

 	double TotalMass = mLength * GetProperties()[DENSITY];

        if ( dimension == 2 ) 
		TotalMass *= GetProperties()[THICKNESS];

        if ( dimension == 3 ) 
		{
		double t_2=GetProperties()[THICKNESS]*GetProperties()[THICKNESS];
		TotalMass *= t_2;
		}

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

 	for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;
            }
        }
	//rMassMatrix*=0.0;
	//KRATOS_WATCH(rMassMatrix)

        KRATOS_CATCH("")
    }

    /**
    * Get the damping matrix of the element.
    * @param rDampMatrix damping matrix
    * @param rCurrentProcessInfo process info
    * TODO: assign the damping matrix
    */
    void LinearElasticTruss::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        int MatSize = number_of_nodes * dimension;
        rDampMatrix.resize(MatSize,MatSize);
        noalias(rDampMatrix)= ZeroMatrix(MatSize,MatSize);
        KRATOS_CATCH("")
    }
    
    /**
    * Get the displacement vector of the element
    * @param values displacement vector
    * @param Step solution step
    * @see GetFirstDerivativesVector(Vector& values, int Step)
    * @see GetSecondDerivativesVector(Vector& values, int Step)
    */
    void LinearElasticTruss::GetValuesVector(Vector& values, int Step)
    {
        unsigned int MatSize = number_of_nodes * dimension;
        if(values.size() != MatSize)    
            values.resize(MatSize);
        for ( unsigned int i=0; i<number_of_nodes; i++) {
            int index = i*dimension;
            values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X,Step);
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y,Step);
            if(dimension == 3)
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z,Step);
        }
        //KRATOS_WATCH( values );
    }

    /**
    * Get the velocity vector of the element
    * @param values velocity vector
    * @param Step solution step
    * @see GetValuesVector(Vector& values, int Step)
    * @see GetSecondDerivativesVector(Vector& values, int Step)
    */
    void LinearElasticTruss::GetFirstDerivativesVector(Vector& values, int Step)
    {
        unsigned int MatSize = number_of_nodes * dimension;
        if(values.size() != MatSize)    
            values.resize(MatSize);
        for (unsigned int i=0; i<number_of_nodes; i++) {
            int index = i*dimension;
            values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
            if(dimension == 3)
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
        }
    }
    
    /**
    * Get the acceleration vector of the element
    * @param values acceleration vector
    * @param Step solution step
    * @see GetValuesVector(Vector& values, int Step)
    * @see GetFirstDerivativesVector(Vector& values, int Step)
    */
    void LinearElasticTruss::GetSecondDerivativesVector(Vector& values, int Step)
    {
        unsigned int MatSize = number_of_nodes * dimension;
        if(values.size() != MatSize)    
            values.resize(MatSize);
        for (unsigned int i=0; i<number_of_nodes; i++) {
            int index = i*dimension;
            values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
            if(dimension == 3)
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
        }
    }
    
    /*
    void LinearElasticTruss::CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
    {
        if( rVariable == TRUSS_STRAIN )
        {
            if(Output.size() != 1)
                        Output.resize(1,false);
            Output[0] = msStrain;
        }
        if( rVariable == TRUSS_STRESS )
        {
            if(Output.size() != 1)
                        Output.resize(1,false);
            Output[0] = msStress;
        }
    }*/
    
    /**
    * Auxiliary function.
    * This calculates the elemental stiffness matrix and the elemental residual vector, when required
    * @param rLeftHandSideMatrix elemental stiffness matrix
    * @param rRightHandSideVector elemental residual vector
    * @param rCurrentProcessInfo process info
    * @param CalculateStiffnessMatrixFlag flag for elemental stiffness matrix
    * @param CalculateResidualVectorFlag flag for elemental residual vector
    */
    void LinearElasticTruss::CalculateAll(MatrixType& rLeftHandSideMatrix, 
                        VectorType& rRightHandSideVector, 
                        ProcessInfo& rCurrentProcessInfo,
                        bool CalculateStiffnessMatrixFlag,
                        bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        unsigned int MatSize = number_of_nodes * dimension;
        //KRATOS_WATCH("LINEAR TRUSSSSSSSSSSSSSSSSSSSS")
        //resizing the LHS tangent stiffness matrix if required
        if (CalculateStiffnessMatrixFlag == true) {
            if(rLeftHandSideMatrix.size1() != MatSize)
                rLeftHandSideMatrix.resize(MatSize,MatSize,false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
        }       

        //resizing the RHS residual vector if required
        if (CalculateResidualVectorFlag == true) {
            if(rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize);
            rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
        }       
        /*
	for (int i=0; i<GetGeometry().size();i++)
		{
		if (GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
			{
			GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z)=0.0;
			GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y)=0.0;
			GetGeometry()[i].Fix(DISPLACEMENT_Y);
			GetGeometry()[i].Fix(DISPLACEMENT_Z);
			}
	
		}
	*/
	double E = GetProperties()[YOUNG_MODULUS];
	double thickness = GetProperties()[THICKNESS];
	double A=thickness*thickness;

	Matrix Aux(MatSize,1);
	for (int i=0; i<MatSize; i++)
		{
		//Aux(i,0)=msX[i];
		Aux(0,0)=msX[0]-msX[3];
		Aux(1,0)=msX[1]-msX[4];
		Aux(2,0)=msX[2]-msX[5];

		Aux(3,0)=msX[3]-msX[0];
		Aux(4,0)=msX[4]-msX[1];
		Aux(5,0)=msX[5]-msX[2];
		}

	//to obtain the cos and sin of the angles we have to divide by L
	Aux/=mLength;

	noalias(rLeftHandSideMatrix)=(E*A/mLength)*prod(Aux,trans(Aux));

	array_1d<double,6> disp_old = ZeroVector(6);	

	disp_old[0] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
	disp_old[1] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
	disp_old[2] = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);

	disp_old[3] = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
	disp_old[4] = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
	disp_old[5] = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z);


       	CalculateAndAdd_ExtForce(rRightHandSideVector, rCurrentProcessInfo);
	noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, disp_old);
	//KRATOS_WATCH(rLeftHandSideMatrix)
	//KRATOS_WATCH(rRightHandSideVector)

        KRATOS_CATCH("")
    }
    
    /**
    * Auxiliary function.
    * This adds the external force vector to the elemental residual vector
    * @param rRightHandSideVector elemental residual vector
    * @param CurrentProcessInfo process info
    */
    void LinearElasticTruss::CalculateAndAdd_ExtForce(VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
	double l=mLength;
	double area=GetProperties()[THICKNESS]*GetProperties()[THICKNESS];
	double ext_force_factor = l*area;
	
        for (unsigned int i=0; i<number_of_nodes; i++) {
            int index = dimension*i;
            for ( int j=0; j<dimension; j++) {
                rRightHandSideVector[index+j] += 0.5*ext_force_factor*GetProperties()[BODY_FORCE][j]; 	
            }           
        }
        KRATOS_CATCH("")
    }
    
    /**
    * Auxiliary function.
    * This substracts the internal force vector to the elemental residual vector
    * @param rRightHandSideVector elemental residual vector
    * @param CurrentProcessInfo process info
    * @param X position vector
    * @param U displacement vector
    * @param weight weighting factor
    */
    /*
    void LinearElasticTruss::CalculateAndMinus_IntForce(VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo, Vector& X, Vector& U, double weight)
    {
        KRATOS_TRY
        for (int i=0; i<dimension; i++) {
            rRightHandSideVector[i] -= weight*( X[i]+U[i]-X[dimension+i]-U[dimension+i] );  
            rRightHandSideVector[dimension+i] -= rRightHandSideVector[i];
        }
        KRATOS_CATCH("")
    }
    */
    /**
    * Auxiliary function.
    * This add the material element stiffness matrix to the elemental stiffness matrix
    * @param rLeftHandSideMatrix elemental stiffness matrix
    * @param A matrix A
    * @param X position vector
    * @param U displacement vector
    * @param weight weighting factor
    */
    /*
    void LinearElasticTruss::CalculateAndAddKm(MatrixType& rLeftHandSideMatrix, const Matrix& A, Vector& X, Vector& U, double weight)
    {
        KRATOS_TRY
        Vector V(X.size());
        noalias(V) = prod(A,(X+U));
        for(unsigned int i=0; i<V.size(); i++)
        {
            for(unsigned int j=0; j<V.size(); j++)
            {
                rLeftHandSideMatrix(i,j)= rLeftHandSideMatrix(i,j) + weight * V[i] * V[j];
            }
        }
            
        KRATOS_CATCH("")
    }  
    */
    /**
    * Auxiliary function.
    * This add the geometrial element stiffness matrix to the elemental stiffness matrix
    * @param rLeftHandSideMatrix elemental stiffness matrix
    * @param A matrix A
    * @param weight weighting factor
    */
    /*
    void LinearElasticTruss::CalculateAndAddKg(MatrixType& rLeftHandSideMatrix, const Matrix& A, double weight)
    {
        KRATOS_TRY
        noalias(rLeftHandSideMatrix) += weight * A;
        KRATOS_CATCH("")
    } 
    */
    /**
    * Auxiliary function.
    * This calculates the GREEN-LAGRANGE strain.
    * @return GREEN-LAGRANGE strain
    * @param A matrix A
    * @param X position vector
    * @param U displacement vector
    * @param weight weighting factor
    */
    /*
    double LinearElasticTruss::CalculateStrain(const Matrix& A, const Vector& X, const Vector& U, double weight)
    {
        KRATOS_TRY
        Vector V = ZeroVector(X.size());
        noalias(V) = prod(A,U);
        double strain = 0;
        for(unsigned int i=0; i<X.size(); i++ ) {
            strain += (X[i]+0.5*U[i])*V[i];
        }
        return strain * weight ;
        KRATOS_CATCH("")
    } 
    */
    
} // Namespace Kratos

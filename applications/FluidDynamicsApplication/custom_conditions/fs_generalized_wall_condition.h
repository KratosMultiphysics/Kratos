/*
 ==============================================================================
 Kratos Fluid Dynamics Application
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 Version 1.0 (Released on march 05, 2007).

 Copyright 2007
 Pooyan Dadvand, Riccardo Rossi
 pooyan@cimne.upc.edu
 rrossi@cimne.upc.edu
 CIMNE (International Center for Numerical Methods in Engineering),
 Gran Capita' s/n, 08034 Barcelona, Spain

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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

#ifndef KRATOS_FS_GENERALIZED_WALL_CONDITION_H
#define KRATOS_FS_GENERALIZED_WALL_CONDITION_H

// System includes
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Implements a generalized wall model accounting for pressure gradients.
/**
 The generalized wall model in "T.-H. Shih, L.A. Povinelli, N.-S. Liu, M.G. Potapczuk
 and J.L. Lumley, A generalized wall function, Tech. Report NASA/TM-1999-209398, 1999"
 is used to calculate wall stress. For distributed problems the
 MetisDivideHeterogeneousInputProcess must be used with SynchronizeConditions = true.
 This ensures the condition belongs to the same partition as its parent element.

 Interface artificial compressibility (IAC) is used on the fluid-structure interface
 to improve the stability of partitioned FSI coupling iterations. This is activated
 by setting the flag INTERFACE to true.

 @see FractionalStep2D, FractionalStep3D
 */
template<unsigned int TDim, unsigned int TNumNodes = TDim>
class FSGeneralizedWallCondition: public Condition
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of FSGeneralizedWallCondition
	KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FSGeneralizedWallCondition);

	typedef Node < 3 > NodeType;

	typedef Properties PropertiesType;

	typedef Geometry<NodeType> GeometryType;

	typedef Geometry<NodeType>::PointType PointType;

	typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

	typedef Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

	typedef Element::WeakPointer ElementWeakPointerType;

	typedef Element::Pointer ElementPointerType;

	typedef std::size_t IndexType;

	typedef std::size_t SizeType;

	typedef std::vector<std::size_t> EquationIdVectorType;

	typedef std::vector< Dof<double>::Pointer > DofsVectorType;

	typedef Vector ShapeFunctionsType;

	/// Type for a matrix containing the shape function gradients
	typedef Matrix ShapeFunctionDerivativesType;

	/// Type for an array of shape function gradient matrices
	typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	/** Admits an Id as a parameter.
	 @param NewId Index of the new condition
	 */
	FSGeneralizedWallCondition(IndexType NewId = 0) : Condition(NewId), mInitializeWasPerformed(false), mpElement()
	{
	}

	/// Constructor using an array of nodes.
	/**
	 @param NewId Index of the new condition
	 @param ThisNodes An array containing the nodes of the new condition
	 */
	FSGeneralizedWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
	: Condition(NewId, ThisNodes), mInitializeWasPerformed(false), mpElement()
	{
	}

	/// Constructor using Geometry.
	/**
	 @param NewId Index of the new condition
	 @param pGeometry Pointer to a geometry object
	 */
	FSGeneralizedWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
	: Condition(NewId, pGeometry), mInitializeWasPerformed(false), mpElement()
	{
	}

	/// Constructor using Properties.
	/**
	 @param NewId Index of the new condition
	 @param pGeometry Pointer to a geometry object
	 @param pProperties Pointer to the condition's properties
	 */
	FSGeneralizedWallCondition(IndexType NewId, GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties)
	: Condition(NewId, pGeometry, pProperties), mInitializeWasPerformed(false), mpElement()
	{
	}

	/// Copy constructor.
	FSGeneralizedWallCondition(FSGeneralizedWallCondition const& rOther) : Condition(rOther),
	mInitializeWasPerformed(rOther.mInitializeWasPerformed), mpElement(rOther.mpElement)
	{
	}

	/// Destructor.
	~FSGeneralizedWallCondition() override
	{}

	///@}
	///@name Operators
	///@{

	/// Copy constructor.
	FSGeneralizedWallCondition& operator=(FSGeneralizedWallCondition const& rOther)
	{
		Condition::operator=(rOther);
		mInitializeWasPerformed = rOther.mInitializeWasPerformed;
		mpElement = rOther.mpElement;
		return *this;
	}

	///@}
	///@name Operations
	///@{

	/// Create a new FSGeneralizedWallCondition object.
	/**
	 @param NewId Index of the new condition
	 @param ThisNodes An array containing the nodes of the new condition
	 @param pProperties Pointer to the condition's properties
	 */
    Condition::Pointer Create(
		IndexType NewId,
		NodesArrayType const& ThisNodes,
		PropertiesType::Pointer pProperties) const override
	{
        return Kratos::make_intrusive<FSGeneralizedWallCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
	}

	/// Create a new FSGeneralizedWallCondition object.
	/**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
		IndexType NewId,
		GeometryType::Pointer pGeom,
		PropertiesType::Pointer pProperties) const override
	{
		return Kratos::make_intrusive<FSGeneralizedWallCondition>(NewId, pGeom, pProperties);
    }

        /// Find the condition's parent element.
	void Initialize() override
	{
		KRATOS_TRY;

		if (this->Is(SLIP))
		{
			const array_1d<double, 3> &rNormal = this->GetValue(NORMAL);
			KRATOS_ERROR_IF(norm_2(rNormal) == 0.0) << "NORMAL must be calculated before using this " << this->Info() << "\n";
		}

		if (mInitializeWasPerformed)
		{
			return;
		}

		mInitializeWasPerformed = true;

		KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 0) << this->Info() << " cannot find parent element\n";

		double EdgeLength;
		array_1d<double, 3> Edge;
		mpElement = this->GetValue(NEIGHBOUR_ELEMENTS)(0);
		GeometryType &rElemGeom = mpElement->GetGeometry();

		Edge = rElemGeom[1].Coordinates() - rElemGeom[0].Coordinates();
		mMinEdgeLength = Edge[0] * Edge[0];
		for (SizeType d = 1; d < TDim; d++)
		{
			mMinEdgeLength += Edge[d] * Edge[d];
		}

		for (SizeType j = 2; j < rElemGeom.PointsNumber(); j++)
		{
			for (SizeType k = 0; k < j; k++)
			{
				Edge = rElemGeom[j].Coordinates() - rElemGeom[k].Coordinates();
				EdgeLength = Edge[0] * Edge[0];

				for (SizeType d = 1; d < TDim; d++)
				{
					EdgeLength += Edge[d] * Edge[d];
				}

				mMinEdgeLength = (EdgeLength < mMinEdgeLength) ? EdgeLength : mMinEdgeLength;
			}
		}
		mMinEdgeLength = sqrt(mMinEdgeLength);
		return;

		KRATOS_CATCH("");
	}

	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
			ProcessInfo& rCurrentProcessInfo) override
	{
		VectorType RHS;
		this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
	}

	/// Calculate wall stress term for all nodes with SLIP set.
	/**
	 @param rLeftHandSideMatrix Left-hand side matrix
	 @param rRightHandSideVector Right-hand side vector
	 @param rCurrentProcessInfo ProcessInfo instance
	 */
	void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override
	{
		KRATOS_TRY;

		if (mInitializeWasPerformed == false)
		{
		        Initialize();
		}

		if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
		{
			// Initialize local contributions
			const SizeType LocalSize = TDim * TNumNodes;

			if (rLeftHandSideMatrix.size1() != LocalSize)
			{
				rLeftHandSideMatrix.resize(LocalSize, LocalSize);
			}
			if (rRightHandSideVector.size() != LocalSize)
			{
				rRightHandSideVector.resize(LocalSize);
			}

			noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
			noalias(rRightHandSideVector) = ZeroVector(LocalSize);

			if (this->Is(SLIP))
			  this->ApplyWallLaw(rLeftHandSideMatrix, rRightHandSideVector);
		}
		else if (rCurrentProcessInfo[FRACTIONAL_STEP] == 5)
		{
			// add IAC penalty to local pressure system
			const SizeType LocalSize = TNumNodes;

			if (rLeftHandSideMatrix.size1() != LocalSize)
			{
				rLeftHandSideMatrix.resize(LocalSize, LocalSize);
			}
			if (rRightHandSideVector.size() != LocalSize)
			{
				rRightHandSideVector.resize(LocalSize);
			}

			noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
			noalias(rRightHandSideVector) = ZeroVector(LocalSize);

			if (this->Is(INTERFACE))
			  this->ApplyIACPenalty(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
		}
		else
		{
			if (rLeftHandSideMatrix.size1() != 0)
			{
				rLeftHandSideMatrix.resize(0,0,false);
			}

			if (rRightHandSideVector.size() != 0)
			{
				rRightHandSideVector.resize(0,false);
			}
		}

		KRATOS_CATCH("");
	}

	/// Check that all data required by this condition is available and reasonable.
	int Check(const ProcessInfo& rCurrentProcessInfo) override
	{
		KRATOS_TRY;

		int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area >= 0

		if (Check != 0)
		{
			return Check;
		}
		else
		{
			// Check that all required variables have been registered
			if(VELOCITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
			if(PRESSURE.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
			if(MESH_VELOCITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
			if(DENSITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
			if(VISCOSITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","");
			if(NORMAL.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"NORMAL Key is 0. Check if the application was correctly registered.","");

			// Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
			for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
			{
				if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
						this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
						this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
				if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
				KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE degree of freedom on node ",this->GetGeometry()[i].Id());
			}

			return Check;
		}

		KRATOS_CATCH("");
	}

	/// Provides the global indices for each one of this condition's local rows.
	/** This determines the equation ID vector for all DOFs.
	 * @param rResult A vector containing the global Id of each row
	 * @param rCurrentProcessInfo the current process info object
	 */
	void EquationIdVector(EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo) override;

	/// Returns a list of the condition's Dofs.
	/**
	 * @param ConditionDofList the list of DOFs
	 * @param rCurrentProcessInfo the current process info instance
	 */
	void GetDofList(DofsVectorType& rConditionDofList,
			ProcessInfo& rCurrentProcessInfo) override;

	/// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z) for each node.
	/**
	 * @param Values Vector of nodal unknowns
	 * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
	 */
	void GetValuesVector(Vector& Values, int Step = 0) override
	{
		const SizeType LocalSize = TDim * TNumNodes;
		unsigned int LocalIndex = 0;

		if (Values.size() != LocalSize)
		{
			Values.resize(LocalSize, false);
		}

		for (unsigned int i = 0; i < TNumNodes; ++i)
		{
			array_1d<double,3>& rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
			for (unsigned int d = 0; d < TDim; ++d)
			{
				Values[LocalIndex++] = rVelocity[d];
			}
		}
	}

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const override
	{
		std::stringstream buffer;
		buffer << "FSGeneralizedWallCondition" << TDim << "D";
		return buffer.str();
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override
	{	rOStream << "FSGeneralizedWallCondition";}

	/// Print object's data.
	void PrintData(std::ostream& rOStream) const override
	{}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	ElementPointerType pGetElement()
	{
		return mpElement->shared_from_this();
	}

	template< class TVariableType >
	void EvaluateInPoint(TVariableType& rResult,
			const Kratos::Variable<TVariableType>& Var,
			const ShapeFunctionsType& rShapeFunc)
	{
		GeometryType& rGeom = this->GetGeometry();

		rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

		for(SizeType i = 1; i < TNumNodes; i++)
		{
			rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
		}
	}

	/// Calculate the pressure gradient using nodal data from the last time step.
	void EvaluateOldPressureGradientInElement(array_1d<double,TDim>& rResult)
	{
	  GeometryType& rElemGeom = pGetElement()->GetGeometry();
	  const SizeType NumNodes = rElemGeom.PointsNumber();
	  ShapeFunctionDerivativesArrayType DN_DX;
	  Vector DetJ;
	  rElemGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ,
		          GeometryData::GI_GAUSS_1);
	  ShapeFunctionDerivativesType& rDN_DX = DN_DX[0];

	  const double& pres = rElemGeom[0].FastGetSolutionStepValue(PRESSURE,1);
	  for (SizeType d = 0; d < TDim; ++d)
	    {
	      rResult[d] = rDN_DX(0,d) * pres;
	    }

	  for (SizeType i = 1; i < NumNodes; i++)
	    {
	      const double& pres = rElemGeom[i].FastGetSolutionStepValue(PRESSURE,1);
	      for (SizeType d = 0; d < TDim; ++d)
		{
		  rResult[d] += rDN_DX(i,d) * pres;
		}
	    }
	}

	/// Calculate input parameters to wall model.
	/**
	 * @param rWallHeight Height of the measurement point above the wall
	 * @param rWallVel Tangential velocity vector at the measurement point
	 * @param rWallGradP Wall pressure gradient in the direction of rWallVel
	 * @param rArea Condition's area
	 */
	void CalculateWallParameters(double& rWallHeight, array_1d<double,3>& rWallVel, double& rWallGradP, double& rArea);

	/// Evaluate the residual of the generalized wall function.
	/**
	 * @param rWallHeight Height of the measurement point above the wall
	 * @param rWallVel Tangential velocity at the measurement point
	 * @param rWallStress Wall stress in direction of rWallVel
	 * @param rWallGradP Wall pressure gradient in the direction of rWallVel
	 */
	double EvaluateWallFunctionResidual(const double& rWallHeight, const double& rWallVel, const double& rWallStress, const double& rWallGradP)
	{
		const ShapeFunctionsType& N = row(this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1),0);
		double rho, nu, func1, func2, Vel1, Vel2, Vel12, YPlus1, YPlus2, sign1, sign2;
		EvaluateInPoint(rho, DENSITY, N);
		EvaluateInPoint(nu, VISCOSITY, N);
		Vel1 = sqrt(fabs(rWallStress) / rho);
		Vel2 = pow(nu * fabs(rWallGradP) / rho, 0.333333);
		Vel12 = Vel1 + Vel2;

		if (Vel12 == 0.0)
		{
			Vel12 = 1.0;
		}

		YPlus1 = Vel1 * rWallHeight / nu;
		YPlus2 = Vel2 * rWallHeight / nu;

		// evaluate func1(YPlus1)
		func1 = 0.0;
		if (YPlus1 <= 5)
		{
			func1 = (1.0 + (0.01 - 0.0029 * YPlus1) * YPlus1) * YPlus1;
		}
		else if (YPlus1 <= 30)
		{
			func1 = -0.872 + (1.465 + (-0.0702 + (0.00166 - 0.00001495 * YPlus1) * YPlus1) * YPlus1) * YPlus1;
		}
		else if (YPlus1 <= 140)
		{
			func1 = 8.6 + (0.1864 + (-0.002006 + (0.00001144 - 0.00000002551 * YPlus1) * YPlus1) * YPlus1) * YPlus1;
		}
		else
		{
			func1 = 2.439 * log(YPlus1) + 5.0;
		}

		// evaluate func2(YPlus2)
		func2 = 0.0;
		if (YPlus2 <= 4)
		{
			func2 = (0.5 - 0.00731 * YPlus2) * YPlus2 * YPlus2;
		}
		else if (YPlus2 <= 15)
		{
			func2 = -15.138 + (8.4688 + (-0.81976 + (0.037292 - 0.00063866 * YPlus2) * YPlus2) * YPlus2) * YPlus2;
		}
		else if (YPlus2 <= 30)
		{
			func2 = 11.925 + (0.934 + (-0.027805 + (0.00046262 - 0.0000031442 * YPlus2) * YPlus2) * YPlus2) * YPlus2;
		}
		else
		{
			func2 = 5.0 * log(YPlus2) + 8.0;
		}

		sign1 = (rWallStress >= 0.0) ? 1.0 : -1.0;
		sign2 = (rWallGradP >= 0.0) ? 1.0 : -1.0;
		return (rWallVel - sign1 * Vel1 * func1 - sign2 * Vel2 * func2) / Vel12;
	}

	/// Calculate the wall stress from the wall function
	/**
	 * @param rWallHeight Height of the measurement point above the wall
	 * @param rWallVel Tangential velocity at the measurement point
	 * @param rWallGradP Wall pressure gradient in the direction of rWallVel
	 */
	double CalculateWallStress(const double& rWallHeight, const double& rWallVel, const double& rWallGradP)
	{
		KRATOS_TRY;

		const SizeType MaxIter = 50;
		const double Eps = 1.0e-08;
		double A =-0.1;
		double B = 0.1;
		double C, Xm, P, Q, R, S, Min1, Min2, ResA, ResB, ResC, Tol;
                double D=0.0;
                double E=0.0;

		for (SizeType k=0; k < 10; ++k)
		{
			A *= 10.0;
			B *= 10.0;
			ResA = EvaluateWallFunctionResidual(rWallHeight,rWallVel,A,rWallGradP);
			ResB = EvaluateWallFunctionResidual(rWallHeight,rWallVel,B,rWallGradP);

			if (ResB * ResA <= 0.0)
			{
				break;
			}
		}

		if (ResB * ResA > 0.0)
		{
		  return 0.;
		  //KRATOS_THROW_ERROR(std::logic_error, "Maximum wall stress search width exceeded","");
		}

		// root finding with Brent's algorithm
		C = B;
		ResC = ResB;

		for (SizeType i=1; i <= MaxIter; ++i)
		{
			if (ResB * ResC > 0.0)
			{
				C = A;
				ResC = ResA;
				D = B - A;
				E = D;
			}

			if (fabs(ResC) < fabs(ResB))
			{
				A = B;
				B = C;
				C = A;
				ResA = ResB;
				ResB = ResC;
				ResC = ResA;
			}

			Tol = 2.0 * Eps * fabs(B) + 0.5e-10;
			Xm = 0.5 * (C - B);

			if (fabs(Xm) <= Tol || ResB == 0.0)
			{
				return B;
			}

			if (fabs(E) >= Tol && fabs(ResA) > fabs(ResB))
			{
				S = ResB / ResA;

				if (A == C)
				{
					P = 2.0 * Xm * S;
					Q = 1.0 - S;
				}
				else
				{
					Q = ResA / ResC;
					R = ResB / ResC;
					P = S * (2.0 * Xm * Q * (Q - R) - (B - A) * (R - 1.0));
					Q = (Q - 1.0) * (R - 1.0) * (S - 1.0);
				}

				if (P > 0.0)
				{
					Q = -Q;
				}

				P = fabs(P);
				Min1 = 3.0 * Xm * Q - fabs(Tol * Q);
				Min2 = fabs(E * Q);

				if (2.0 * P < std::min(Min1, Min2))
				{
					E = D;
					D = P / Q;
				}
				else
				{
					D = Xm;
					E = D;
				}
			}
			else
			{
				D = Xm;
				E = D;
			}

			A = B;
			ResA = ResB;

			if (fabs(D) > Tol)
			{
				B += D;
			}
			else
			{
				if (Xm >= 0.0)
				{
					B += fabs(Tol);
				}
				else
				{
					B -= fabs(Tol);
				}
			}

			ResB = EvaluateWallFunctionResidual(rWallHeight,rWallVel,B,rWallGradP);
		}

		return 0.0;

		KRATOS_CATCH("");
	}

	/// Check angles between condition normal and node normals
	/**
	@param AngleCosine Cosine of minimum angle between two normals at a corner
	 */
	bool IsCorner(double AngleCosine)
	{
		GeometryType& rGeometry = this->GetGeometry();
		const array_1d<double, 3>& rNormal = this->GetValue(NORMAL);
		const double Area = norm_2(rNormal);

		for (SizeType iNode = 0; iNode < rGeometry.PointsNumber(); ++iNode)
		{
			const array_1d<double, 3>& rNodeNormal = rGeometry[iNode].FastGetSolutionStepValue(NORMAL);
			if (inner_prod(rNormal, rNodeNormal) < AngleCosine * Area * norm_2(rNodeNormal))
				return true;
		}
		return false;
	}

	/// Compute the wall stress and add corresponding terms to the system contributions.
	/**
	 @param rLocalMatrix Local system matrix
	 @param rLocalVector Local right hand side
	 */
	void ApplyWallLaw(MatrixType& rLocalMatrix, VectorType& rLocalVector)
	{
		const unsigned int BlockSize = TDim;
		double WallHeight, Area, WallVelMag, tmp, WallStress, WallGradP, WallForce;
		array_1d<double,3> WallVel;
		GeometryType& rGeometry = this->GetGeometry();

		CalculateWallParameters(WallHeight, WallVel, WallGradP, Area);
		WallVelMag = norm_2(WallVel);

		if (IsCorner(0.966) == false)
		{
			WallStress = CalculateWallStress(WallHeight, WallVelMag, WallGradP);
			WallForce = (Area / static_cast<double>(TDim)) * WallStress;

			for(SizeType i=0; i < rGeometry.PointsNumber(); ++i)
			{
				const NodeType& rNode = rGeometry[i];
				if(rNode.GetValue(Y_WALL) != 0.0 && rNode.Is(SLIP))
				{
					WallVel = rNode.FastGetSolutionStepValue(VELOCITY,1) - rNode.FastGetSolutionStepValue(MESH_VELOCITY,1);
					tmp = norm_2(WallVel);
					WallVel /= (tmp != 0.0) ? tmp : 1.0;

					for (unsigned int d=0; d < TDim; d++)
					{
						unsigned int k = i * BlockSize + d;
						rLocalVector[k] -= WallVel[d] * WallForce;
					}
				}
			}
		}
	}

	/// Apply an IAC penalty term
	/**
	 @param rLeftHandSideMatrix Left-hand side matrix
	 @param rRightHandSideVector Right-hand side vector
	 @param rCurrentProcessInfo ProcessInfo instance
	 */
	void ApplyIACPenalty(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo)
	{
		GeometryType& rGeometry = this->GetGeometry();
		const array_1d<double,3>& rNormal = this->GetValue(NORMAL);
		const double Area = norm_2(rNormal);
		const double DiagonalTerm = Area / static_cast<double>(TNumNodes)
		/ (rCurrentProcessInfo[DENSITY] * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

		for(SizeType iNode=0; iNode < rGeometry.PointsNumber(); ++iNode)
		{
			rLeftHandSideMatrix(iNode,iNode) += DiagonalTerm;
		}
	}

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{
	bool mInitializeWasPerformed;
	double mMinEdgeLength;
	ElementWeakPointerType mpElement;

	///@}
	///@name Serialization
	///@{

	friend class Serializer;

	void save(Serializer& rSerializer) const override
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
	}

	void load(Serializer& rSerializer) override
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
	}

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	///@}

}; // Class FSGeneralizedWallCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator >>(std::istream& rIStream,
		FSGeneralizedWallCondition<TDim, TNumNodes>& rThis)
{
	return rIStream;
}

/// output stream function
template<unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator <<(std::ostream& rOStream,
		const FSGeneralizedWallCondition<TDim, TNumNodes>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_FS_GENERALIZED_WALL_CONDITION_H

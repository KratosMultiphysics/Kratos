//   
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_2d_modified.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

	namespace
	{
    	template <class TMatrixType>
    	void FillConsistentMassMatrix(TMatrixType& rMassMatrix, const double Area)
    	{
        	const double factor = Area / 12.0;
        	for (unsigned int i = 0; i < 3; ++i)
        	{
            	for (unsigned int j = 0; j < 3; ++j)
            	{
               	 rMassMatrix(i, j) = factor;
            	}
            	rMassMatrix(i, i) += factor;
        	}
    	}
	}

	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM22DModified::MonolithicPFEM22DModified(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{}

	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM22DModified::MonolithicPFEM22DModified(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{}

	Element::Pointer MonolithicPFEM22DModified::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new MonolithicPFEM22DModified(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

    Element::Pointer MonolithicPFEM22DModified::Create(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new MonolithicPFEM22DModified(NewId, pGeometry, pProperties));
    }

	MonolithicPFEM22DModified::~MonolithicPFEM22DModified()
	{
	}


	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const int TDim = 2;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = TDim+1;
		const GeometryType& rGeom = this->GetGeometry();

		if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

		KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(DELTA_TIME))
            << "DELTA_TIME is not set in ProcessInfo." << std::endl;
        const double dt = rCurrentProcessInfo[DELTA_TIME];
        KRATOS_ERROR_IF(dt <= 0.0)
            << "DELTA_TIME must be positive, got " << dt << std::endl;
        const double inv_dt = 1.0 / dt;

        // theta = 1.0 -> BDF1 (Backward Euler), 0.5 -> Crank-Nicolson, 0.0 -> Forward Euler
        const double theta = rCurrentProcessInfo.Has(TIMEDIS_THETA)
            ? rCurrentProcessInfo[TIMEDIS_THETA]
            : 1.0;

        const double diffusivity = this->GetProperties().Has(CONDUCTIVITY)
            ? this->GetProperties()[CONDUCTIVITY]
            : 1.0;

		double Area;
        BoundedMatrix<double, 3, 2> DN_DX;
        array_1d<double, 3> N;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

		// Consistent mass matrix M_ij = (Area / 12) * (1 + delta_ij)
        BoundedMatrix<double, 3, 3> M;
        FillConsistentMassMatrix(M, Area);

		// Diffusion (stiffness) matrix K_ij = Area * k * (grad N_i . grad N_j)
        const BoundedMatrix<double, 3, 3> K = Area * diffusivity * prod(DN_DX, trans(DN_DX));

        array_1d<double, 3> T_curr, T_prev;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            const auto& rNode = this->GetGeometry()[i];
            T_curr[i] = rNode.FastGetSolutionStepValue(TEMPERATURE, 0);
            T_prev[i] = (rNode.GetBufferSize() > 1)
                ? rNode.FastGetSolutionStepValue(TEMPERATURE, 1)
                : T_curr[i];
        }

        array_1d<double, 3> f_curr, f_prev;
        this->GetNodalSource(f_curr, 0);
        this->GetNodalSource(f_prev, 1);

        const array_1d<double, 3> f_theta = theta * f_curr + (1.0 - theta) * f_prev;
        const array_1d<double, 3> T_theta = theta * T_curr + (1.0 - theta) * T_prev;

        // Theta-method residual form: (1/dt * M + theta * K) * dT = M * f_theta - (1/dt) * M * (T_curr - T_prev) - K * T_theta
        noalias(rLeftHandSideMatrix) = inv_dt * M + theta * K;
        noalias(rRightHandSideVector) = prod(M, f_theta)
            - inv_dt * prod(M, T_curr - T_prev)
            - prod(K, T_theta);

		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the
	// fractional step procedure
	void MonolithicPFEM22DModified::InitializeSolutionStep(const ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
	{
		const int TDim = 2;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = TDim+1;
		const GeometryType& rGeom = this->GetGeometry();

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const
	{           
		const int TDim = 2;
		
		const SizeType NumNodes = TDim+1;
        const SizeType LocalSize = TDim+1;
        const GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
            ElementalDofList.resize(LocalSize);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            ElementalDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
        }

	}


	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::GetNodalSource(array_1d<double, 3>& rSource, const IndexType StepIndex) const
    {
        const GeometryType& rGeom = this->GetGeometry();
        rSource = ZeroVector(3);

        bool found_nodal = false;
        for (unsigned int i = 0; i < 3; ++i)
        {
            const IndexType step = (rGeom[i].GetBufferSize() > StepIndex) ? StepIndex : 0;
            if (rGeom[i].SolutionStepsDataHas(HEAT_FLUX))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(HEAT_FLUX, step);
                found_nodal = true;
            }
            else if (rGeom[i].SolutionStepsDataHas(BODY_FORCE))
            {
                rSource[i] = rGeom[i].FastGetSolutionStepValue(BODY_FORCE, step)[0];
                found_nodal = true;
            }
        }

        // If no nodal source variable is present, fall back to property-level HEAT_FLUX
        if (!found_nodal)
        {
            double prop_source = 0.0;
            if (this->GetProperties().Has(HEAT_FLUX))
                prop_source = this->GetProperties()[HEAT_FLUX];

            for (unsigned int i = 0; i < 3; ++i)
                rSource[i] = prop_source;
        }
    }


	void MonolithicPFEM22DModified::AddViscousTerm(BoundedMatrix<double, 13, 13 > & output,
										  BoundedMatrix<double, (3), 2 >& rShapeDeriv,
										  array_1d<double,3>&  distances,
										  std::vector< Matrix >& gauss_gradients,
										  array_1d<double,3>&  viscosities,
										  array_1d<double,3>&  signs,
										  array_1d<double,3>&  volumes ,
										  const unsigned int ndivisions)
	{
		BoundedMatrix<double, 9, 9 > ExtendedDampMatrix=ZeroMatrix(9,9);
		//BoundedMatrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		BoundedMatrix<double, 9,3 > B_matrix = ZeroMatrix(9,3);

		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			//standard shape functions:
			for (unsigned int i=0; i!=(3); i++) //i node
			{
				for (unsigned int j=0; j!=(2); j++) //x,y,z
					B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);

				//useful for both 2d and 3d:
				//relating 12 and 21 stresses
				B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
				B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);

			}

			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(3*(2)+j,j)= gauss_gradients[division](0,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(3*(2)+0,2)=gauss_gradients[division](0,1);
			B_matrix(3*(2)+1,2)=gauss_gradients[division](0,0);

			//the jump in x direction
			B_matrix(4*(2)+0,0)= gauss_gradients[division](1,0);
			//relating 12 and 21 stresses
			B_matrix(4*(2)+0,2)=gauss_gradients[division](1,1);

			BoundedMatrix<double, 3 , 9  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}

		//now we put it all together in the big matrix:
		for (unsigned int i=0; i!=(4); i++) //4 nodes + 1dof in the new virtual node
			for (unsigned int j=0; j!=(4); j++) //4 nodes + 1dof in the new virtual node
				for (unsigned int k=0; k!=(2); k++) //x,y,(z)
					for (unsigned int l=0; l!=(2); l++) //x,y,(z)
						 output(i*3+k,j*3+l) += ExtendedDampMatrix(i*2+k,j*2+l);

		output(11,11) += ExtendedDampMatrix(8,8);
		for (unsigned int i=0; i!=(4); i++) // 3 nodes + 1 new for gradient.
			for (unsigned int k=0; k!=(2); k++) //x,y,(z)
			{
				output(i*3+k,11) += ExtendedDampMatrix(i*2+k,8);
				output(11,i*3+k) += ExtendedDampMatrix(8,i*2+k);
			}
	}


    void MonolithicPFEM22DModified::AddViscousTerm(BoundedMatrix<double, 12, 12 > & output,
                                        BoundedMatrix<double, (3), 2 >& rShapeDeriv,
                                        array_1d<double,3>&  distances,
                                        std::vector< Matrix >& gauss_gradients,
                                        array_1d<double,3>&  viscosities,
                                        array_1d<double,3>&  signs,
                                        array_1d<double,3>&  volumes ,
                                        const unsigned int ndivisions)
	{
		BoundedMatrix<double, 8, 8 > ExtendedDampMatrix=ZeroMatrix(8,8);
		//BoundedMatrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		BoundedMatrix<double, 8,3 > B_matrix = ZeroMatrix(8,3);

		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			//standard shape functions:
			for (unsigned int i=0; i!=(3); i++) //i node
			{
				for (unsigned int j=0; j!=(2); j++) //x,y,z
					B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);

				//useful for both 2d and 3d:
				//relating 12 and 21 stresses
				B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
				B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);

			}

			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(3*(2)+j,j)= gauss_gradients[division](0,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(3*(2)+0,2)=gauss_gradients[division](0,1);
			B_matrix(3*(2)+1,2)=gauss_gradients[division](0,0);


			BoundedMatrix<double, 3 , 8  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}

		//now we put it all together in the big matrix:
		for (unsigned int i=0; i!=(4); i++) //4 nodes + 1dof in the new virtual node
			for (unsigned int j=0; j!=(4); j++) //4 nodes + 1dof in the new virtual node
				for (unsigned int k=0; k!=(2); k++) //x,y,(z)
					for (unsigned int l=0; l!=(2); l++) //x,y,(z)
						 output(i*3+k,j*3+l) += ExtendedDampMatrix(i*2+k,j*2+l);
	}


	//with matrixtype, constant coefficient
	void MonolithicPFEM22DModified::AddViscousTerm(MatrixType& rDampMatrix,
                         const BoundedMatrix<double, 3, 2 >& rShapeDeriv,
                         double& Viscosity,const double Area)
	{
		BoundedMatrix<double, 6,3 > B_matrix = ZeroMatrix(6,3);
		for (unsigned int i=0; i!=(3); i++) //i node
		{
			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
			B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);
		}

		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		C_matrix*= Viscosity*Area;

		BoundedMatrix<double, 3 , 6  > temp_matrix = prod(C_matrix,trans(B_matrix));
		BoundedMatrix<double, 6 , 6  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=3; i++) //i node
		{
			for (unsigned int j=0; j!=3; j++) //j neighbour
			{
				for (unsigned int k=0; k!=2; k++) //xyz
					for (unsigned int l=0; l!=2; l++) //xyz
						rDampMatrix(i*3+k,j*3+l)+=viscosity_matrix(i*2+k,j*2+l);
			}
		}
	}


	template<class T>
	bool MonolithicPFEM22DModified::InvertMatrix(const T& input, T& inverse)
	{
		typedef permutation_matrix<std::size_t> pmatrix;

		// create a working copy of the input
		T A(input);

		// create a permutation matrix for the LU-factorization
		pmatrix pm(A.size1());

		// perform LU-factorization
		int res = lu_factorize(A, pm);
		if (res != 0)
			return false;

		// create identity matrix of "inverse"
		inverse.assign(identity_matrix<double> (A.size1()));

		// backsubstitute to get the inverse
		lu_substitute(A, pm, inverse);

		return true;
	}

} // Namespace Kratos
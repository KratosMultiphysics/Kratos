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


	void MonolithicPFEM22DModified::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
			case 10: //calculating pressure projection. notthing returned. saving data in PRESS_PROJ, PRESS_PROJ_NO_RO , NODAL_MASS and NODAL_AREA
			{
				this->CalculatePressureProjection(rCurrentProcessInfo);
				break;
			}
			default:
			{
				KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
		}

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_THROW_ERROR(std::logic_error, "MonolithicPFEM22DModified::CalculateLocalSystem has yet to be implemented", "");

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
		if(CurrentProcessInfo[FRACTIONAL_STEP]<100 )
		{
		const int TDim = 2;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = (TDim+1)*(TDim+1);
		const GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE).EquationId();
		}
		}
		else
		{
			const int TDim = 2;

			const SizeType NumNodes = TDim+1;
			const SizeType LocalSize = (TDim+1);
			const GeometryType& rGeom = this->GetGeometry();

			SizeType LocalIndex = 0;

			if (rResult.size() != LocalSize)
				rResult.resize(LocalSize, false);

			for (SizeType i = 0; i < NumNodes; ++i)
			{
				if(CurrentProcessInfo[FRACTIONAL_STEP]==100)
					rResult[LocalIndex++] = rGeom[i].GetDof(DISTANCE).EquationId();
				if(CurrentProcessInfo[FRACTIONAL_STEP]==101)
					rResult[LocalIndex++] = rGeom[i].GetDof(PROJECTED_VELOCITY_X).EquationId();
				if(CurrentProcessInfo[FRACTIONAL_STEP]==102)
					rResult[LocalIndex++] = rGeom[i].GetDof(PROJECTED_VELOCITY_Y).EquationId();
			}
		}
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& CurrentProcessInfo) const
	{                                                                                                                                                                                                                   
		if(CurrentProcessInfo[FRACTIONAL_STEP]<100 )
		{
		const int TDim = 2;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = (TDim+1)*(TDim+1);
		const GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
		}
		}
		else
		{
			const int TDim = 2;

			const SizeType NumNodes = TDim+1;
			const SizeType LocalSize = (TDim+1);
			const GeometryType& rGeom = this->GetGeometry();

			if (ElementalDofList.size() != LocalSize)
				ElementalDofList.resize(LocalSize);

			SizeType LocalIndex = 0;

			for (SizeType i = 0; i < NumNodes; ++i)
			{
				if(CurrentProcessInfo[FRACTIONAL_STEP]==100)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(DISTANCE);
				if(CurrentProcessInfo[FRACTIONAL_STEP]==101)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PROJECTED_VELOCITY_X);
				if(CurrentProcessInfo[FRACTIONAL_STEP]==102)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PROJECTED_VELOCITY_Y);
			}
		}


	}


	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM22DModified::CalculatePressureProjection(const ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int TDim=2;

		double Area;
		Geometry<Node >& geom = this->GetGeometry();
		BoundedMatrix<double, (TDim+1), TDim > DN_DX;
		array_1d<double, (TDim+1) > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double factor = 1.0/ (1.0 + double (TDim) );

		//const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);

		//if( (number_of_particles_in_elem>0))
		if(true)
		{

            const double density_air = CurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
            const double density_water = CurrentProcessInfo[DENSITY_WATER];

            array_1d<double,TDim+1>  pressures = ZeroVector(TDim+1); //to calculate the deformation Gradient F. Dimension = velocity dofs
            BoundedMatrix<double,TDim+1, TDim > coords; //coordinates of the nodes
            bool has_negative_node=false;
            bool has_positive_node=false;
            double element_mean_distance=0.0;
            for(unsigned int iii = 0; iii<(TDim+1); iii++)
            {
                //saving everything
                pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

                //
                const array_1d<double, 3 > & xyz = this->GetGeometry()[iii].Coordinates();
                for (unsigned int j = 0; j < TDim; j++)
                    coords(iii, j) = xyz[j];
                //
                if (this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)<0.0)
                    has_negative_node=true;
                else
                    has_positive_node=true;

                element_mean_distance+=factor*this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
            }

            bool split_element=false;
            if (has_positive_node && has_negative_node)
                split_element=true;

            if (split_element==false) //we only calculate if we have a pure fluid element
            {


                BoundedMatrix<double, (TDim+1), (TDim+1)*TDim > G_matrix; //(gradient)
                noalias(G_matrix) = ZeroMatrix((TDim+1), (TDim+1)*TDim);

                //const double water_fraction = -( element_mean_distance - 1.0) *0.5;
                //const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
                //double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
                double density = density_air;
                if (has_negative_node==true)
                    density= density_water;

                for (unsigned int i = 0; i < (TDim+1); i++) //loop in shape functions (row)
                {
                    for (unsigned int j = 0; j < (TDim+1) ; j++) //loop in shape functions (column)
                    {
                        for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
                        {
                            G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
                        }
                    }
                }
                G_matrix /=density;

                //now we save the data:
                for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
                {
                    geom[i].SetLock();
                    array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);

                    for (unsigned int j = 0; j < (TDim+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
                    {

                        for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
                        {
                            current_press_proj[k] += G_matrix(j,i*(TDim)+k)*(pressures(j));///-old_pressures(j)); //gamma=0!
                        }
                    }

                    geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor;

                    geom[i].UnSetLock();

                }

            } //closing the useful elements
            else// we add some small addition to the area so that we do not have a division by zero.
            {
                for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
                {
                    geom[i].SetLock();
                    geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor*0.000000000000001;
                    geom[i].UnSetLock();
                }
            }

		} //closing the if(is_inactive==false)
		else// we add some small addition to the area so that we do not have a division by zero.
		{
			for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
			{
				geom[i].SetLock();
				geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor*0.000000000000001;
				geom[i].UnSetLock();
			}
		}

		KRATOS_CATCH("");
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

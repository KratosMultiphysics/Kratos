// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//                   
//                   
//
// System includes

// External includes

// Project includes
#include "custom_elements/gnl_beam_lumpe_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"



namespace Kratos
{
	GnlBeamLumpe3D2N::GnlBeamLumpe3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
	}

	GnlBeamLumpe3D2N::GnlBeamLumpe3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer GnlBeamLumpe3D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new GnlBeamLumpe3D2N(
			NewId, rGeom.Create(rThisNodes), pProperties));
	}

	GnlBeamLumpe3D2N::~GnlBeamLumpe3D2N() {}

	void GnlBeamLumpe3D2N::EquationIdVector(EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo) {
		if (rResult.size() != msElementSize) rResult.resize(msElementSize);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msNumberOfNodes * msDimension;
			rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y)
				.EquationId();
			rResult[index + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Z)
				.EquationId();
			rResult[index + 2] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X)
				.EquationId();

			rResult[index + 3] = this->GetGeometry()[i].GetDof(ROTATION_Y)
				.EquationId();
			rResult[index + 4] = this->GetGeometry()[i].GetDof(ROTATION_Z)
				.EquationId();
			rResult[index + 5] = this->GetGeometry()[i].GetDof(ROTATION_X)
				.EquationId();
		}

    }
    
    void GnlBeamLumpe3D2N::GetDofList(DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo) {

		if (rElementalDofList.size() != msElementSize) {
			rElementalDofList.resize(msElementSize);
		}

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msNumberOfNodes * msDimension;
			rElementalDofList[index] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Y);
			rElementalDofList[index + 1] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Z);
			rElementalDofList[index + 2] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_X);

			rElementalDofList[index + 3] = this->GetGeometry()[i]
				.pGetDof(ROTATION_Y);
			rElementalDofList[index + 4] = this->GetGeometry()[i]
				.pGetDof(ROTATION_Z);
			rElementalDofList[index + 5] = this->GetGeometry()[i]
				.pGetDof(ROTATION_X);
		}
    }
    
    void GnlBeamLumpe3D2N::Initialize() {}

	void GnlBeamLumpe3D2N::GetValuesVector(Vector& rValues, int Step) {
        
        KRATOS_TRY
        if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

        for (int i = 0; i < msNumberOfNodes; ++i)
        {
            int index = i * msDimension * 2;
            rValues[index] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
            rValues[index + 1] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
            rValues[index + 2] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(DISPLACEMENT_X, Step);

            rValues[index + 3] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ROTATION_Y, Step);
            rValues[index + 4] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ROTATION_Z, Step);
            rValues[index + 5] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ROTATION_X, Step);
        }
        KRATOS_CATCH("")
    }
        
    void GnlBeamLumpe3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
    {

        KRATOS_TRY
        if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

        for (int i = 0; i < msNumberOfNodes; ++i)
        {
            int index = i * msDimension * 2;
            rValues[index] = this->GetGeometry()[i].
                FastGetSolutionStepValue(VELOCITY_Y, Step);
            rValues[index + 1] = this->GetGeometry()[i].
                FastGetSolutionStepValue(VELOCITY_Z, Step);
            rValues[index + 2] = this->GetGeometry()[i].
                FastGetSolutionStepValue(VELOCITY_X, Step);

            rValues[index + 3] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_VELOCITY_Y, Step);
            rValues[index + 4] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_VELOCITY_Z, Step);
            rValues[index + 5] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_VELOCITY_X, Step);
        }

        KRATOS_CATCH("")
    }

    void GnlBeamLumpe3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
    {

        KRATOS_TRY
        if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

        for (int i = 0; i < msNumberOfNodes; ++i)
        {
            int index = i * msDimension * 2;

            rValues[index] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ACCELERATION_Y, Step);
            rValues[index + 1] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ACCELERATION_Z, Step);
            rValues[index + 2] = this->GetGeometry()[i]
                .FastGetSolutionStepValue(ACCELERATION_X, Step);

            rValues[index + 3] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_ACCELERATION_Y, Step);
            rValues[index + 4] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z, Step);
            rValues[index + 5] = this->GetGeometry()[i].
                FastGetSolutionStepValue(ANGULAR_ACCELERATION_X, Step);
        }
        KRATOS_CATCH("")
    }

    void GnlBeamLumpe3D2N::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) 
        {
            this->CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
            this->CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

            bounded_matrix<double,msDimension,msDimension> RotationMatrix0 = this->RotationMatrix0();
            KRATOS_WATCH(RotationMatrix0);

            bounded_vector<double,3> psi_k = ZeroVector(3);
            psi_k[1] = -0.5235987756;
            bounded_matrix<double,msDimension,msDimension> RR = this->RotationMatrix(psi_k);
            KRATOS_WATCH(RR);            

        }

    void GnlBeamLumpe3D2N::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) 
        {
            // testing (this is linear analysis) and already new K(k+1)!
            Vector Nodal_Deformation = ZeroVector(msElementSize);
            this->GetValuesVector(Nodal_Deformation);

            Matrix LeftHandSide = ZeroMatrix(msElementSize,msElementSize);
            this->CalculateLeftHandSide(LeftHandSide,rCurrentProcessInfo);


            rRightHandSideVector = ZeroVector(msElementSize);
            rRightHandSideVector -= prod(LeftHandSide,Nodal_Deformation);

        }

    void GnlBeamLumpe3D2N::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) 
        {
            rLeftHandSideMatrix = ZeroMatrix(msElementSize,msElementSize);
            rLeftHandSideMatrix = this->StiffnessMatrix();

            bounded_matrix<double,msDimension,msDimension> RotationMatrix = this->RotationMatrix0();
            bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = ZeroMatrix(msElementSize,msElementSize);
            this->AssembleTransformationMatrix(RotationMatrix,TransformationMatrix);

            rLeftHandSideMatrix = prod(rLeftHandSideMatrix,TransformationMatrix);
            rLeftHandSideMatrix = prod(Matrix(trans(TransformationMatrix)),rLeftHandSideMatrix);
            
        }
  
    //////////////////// custom functions ////////////////////
    
    double GnlBeamLumpe3D2N::CalculateReferenceLength() 
    {  
        KRATOS_TRY;
        const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
        const double L = sqrt(dx*dx + dy*dy + dz*dz);
        return L;
        KRATOS_CATCH("")
    }

    bounded_matrix<double,GnlBeamLumpe3D2N::msLocalSize,
    GnlBeamLumpe3D2N::msElementSize> GnlBeamLumpe3D2N::Mode_TransformationMatrix()
    {
        bounded_matrix<double,msLocalSize,msElementSize>
         A = ZeroMatrix(msLocalSize, msElementSize);

        const double length_k_inv = 1.00 / this->GetGeometry().Length();
        // LINEAR: const double length_k_inv = 1.00 / this->CalculateReferenceLength

        A(0,1) = -length_k_inv;
        A(0,3) = 1.00;
        A(0,7) = length_k_inv;
        A(1,0) = length_k_inv;
        A(1,4) = 1.00;
        A(1,6) = -length_k_inv;
        A(2,2) = -1.00;
        A(2,8) = 1.00;
        A(3,1) = -length_k_inv;
        A(3,7) = length_k_inv;
        A(3,9) = 1.00;
        A(4,0) = length_k_inv;
        A(4,6) = -length_k_inv;
        A(4,10) = 1.00;
        A(5,5) = -1.00;
        A(5,11) = 1.00;
        return A;
    }

    bounded_matrix<double,GnlBeamLumpe3D2N::msLocalSize,
    GnlBeamLumpe3D2N::msLocalSize> GnlBeamLumpe3D2N::Delta_StiffnessMatrix()
    {
        bounded_matrix<double,msLocalSize,msLocalSize>
         K_delta = ZeroMatrix(msLocalSize, msLocalSize);

        const double L_ref = this->CalculateReferenceLength();


        K_delta(0,0) = E*Iy*((3/alphaz)+1)/L_ref;
        K_delta(0,3) = E*Iy*((3/alphaz)-1)/L_ref;

        K_delta(1,1) = E*Iz*((3/alphay)+1)/L_ref;
        K_delta(1,4) = E*Iz*((3/alphay)-1)/L_ref;

        K_delta(2,2) = E*A/L_ref;

        K_delta(3,0) = E*Iy*((3/alphaz)-1)/L_ref;
        K_delta(3,3) = E*Iy*((3/alphaz)+1)/L_ref;

        K_delta(4,1) = E*Iz*((3/alphay)-1)/L_ref;
        K_delta(4,4) = E*Iz*((3/alphay)+1)/L_ref;

        K_delta(5,5) = G*It/(L_ref*mu);    

        return K_delta;
    }

    bounded_matrix<double,GnlBeamLumpe3D2N::msElementSize,
    GnlBeamLumpe3D2N::msElementSize> GnlBeamLumpe3D2N::StiffnessMatrix()
    {
        bounded_matrix<double,msElementSize,msElementSize>
         K = ZeroMatrix(msElementSize, msElementSize);
        bounded_matrix<double,msLocalSize,msElementSize> A = this->Mode_TransformationMatrix();
        bounded_matrix<double,msLocalSize,msElementSize> A_transpose = Matrix(trans(A));
        bounded_matrix<double,msLocalSize,msLocalSize> K_delta = this->Delta_StiffnessMatrix();
        K = prod(K_delta,A);
        K = prod(A_transpose,K);
        return K;
    }

    bounded_matrix<double,GnlBeamLumpe3D2N::msDimension,GnlBeamLumpe3D2N::msDimension>
    GnlBeamLumpe3D2N::RotationMatrix(bounded_vector<double,msDimension> Psi_k)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        const double Norm_Psi = MathUtils<double>::Norm(Psi_k);

        identity_matrix<double> Eye_Matrix (msDimension);
        
        if (Norm_Psi < numerical_limit) return Eye_Matrix;

        else
        {
            bounded_matrix<double,msDimension,msDimension> RotationMatrix =
            ZeroMatrix(msDimension,msDimension);

            bounded_matrix<double, msDimension,msDimension> S =
             ZeroMatrix(msDimension,msDimension);

            S(0,1) = Psi_k[2];
            S(0,2) = -Psi_k[1];
            S(1,0) = -Psi_k[2];
            S(1,2) = Psi_k[0];
            S(2,0) = Psi_k[1];
            S(2,1) = -Psi_k[0];

            bounded_matrix<double, msDimension,msDimension> S2 = prod(S,S);

            RotationMatrix = Eye_Matrix;
            RotationMatrix += S * (std::sin(Norm_Psi)/Norm_Psi);            
            RotationMatrix += S2 * ((1.00 - std::cos(Norm_Psi)) / (Norm_Psi*Norm_Psi));           
            return RotationMatrix;
        }    

    }

    bounded_matrix<double,GnlBeamLumpe3D2N::msDimension,GnlBeamLumpe3D2N::msDimension>
    GnlBeamLumpe3D2N::RotationMatrix0()
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        bounded_matrix<double,msDimension,msDimension> RotationMatrix =
        ZeroMatrix(msDimension,msDimension);   

		array_1d<double, msDimension> DirectionVectorG3 = ZeroVector(msDimension);
		array_1d<double, msDimension> DirectionVectorG2 = ZeroVector(msDimension);
		array_1d<double, msDimension> DirectionVectorG1 = ZeroVector(msDimension);
		array_1d<double, msLocalSize> ReferenceCoordinates = ZeroVector(msLocalSize);

		ReferenceCoordinates[0] = this->GetGeometry()[0].Y0();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Z0();
		ReferenceCoordinates[2] = this->GetGeometry()[0].X0();
		ReferenceCoordinates[3] = this->GetGeometry()[1].Y0();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Z0();
        ReferenceCoordinates[5] = this->GetGeometry()[1].X0();


        for (unsigned int i = 0; i < msDimension; ++i)
		{
			DirectionVectorG3[i] = (ReferenceCoordinates[i + msDimension]
				- ReferenceCoordinates[i]);
		}
        
        // take user defined local axis 2 from GID input
        if (this->Has(LOCAL_AXIS_2)) 
        {
            double VectorNorm = MathUtils<double>::Norm(DirectionVectorG3);
            if (VectorNorm > numerical_limit) DirectionVectorG3 /= VectorNorm;

            array_1d<double, msLocalSize> TempLocalAxis2  = this->GetValue(LOCAL_AXIS_2);
            DirectionVectorG1[0] = TempLocalAxis2[1];
            DirectionVectorG1[1] = TempLocalAxis2[2];
            DirectionVectorG1[2] = TempLocalAxis2[0];

            VectorNorm = MathUtils<double>::Norm(DirectionVectorG1);
            if (VectorNorm > numerical_limit) DirectionVectorG1 /= VectorNorm;

            DirectionVectorG2[0] = DirectionVectorG3[1]*DirectionVectorG1[2]-DirectionVectorG3[2]*DirectionVectorG1[1];
            DirectionVectorG2[1] = DirectionVectorG3[2]*DirectionVectorG1[0]-DirectionVectorG3[0]*DirectionVectorG1[2];
            DirectionVectorG2[2] = DirectionVectorG3[0]*DirectionVectorG1[1]-DirectionVectorG3[1]*DirectionVectorG1[0];

            VectorNorm = MathUtils<double>::Norm(DirectionVectorG2);
            if (VectorNorm > numerical_limit) DirectionVectorG2 /= VectorNorm;
            else KRATOS_ERROR << "LOCAL_AXIS_3 has length 0 for element " << this->Id() << std::endl;

            for (int i = 0; i < msDimension; ++i)
            {
                RotationMatrix(0, i) = DirectionVectorG1[i];
                RotationMatrix(1, i) = DirectionVectorG2[i];
                RotationMatrix(2, i) = DirectionVectorG3[i];
            }
        }
        else KRATOS_ERROR << "Local Coordinate System not defined !!" << this->Id() << std::endl;
        return RotationMatrix;
    }

    void GnlBeamLumpe3D2N::AssembleTransformationMatrix(Matrix RotationMatrix, bounded_matrix<double,
        GnlBeamLumpe3D2N::msElementSize,GnlBeamLumpe3D2N::msElementSize>& TransformationMatrix)
        {
        for (unsigned int kk = 0; kk < msElementSize; kk += msDimension)
        {
            for (int i = 0; i<msDimension; ++i)
            {
                for (int j = 0; j<msDimension; ++j)
                {
                    TransformationMatrix(i + kk, j + kk) = RotationMatrix(i, j);
                }
            }
        }
        }

} // namespace Kratos.



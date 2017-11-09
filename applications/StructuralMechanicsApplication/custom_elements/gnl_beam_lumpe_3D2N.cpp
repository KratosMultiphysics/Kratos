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

//TODO::    - find G1 for R_0(psi) --> RotationMatrix0()
//          - do eigenwert analysis with eigenwert +1 --> UpdateIncrementRotation_i()

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


            // testing !!
            bounded_vector<double,3> psi_1 = ZeroVector(3);
            bounded_vector<double,3> psi_2 = ZeroVector(3);
            psi_1[0]=1;
            psi_1[1]=2;
            psi_1[2]=3;
            psi_2[0]=5;
            psi_2[1]=7;
            psi_2[2]=8;
 
            double innerprod = inner_prod(psi_1,psi_2);
            KRATOS_WATCH(innerprod);

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

            bounded_matrix<double,msDimension,msDimension> RotationMatrix = this->RotationMatrixR();
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
    GnlBeamLumpe3D2N::RotationMatrix0(
        bounded_vector<double,msDimension> G_3_0,
        bounded_vector<double,msDimension> G_1_0)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        /////// Calculate euler angles (see Lumpe beam p.19 ff.)
        // Get current nodal basevector G3
        bounded_vector<double,msDimension> G_3_i_k = this->CurrentBaseVector3();
        bounded_vector<double,msDimension> G_1_head = MathUtils<double>::CrossProduct(G_3_0,G_3_i_k);

        double VectorNorm = MathUtils<double>::Norm(G_1_head);
        if (VectorNorm > numerical_limit) G_1_head /= VectorNorm;

        //get euler angles
        const double cos_psi_P = inner_prod(G_1_head,G_1_0);
        const double sin_psi_P = sqrt(1.00 - (cos_psi_P*cos_psi_P));

        const double cos_psi_N = inner_prod(G_3_i_k,G_3_0);
        const double sin_psi_N = sqrt(1.00 - (cos_psi_N*cos_psi_N));     
        
        
        ///// get G1 or G2 ???? HOW ???
        const double cos_psi_E = 0.00; 
        const double sin_psi_E = sqrt(1.00 - (cos_psi_E*cos_psi_E));  


        //calculate R0
        bounded_matrix<double,msDimension,msDimension> R_psi_p = ZeroMatrix(msDimension,msDimension);
        bounded_matrix<double,msDimension,msDimension> R_psi_n = ZeroMatrix(msDimension,msDimension);
        bounded_matrix<double,msDimension,msDimension> R_psi_e = ZeroMatrix(msDimension,msDimension);
        bounded_matrix<double,msDimension,msDimension> R_0 = ZeroMatrix(msDimension,msDimension);

        R_psi_p(0,0) = cos_psi_P;
        R_psi_p(0,1) = sin_psi_P;
        R_psi_p(1,0) = -sin_psi_P;
        R_psi_p(1,1) = cos_psi_P;
        R_psi_p(2,2) = 1.00;

        R_psi_n(0,0) = 1.00;
        R_psi_n(1,1) = cos_psi_N;
        R_psi_n(1,2) = sin_psi_N;
        R_psi_n(2,1) = -sin_psi_N;
        R_psi_n(2,2) = cos_psi_N;

        R_psi_e(0,0) = cos_psi_E;
        R_psi_e(0,1) = sin_psi_E;
        R_psi_e(1,0) = -sin_psi_E;
        R_psi_e(1,1) = cos_psi_E;
        R_psi_e(2,2) = 1.00;     

        R_0 = prod(R_psi_n,R_psi_p);
        R_0 = prod(R_psi_n,R_psi_e);

        return R_0;

    }
    //RotationMatrixR :: G_0 = R_R * G_R
    bounded_matrix<double,GnlBeamLumpe3D2N::msDimension,GnlBeamLumpe3D2N::msDimension>
    GnlBeamLumpe3D2N::RotationMatrixR()
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
    GnlBeamLumpe3D2N::msElementSize,GnlBeamLumpe3D2N::msElementSize>&TransformationMatrix)
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

    bounded_vector<double,GnlBeamLumpe3D2N::msLocalSize> 
    GnlBeamLumpe3D2N::UpdateIncrementDeformation()
    {
        // Get v_0
        Vector Nodal_Deformation = ZeroVector(msElementSize);
        this->GetValuesVector(Nodal_Deformation);

        bounded_matrix<double,msDimension,msDimension> R_R = this->RotationMatrixR();
        bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix_R = ZeroMatrix(msElementSize,msElementSize);
        this->AssembleTransformationMatrix(R_R,TransformationMatrix_R);    
        Vector Nodal_Deformation_0 = prod(TransformationMatrix_R,Nodal_Deformation); /// ????


        // Get G_0
        bounded_vector<double,msDimension> G_1_0 = ZeroVector(msDimension);
        bounded_vector<double,msDimension> G_2_0 = ZeroVector(msDimension);
        bounded_vector<double,msDimension> G_3_0 = ZeroVector(msDimension);

        for (int i = 0;i < msDimension;++i)
        {
            G_1_0[i] = R_R(0, i);
            G_2_0[i] = R_R(1, i);
            G_3_0[i] = R_R(2, i);
        }


        // UpdateIncrementDef
        bounded_vector<double,msLocalSize> IncrementDeformation = ZeroVector(msLocalSize);
        // dv = [d_phi_1_a,d_phi_2_a,d_u_3_b,d_phi_1_b,d_phi_2_b,d_phi_3_b]
        bounded_vector<double,msDimension> Phi_i = ZeroVector(msDimension);
        bounded_vector<double,msDimension> d_Phi_i = ZeroVector(msDimension);
        IncrementDeformation[2] = this->GetGeometry().Length() - this->CalculateReferenceLength();

        for (int i = 0; i < msNumberOfNodes; ++i)
        {
            int index = i*msLocalSize;
            
            for (int j = 0;j< msDimension;++j)
            {
                Phi_i[j] = Nodal_Deformation_0[index+msDimension+j];
            }
            d_Phi_i = this->UpdateIncrementRotation_i(Phi_i,G_1_0,G_2_0,G_3_0);

            if (i == 0) for (int k=0;k<2;++k) IncrementDeformation[k] = d_Phi_i[k];
            if (i == 1) for (int k=0;k<3;++k) IncrementDeformation[3+k] = d_Phi_i[k];
        }

        return IncrementDeformation;
    }


    bounded_vector<double,GnlBeamLumpe3D2N::msDimension>
    GnlBeamLumpe3D2N::UpdateIncrementRotation_i(bounded_vector<double,
        GnlBeamLumpe3D2N::msDimension> Psi_i,
        bounded_vector<double,msDimension> G_1_0,
        bounded_vector<double,msDimension> G_2_0,
        bounded_vector<double,msDimension> G_3_0)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        bounded_matrix<double,msDimension,msDimension> R_i = this->RotationMatrix(Psi_i);
        bounded_matrix<double,msDimension,msDimension> R_0 = this->RotationMatrix0(G_3_0,G_1_0);

        bounded_matrix<double,msDimension,msDimension> R_d = prod(R_i,Matrix(trans(R_0)));
        
        
        ///EIGENWERT ANALYIS OF R_Dwith eigenwert = 1.00
        bounded_vector<double,msDimension> Phi_d = ZeroVector(msDimension);
        const double VectorNorm = MathUtils<double>::Norm(Phi_d);
        if (VectorNorm > numerical_limit) Phi_d /= VectorNorm;
        
        //NORM OF PHI_D
        double norm_phi = 0.00;
        for (int i=0;i<msDimension;++i) norm_phi += R_i(i,i);
        norm_phi = (0.50 * norm_phi) - 1.00;
        norm_phi = std::acos(norm_phi);

        Phi_d *= norm_phi;

        return Phi_d;
    }

    bounded_vector<double,GnlBeamLumpe3D2N::msDimension> GnlBeamLumpe3D2N::CurrentBaseVector3()
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
		array_1d<double, msDimension> DirectionVectorG3 = ZeroVector(msDimension);
		array_1d<double, msLocalSize> CurrentCoordinates = ZeroVector(msLocalSize);

		CurrentCoordinates[0] = this->GetGeometry()[0].Y();
		CurrentCoordinates[1] = this->GetGeometry()[0].Z();
		CurrentCoordinates[2] = this->GetGeometry()[0].X();
		CurrentCoordinates[3] = this->GetGeometry()[1].Y();
		CurrentCoordinates[4] = this->GetGeometry()[1].Z();
        CurrentCoordinates[5] = this->GetGeometry()[1].X();


        for (unsigned int i = 0; i < msDimension; ++i)
		{
			DirectionVectorG3[i] = (CurrentCoordinates[i + msDimension]
				- CurrentCoordinates[i]);
        }
        
        const double VectorNorm = MathUtils<double>::Norm(DirectionVectorG3);
        if (VectorNorm > numerical_limit) DirectionVectorG3 /= VectorNorm;

        return DirectionVectorG3;
    }
} // namespace Kratos.



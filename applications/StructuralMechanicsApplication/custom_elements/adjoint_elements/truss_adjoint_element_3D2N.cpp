// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Long Chen
//                   
//                   
//
#include "custom_elements/adjoint_elements/truss_adjoint_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{
	TrussAdjointElement3D2N::TrussAdjointElement3D2N(IndexType NewId, 
									GeometryType::Pointer pGeometry,
									bool rLinear)
									: TrussElement3D2N(NewId, pGeometry)
	{
		this->mIsLinearElement = rLinear;
	}

	TrussAdjointElement3D2N::TrussAdjointElement3D2N(IndexType NewId,
									GeometryType::Pointer pGeometry,
									PropertiesType::Pointer pProperties,
									bool rLinear) 
									: TrussElement3D2N(NewId, pGeometry, pProperties)
	{
		this->mIsLinearElement = rLinear;
	}

	Element::Pointer TrussAdjointElement3D2N::Create(IndexType NewId,
									NodesArrayType const& rThisNodes,
									PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new TrussAdjointElement3D2N(
			NewId, rGeom.Create(rThisNodes), pProperties, this->mIsLinearElement));
	}

    Element::Pointer TrussAdjointElement3D2N::Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(
            new TrussAdjointElement3D2N(NewId, pGeom, pProperties, this->mIsLinearElement));
        KRATOS_CATCH("")
    }



	TrussAdjointElement3D2N::~TrussAdjointElement3D2N(){}

	void TrussAdjointElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
									ProcessInfo& rCurrentProcessInfo){

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		if (rResult.size() != local_size) rResult.resize(local_size);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * 3;
			rResult[index] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X)
				.EquationId();
			rResult[index + 1] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y)
				.EquationId();
			rResult[index + 2] = this->GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z)
				.EquationId();
		}

	}
	void TrussAdjointElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
									ProcessInfo& rCurrentProcessInfo){

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		if (rElementalDofList.size() != local_size) {
			rElementalDofList.resize(local_size);
		}

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * 3;
			rElementalDofList[index] = this->GetGeometry()[i]
				.pGetDof(ADJOINT_DISPLACEMENT_X);
			rElementalDofList[index + 1] = this->GetGeometry()[i]
				.pGetDof(ADJOINT_DISPLACEMENT_Y);
			rElementalDofList[index + 2] = this->GetGeometry()[i]
				.pGetDof(ADJOINT_DISPLACEMENT_Z);
		}
	}

    double TrussAdjointElement3D2N::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
    {
        KRATOS_TRY;

        if (this->GetProperties().Has(rDesignVariable))
        {
            const double variable_value = this->GetProperties()[rDesignVariable];
            return variable_value;
        }
        else
            return 1.0;

        KRATOS_CATCH("");
    }

    double TrussAdjointElement3D2N::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double, 3>>& rDesignVariable)
    {
        KRATOS_TRY;

        if (rDesignVariable == SHAPE_SENSITIVITY)
        {
            double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
            double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
            double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
            double L = sqrt(dx*dx + dy*dy + dz*dz);
            return L;
        }
        else
            return 1.0;

        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // define working variables
        Vector RHS_undist;
        Vector RHS_dist;
        ProcessInfo testProcessInfo = rCurrentProcessInfo;

        // Compute RHS before disturbing
        this->CalculateRightHandSide(RHS_undist, testProcessInfo);

        rOutput.resize(1, RHS_undist.size());

        // Get disturbance measure
        double delta = this->GetValue(DISTURBANCE_MEASURE);
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        if (this->GetProperties().Has(rDesignVariable))
        {
            // Save properties and its pointer
            Properties& r_global_property = this->GetProperties();
            Properties::Pointer p_global_properties = this->pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(new Properties(r_global_property));
            this->SetProperties(p_local_property);

            // Disturb the design variable
            const double current_property_value = this->GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

            // compute RHS after disturbance
            this->CalculateRightHandSide(RHS_dist, testProcessInfo);

            rOutput.resize(1, RHS_dist.size());

            // Compute derivative of RHS w.r.t. design variable with finite differences
            RHS_dist -= RHS_undist;
            RHS_dist /= delta;
            for (unsigned int i = 0; i < RHS_dist.size(); ++i)
                rOutput(0, i) = RHS_dist[i];

            // Give element original properties back
            this->SetProperties(p_global_properties);

            // Compute RHS again in order to ensure that changed member variables like mLHS get back to their origin values
            this->CalculateRightHandSide(RHS_dist, testProcessInfo);
        }
        else
            rOutput.clear();


        KRATOS_CATCH("");

    }

    void TrussAdjointElement3D2N::CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable, Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // define working variables
        Vector RHS_undist;
        Vector RHS_dist;
        ProcessInfo testProcessInfo = rCurrentProcessInfo;

        // Get disturbance measure
        double delta = this->GetValue(DISTURBANCE_MEASURE);
        
        double correction_factor = this->GetDisturbanceMeasureCorrectionFactor(rDesignVariable);
        delta *= correction_factor;

        std::cout << "#######################################" << std::endl;
        KRATOS_WATCH(this->GetId());

        if (rDesignVariable == SHAPE_SENSITIVITY)
        {
            const int number_of_nodes = GetGeometry().PointsNumber();
            const int dimension = this->GetGeometry().WorkingSpaceDimension();
            const int local_size = number_of_nodes*dimension * 1;

            rOutput.resize(dimension*number_of_nodes, local_size);

            // compute RHS before disturbing
            //this->mIsLinearElement = true;

            std::cout << "calculate RightHandside before disturbance!" << std::endl;

            this->CalculateRightHandSide(RHS_undist, testProcessInfo);

            //this->CalculateRightHandSideUnDist(RHS_undist, testProcessInfo);

            KRATOS_WATCH(RHS_undist);


            //this->mIsLinearElement = false;
            for (int j = 0; j < number_of_nodes; j++)
            {
                KRATOS_WATCH(j);
                

                //begin: derive w.r.t. x-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].X0() += delta;


                // compute RHS after disturbance
                
                this->CalculateRightHandSide(RHS_dist, testProcessInfo);
                //this->CalculateRightHandSideUnDist(RHS_dist, testProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                RHS_dist -= RHS_undist;
                RHS_dist /= delta;
                for (unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput((0 + j*dimension), i) = RHS_dist[i];

                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].X0() -= delta;

                //end: derive w.r.t. x-coordinate-----------------------------------------------------
                
                //begin: derive w.r.t. y-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].Y0() += delta;
                std::cout << "#######" << std::endl;
                std::cout << "Y: compute RHS after disturbance" << std::endl;

                double y0 = this->GetGeometry()[j].Y0();

                KRATOS_WATCH(y0);

                // compute RHS after disturbance
                this->CalculateRightHandSide(RHS_dist, testProcessInfo);
                //this->CalculateRightHandSideUnDist(RHS_dist, testProcessInfo);
                KRATOS_WATCH(RHS_dist);

                std::cout << "Y: compute dRHS" << std::endl;

                //compute derivative of RHS w.r.t. design variable with finite differences
                RHS_dist -= RHS_undist;
                RHS_dist /= delta;
                KRATOS_WATCH(RHS_dist);
                for (unsigned int i = 0; i < RHS_dist.size(); i++)
                    rOutput((1 + j*dimension), i) = RHS_dist[i];

                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].Y0() -= delta;

                //end: derive w.r.t. y-coordinate-----------------------------------------------------
               
                //begin: derive w.r.t. z-coordinate---------------------------------------------------

                // disturb the design variable
                this->GetGeometry()[j].Z0() += delta;
                // compute RHS after disturbance
                this->CalculateRightHandSide(RHS_dist, testProcessInfo);
                //this->CalculateRightHandSideUnDist(RHS_dist, testProcessInfo);
                //compute derivative of RHS w.r.t. design variable with finite differences
                RHS_dist -= RHS_undist;
                RHS_dist /= delta;
                for (unsigned int i = 0; i < RHS_dist.size(); i++)
                {
                    rOutput((2 + j*dimension), i) = RHS_dist[i];  //LChen tmp
                    //rOutput((2 + j*dimension), i) = 0.0;
                }                   

                //KRATOS_WATCH(rOutput);


                // Reset pertubed vector
                RHS_dist = Vector(0);

                // undisturb the design variable
                this->GetGeometry()[j].Z0() -= delta;

                // Compute RHS again in order to ensure that changed member variables like mLHS get back their origin values
                this->CalculateRightHandSide(RHS_dist, testProcessInfo);

                //end: derive w.r.t. z-coordinate-----------------------------------------------------


               

            }// end loop over element nodes

        }
        else
            KRATOS_ERROR << "Unsupported design variable!" << std::endl;

        KRATOS_CATCH("");
    }

    /*
    void TrussAdjointElement3D2N::CalculateRightHandSideUnDist(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        const int NumNodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const int LocalSize = NumNodes * dimension;

        rRightHandSideVector = ZeroVector(LocalSize);

        VectorType InternalForces = ZeroVector(LocalSize);
        //this->UpdateInternalForces(InternalForces);

        //###############################

        //const int NumNodes = this->GetGeometry().PointsNumber();
        //const int dimension = this->GetGeometry().WorkingSpaceDimension();
        //const int LocalSize = NumNodes * dimension;

        MatrixType TransformationMatrix = ZeroMatrix(LocalSize, LocalSize);
        this->CreateTransformationMatrix(TransformationMatrix);
        //const double InternalStrainGL = this->CalculateGreenLagrangeStrain();


        const double l = this->CalculateReferenceLength();  // only for the undisturbed calculation
        const double L0 = this->CalculateReferenceLength();
        const double E = this->GetProperties()[YOUNG_MODULUS];
        const double A = this->GetProperties()[CROSS_AREA];

        const double InternalStrainGL = ((l * l - L0 * L0) / (2.00 * L0 * L0));
        //const double InternalStrainGL = this->CalculateGreenLagrangeStrain();
        KRATOS_WATCH(InternalStrainGL);

        double S_pre = 0.00;
        if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
            S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        const double N = ((E*InternalStrainGL + S_pre) * l * A) / L0;

        double Spk2 = E*InternalStrainGL;
        KRATOS_WATCH(Spk2);
        KRATOS_WATCH(l);
        KRATOS_WATCH(L0);
        KRATOS_WATCH(A);
        KRATOS_WATCH(N);

        if (N < 0.00) this->mIsCompressed = true;
        else this->mIsCompressed = false;

        //internal force vectors
        VectorType f_local = ZeroVector(LocalSize);
        f_local[0] = -1.00 * N;
        f_local[3] = 1.00 * N;

        InternalForces = ZeroVector(LocalSize);

        std::cout << "before transformation" << std::endl;

        KRATOS_WATCH(TransformationMatrix);
        KRATOS_WATCH(f_local);

        noalias(InternalForces) = prod(TransformationMatrix, f_local);
        std::cout << "after transformation" << std::endl;



        //################################





        rRightHandSideVector -= InternalForces;
        //KRATOS_WATCH(rRightHandSideVector);
        if (this->mIsLinearElement == true)
        {
            Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
            this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
            Vector NodalDeformation = ZeroVector(LocalSize);
            this->GetValuesVector(NodalDeformation);
            rRightHandSideVector = ZeroVector(LocalSize);
            rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
        }

        //add bodyforces 
        rRightHandSideVector += this->CalculateBodyForces();

        //VectorType BodyForces = this->CalculateBodyForces();
        //KRATOS_WATCH(InternalForces);
        //KRATOS_WATCH(BodyForces);



        if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
            rRightHandSideVector = ZeroVector(LocalSize);
        }

        KRATOS_CATCH("")

    }

    void TrussAdjointElement3D2N::CalculateRightHandSideDist(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        const int NumNodes = this->GetGeometry().PointsNumber();
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        const int LocalSize = NumNodes * dimension;

        rRightHandSideVector = ZeroVector(LocalSize);

        VectorType InternalForces = ZeroVector(LocalSize);
        //this->UpdateInternalForces(InternalForces);

        //###############################

        //const int NumNodes = this->GetGeometry().PointsNumber();
        //const int dimension = this->GetGeometry().WorkingSpaceDimension();
        //const int LocalSize = NumNodes * dimension;

        MatrixType TransformationMatrix = ZeroMatrix(LocalSize, LocalSize);
        this->CreateTransformationMatrix(TransformationMatrix);
        //const double InternalStrainGL = this->CalculateGreenLagrangeStrain();


        const double l = this->CalculateReferenceLength();  // only for the undisturbed calculation
        const double L0 = this->CalculateReferenceLength();
        const double E = this->GetProperties()[YOUNG_MODULUS];
        const double A = this->GetProperties()[CROSS_AREA];

        const double InternalStrainGL = ((l * l - L0 * L0) / (2.00 * L0 * L0));
        KRATOS_WATCH(InternalStrainGL);

        double S_pre = 0.00;
        if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
            S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        const double N = ((E*InternalStrainGL + S_pre) * l * A) / L0;

        double Spk2 = E*InternalStrainGL;
        KRATOS_WATCH(Spk2);
        KRATOS_WATCH(l);
        KRATOS_WATCH(L0);
        KRATOS_WATCH(A);
        KRATOS_WATCH(N);

        if (N < 0.00) this->mIsCompressed = true;
        else this->mIsCompressed = false;

        //internal force vectors
        VectorType f_local = ZeroVector(LocalSize);
        f_local[0] = -1.00 * N;
        f_local[3] = 1.00 * N;

        InternalForces = ZeroVector(LocalSize);

        std::cout << "before transformation" << std::endl;

        KRATOS_WATCH(TransformationMatrix);
        KRATOS_WATCH(f_local);

        noalias(InternalForces) = prod(TransformationMatrix, f_local);
        std::cout << "after transformation" << std::endl;



        //################################





        rRightHandSideVector -= InternalForces;
        //KRATOS_WATCH(rRightHandSideVector);
        if (this->mIsLinearElement == true)
        {
            Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
            this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
            Vector NodalDeformation = ZeroVector(LocalSize);
            this->GetValuesVector(NodalDeformation);
            rRightHandSideVector = ZeroVector(LocalSize);
            rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
        }

        //add bodyforces 
        rRightHandSideVector += this->CalculateBodyForces();

        //VectorType BodyForces = this->CalculateBodyForces();
        //KRATOS_WATCH(InternalForces);
        //KRATOS_WATCH(BodyForces);



        if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
            rRightHandSideVector = ZeroVector(LocalSize);
        }

        KRATOS_CATCH("")
    }
    */



    void TrussAdjointElement3D2N::Calculate(const Variable<Vector >& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::Calculate(const Variable<Matrix >& rVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable, Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
        const Variable<Vector>& rStressVariable, Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::CalculateStressDesignVariableDerivative(const Variable<array_1d<double, 3>>& rDesignVariable,
        const Variable<Vector>& rStressVariable,
        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }

    void TrussAdjointElement3D2N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        std::cout << "this is not implemented yet!" << std::endl;
        KRATOS_CATCH("");
    }


    /*
	void TrussAdjointElement3D2N::Initialize() {

		KRATOS_TRY
		KRATOS_CATCH("")
	}

	TrussAdjointElement3D2N::MatrixType TrussAdjointElement3D2N::CreateElementStiffnessMatrix(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		const double E = this->GetProperties()[YOUNG_MODULUS];
		double A = this->GetProperties()[CROSS_AREA];

		double S_pre = 0.00;
		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
			S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
		}

		MatrixType LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);

		// du... delta displacement in x-direction
		// dv... delta displacement in y-direction
		// dw... delta displacement in z-direction
		// L... inital member length
		// l... deformed member length
		// e_gl... green_lagrange strain

	    double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
	    double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
	    double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = this->CalculateReferenceLength();
		const double l = this->CalculateCurrentLength();
	    double e_gL = (l*l - L*L) / (2.00 * L*L);
		const double L2 = L*L;
		const double L4 = L2*L2;


		if (this->mIsLinearElement == true)
		{
			du = 0.00;
			dv = 0.00;
			dw = 0.00;
			e_gL = 0.00;
		}


		const double K_1 = e_gL*E + S_pre;

		//if cable + compressed -> no contribution to global K
		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) A = 0;

		LocalStiffnessMatrix(0, 0) = A*L*(K_1 / L2 + E*(dx + du)*(dx + du) / L4); 
		LocalStiffnessMatrix(3, 3) = LocalStiffnessMatrix(0, 0); 

		LocalStiffnessMatrix(1, 1) = A*L*(K_1 / L2 + E*(dy + dv)*(dy + dv) / L4); 
		LocalStiffnessMatrix(4, 4) = LocalStiffnessMatrix(1, 1); 

		LocalStiffnessMatrix(2, 2) = A*L*(K_1 / L2 + E*(dz + dw)*(dz + dw) / L4); 
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(2, 2); 

		LocalStiffnessMatrix(0, 1) = A*L*((dx + du)*(dy + dv)*E / L4);		
		LocalStiffnessMatrix(1, 0) = LocalStiffnessMatrix(0, 1);			 

		LocalStiffnessMatrix(0, 2) = A*L*((dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(2, 0) = LocalStiffnessMatrix(0, 2); 

		LocalStiffnessMatrix(0, 3) = A*L*(-K_1 / L2 - E*(dx + du)*(dx + du) / L4); 
		LocalStiffnessMatrix(3, 0) = LocalStiffnessMatrix(0, 3); 

		LocalStiffnessMatrix(0, 4) = A*L*((-1.00)*(dx + du)*(dy + dv)*E / L4); 
		LocalStiffnessMatrix(4, 0) = LocalStiffnessMatrix(0, 4); 

		LocalStiffnessMatrix(0, 5) = A*L*((-1.00)*(dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 0) = LocalStiffnessMatrix(0, 5); 

		LocalStiffnessMatrix(1, 2) = A*L*((dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(2, 1) = LocalStiffnessMatrix(1, 2); 

		LocalStiffnessMatrix(1, 3) = A*L*((-1.00)*(dy + dv)*(dx + du)*E / L4); 
		LocalStiffnessMatrix(3, 1) = LocalStiffnessMatrix(1, 3); 

		LocalStiffnessMatrix(1, 4) = A*L*(-K_1 / L2 - E*(dy + dv)*(dy + dv) / L4);  
		LocalStiffnessMatrix(4, 1) = LocalStiffnessMatrix(1, 4); 

		LocalStiffnessMatrix(1, 5) = A*L*((-1.00)*(dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 1) = LocalStiffnessMatrix(1, 5); 

		LocalStiffnessMatrix(2, 3) = A*L*((-1.00)*(dw + dz)*(dx + du)*E / L4); 
		LocalStiffnessMatrix(3, 2) = LocalStiffnessMatrix(2, 3); 

		LocalStiffnessMatrix(2, 4) = A*L*((-1.00)*(dw + dz)*(dy + dv)*E / L4);  
		LocalStiffnessMatrix(4, 2) = LocalStiffnessMatrix(2, 4); 

		LocalStiffnessMatrix(2, 5) = A*L*(-K_1 / L2 - E*(dz + dw)*(dz + dw) / L4); 
		LocalStiffnessMatrix(5, 2) = LocalStiffnessMatrix(2, 5); 

		LocalStiffnessMatrix(3, 4) = A*L*((dx + du)*(dy + dv)*E / L4); 
		LocalStiffnessMatrix(4, 3) = LocalStiffnessMatrix(3, 4); 

		LocalStiffnessMatrix(3, 5) = A*L*((dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 3) = LocalStiffnessMatrix(3, 5); 

		LocalStiffnessMatrix(4, 5) = A*L*((dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 4) = LocalStiffnessMatrix(4, 5); 

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

		if (rDampingMatrix.size1() != MatSize)
		{
			rDampingMatrix.resize(MatSize, MatSize, false);
		}

		noalias(rDampingMatrix) = ZeroMatrix(MatSize, MatSize);

		MatrixType StiffnessMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		MatrixType MassMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

		double alpha = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_ALPHA))
		{
			alpha = this->GetProperties()[RAYLEIGH_ALPHA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
		{
			alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
		}

		double beta = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_BETA))
		{
			beta = this->GetProperties()[RAYLEIGH_BETA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
		{
			beta = rCurrentProcessInfo[RAYLEIGH_BETA];
		}

		rDampingMatrix += alpha * MassMatrix;
		rDampingMatrix += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}

		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		const double TotalMass = A * L * rho;

		Vector LumpFact = ZeroVector(number_of_nodes);

		LumpFact = this->GetGeometry().LumpingFactors(LumpFact);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			double temp = LumpFact[i] * TotalMass;

			for (int j = 0; j < dimension; ++j)
			{
				int index = i *dimension + j;

				rMassMatrix(index, index) = temp;
			}
		}
		KRATOS_CATCH("")
	}

	TrussAdjointElement3D2N::VectorType TrussAdjointElement3D2N::CalculateBodyForces(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

		//getting shapefunctionvalues 
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		//creating necessary values 
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		double TotalMass = A * L * rho;
		VectorType BodyForcesNode = ZeroVector(dimension);
		VectorType BodyForcesGlobal = ZeroVector(MatSize);

		//assemble global Vector
		for (int i = 0; i < number_of_nodes; ++i) {
			BodyForcesNode = TotalMass*this->GetGeometry()[i]
				.FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0,i);

			for (int j = 0; j < dimension; ++j) {
				BodyForcesGlobal[(i*dimension) + j] = BodyForcesNode[j];
			}
		}
		
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	} */

	void TrussAdjointElement3D2N::GetValuesVector(Vector& rValues, int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension;
			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Z, Step);
		}
		KRATOS_CATCH("")
	}

    /*
	void TrussAdjointElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension;
			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(VELOCITY_Z, Step);
		}
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::GetSecondDerivativesVector(Vector& rValues,int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension;
			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Z, Step);
		}

		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
										VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		//calculate internal forces
		VectorType InternalForces = ZeroVector(LocalSize);
		this->UpdateInternalForces(InternalForces);
		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix();

		//create+compute RHS
		rRightHandSideVector = ZeroVector(LocalSize);
		//update Residual
		rRightHandSideVector -= InternalForces;
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();


		if (this->mIsLinearElement == true)
		{
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation, 0);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
			rRightHandSideVector += this->CalculateBodyForces();
		}

		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		rRightHandSideVector = ZeroVector(LocalSize);

		VectorType InternalForces = ZeroVector(LocalSize);
		this->UpdateInternalForces(InternalForces);
		rRightHandSideVector -= InternalForces;


		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
		}

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();

		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}

		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix();
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateOnIntegrationPoints(
										const Variable<double>& rVariable,
										std::vector<double>& rOutput,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		const GeometryType::IntegrationPointsArrayType& integration_points =
			GetGeometry().IntegrationPoints();

		if (rOutput.size() != integration_points.size()) {
			rOutput.resize(integration_points.size());
		}
		if (rVariable == TRUSS_PRESTRESS_PK2) {
			rOutput[0] = 0.00;
			if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
				rOutput[0] = this->GetProperties()[TRUSS_PRESTRESS_PK2];
			}
		}
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::CalculateOnIntegrationPoints(
										const Variable<Vector>& rVariable,
										std::vector<Vector>& rOutput,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const GeometryType::IntegrationPointsArrayType& integration_points =
			GetGeometry().IntegrationPoints();
		if (rOutput.size() != integration_points.size()) {
			rOutput.resize(integration_points.size());
		}
		if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR)
		{
			Vector Strain = ZeroVector(dimension);
			Strain[0] = this->CalculateGreenLagrangeStrain();
			Strain[1] = 0.00;
			Strain[2] = 0.00;
			rOutput[0] = Strain;
		}
		KRATOS_CATCH("")
	}

	void TrussAdjointElement3D2N::GetValueOnIntegrationPoints(
										const Variable<double>& rVariable,
										std::vector<double>& rValues,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues,
										rCurrentProcessInfo);
		KRATOS_CATCH("")
	}
	void TrussAdjointElement3D2N::GetValueOnIntegrationPoints(
										const Variable<Vector>& rVariable,
										std::vector<Vector>& rValues,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues,
										rCurrentProcessInfo);
		KRATOS_CATCH("")
	}


	bool TrussAdjointElement3D2N::ReturnIfIsCable()
	{
		KRATOS_TRY;
		bool IsCable = false;
		if (this->GetProperties().Has(TRUSS_IS_CABLE) == true) {
			IsCable = this->GetProperties()[TRUSS_IS_CABLE];
		}
		return IsCable;
		KRATOS_CATCH("")
	}

	
	double TrussAdjointElement3D2N::CalculateGreenLagrangeStrain(){

		KRATOS_TRY
		const double l = this->CalculateCurrentLength();
		const double L = this->CalculateReferenceLength();
		const double e = ((l * l - L * L) / (2.00 * L * L));
		return e;
		KRATOS_CATCH("")
	}

	double TrussAdjointElement3D2N::CalculateReferenceLength(){

		KRATOS_TRY
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}
	double TrussAdjointElement3D2N::CalculateCurrentLength(){

		KRATOS_TRY
		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) 
			+ (dw + dz)*(dw + dz));
		return l;
		KRATOS_CATCH("")
	}
	void TrussAdjointElement3D2N::UpdateInternalForces(VectorType& rinternalForces){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		MatrixType TransformationMatrix = ZeroMatrix(LocalSize, LocalSize);
		this->CreateTransformationMatrix(TransformationMatrix);
		const double InternalStrainGL = this->CalculateGreenLagrangeStrain();
		const double l = this->CalculateCurrentLength();
		const double L0 = this->CalculateReferenceLength();
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double A = this->GetProperties()[CROSS_AREA];

		double S_pre = 0.00;
		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
			S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
		}

		const double N = ((E*InternalStrainGL + S_pre) * l * A) / L0;

		if (N < 0.00) this->mIsCompressed = true;
		else this->mIsCompressed = false;

		//internal force vectors
		VectorType f_local = ZeroVector(LocalSize);
		f_local[0] = -1.00 * N;
		f_local[3] = 1.00 * N;
		
		rinternalForces = ZeroVector(LocalSize);
		noalias(rinternalForces) = prod(TransformationMatrix, f_local);
		KRATOS_CATCH("");
	}

	void TrussAdjointElement3D2N::CreateTransformationMatrix(Matrix& rRotationMatrix){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		//1st calculate transformation matrix
		Vector DirectionVectorX = ZeroVector(dimension);
		Vector DirectionVectorY = ZeroVector(dimension);
		Vector DirectionVectorZ = ZeroVector(dimension);
		Vector ReferenceCoordinates = ZeroVector(local_size);
		Vector GlobalZ = ZeroVector(dimension);
		GlobalZ[2] = 1.0;

		ReferenceCoordinates[0] = this->GetGeometry()[0].X();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z();

		for (int i = 0; i < dimension; ++i)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + dimension] -
				ReferenceCoordinates[i]);
		}
		// local x-axis (e1_local) is the beam axis  (in GID is e3_local)
		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
		if (VectorNorm != 0) DirectionVectorX /= VectorNorm;

		if (DirectionVectorX[2] == 1.00) {
			DirectionVectorY[1] = 1.0;
			DirectionVectorZ[0] = -1.0;
		}

		if (DirectionVectorX[2] == -1.00) {
			DirectionVectorY[1] = 1.0;
			DirectionVectorZ[0] = 1.0;
		}

		if (fabs(DirectionVectorX[2]) != 1.00) {

			DirectionVectorY = MathUtils<double>::CrossProduct(GlobalZ,
				DirectionVectorX);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorY);
			if (VectorNorm != 0) DirectionVectorY /= VectorNorm;

			DirectionVectorZ = MathUtils<double>::CrossProduct(DirectionVectorX,
				DirectionVectorY);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
			if (VectorNorm != 0) DirectionVectorZ /= VectorNorm;
		}

		//2nd fill big rotation matrix
		MatrixType CurrentCS = ZeroMatrix(dimension, dimension);
		for (int i = 0; i < dimension; ++i) {
			CurrentCS(i, 0) = DirectionVectorX[i];
			CurrentCS(i, 1) = DirectionVectorY[i];
			CurrentCS(i, 2) = DirectionVectorZ[i];
		}

		rRotationMatrix = ZeroMatrix(local_size, local_size);
		if (rRotationMatrix.size1() != local_size) {
			rRotationMatrix.resize(local_size, local_size, false);
		}
		//Building the rotation matrix for the local element matrix
		for (unsigned int kk = 0; kk < local_size; kk += dimension)
		{
			for (int i = 0; i<dimension; ++i)
			{
				for (int j = 0; j<dimension; ++j)
				{
					rRotationMatrix(i + kk, j + kk) = CurrentCS(i, j);
				}
			}
		}
		KRATOS_CATCH("")

	}	



	void TrussAdjointElement3D2N::AddExplicitContribution(const VectorType& rRHSVector,
		const Variable<VectorType>& rRHSVariable,
		Variable<array_1d<double, 3> >& rDestinationVariable,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		
		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
		{

			for (int i = 0; i< number_of_nodes; ++i)
			{
				int index = dimension * i;

				GetGeometry()[i].SetLock();

				array_1d<double, 3 > &ForceResidual = 
					GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

				for (int j = 0; j<dimension; ++j)
				{
					ForceResidual[j] += rRHSVector[index + j];
				}

				GetGeometry()[i].UnSetLock();
			}
		}
		KRATOS_CATCH("")
	}

    */
    int  TrussAdjointElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

            if (this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().PointsNumber() != 2)
            {
                KRATOS_THROW_ERROR(std::invalid_argument,
                    "The truss element works only in 3D and with 2 noded elements", "")
            }
        //verify that the variables are correctly initialized
        if (VELOCITY.Key() == 0)
            KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << "" << std::endl;
        if (DISPLACEMENT.Key() == 0)
            KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << "" << std::endl;
        if (ACCELERATION.Key() == 0)
            KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << "" << std::endl;
        if (DENSITY.Key() == 0)
            KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered" << "" << std::endl;
        if (CROSS_AREA.Key() == 0)
            KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is correctly registered" << "" << std::endl;
        //verify that the dofs exist
        for (unsigned int i = 0; i<this->GetGeometry().PointsNumber(); ++i)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            if (this->GetGeometry()[i].HasDofFor(ADJOINT_DISPLACEMENT_X) == false || 
                this->GetGeometry()[i].HasDofFor(ADJOINT_DISPLACEMENT_Y) == false || 
                this->GetGeometry()[i].HasDofFor(ADJOINT_DISPLACEMENT_Z) == false)
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
        }



        if (this->GetProperties().Has(CROSS_AREA) == false ||
            this->GetProperties()[CROSS_AREA] == 0)
        {
            KRATOS_ERROR << "CROSS_AREA not provided for this element" << this->Id() << std::endl;
        }

        if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
            this->GetProperties()[YOUNG_MODULUS] == 0)
        {
            KRATOS_ERROR << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
        }
        if (this->GetProperties().Has(DENSITY) == false)
        {
            KRATOS_ERROR << "DENSITY not provided for this element" << this->Id() << std::endl;
        }

        return 0;

        KRATOS_CATCH("")
    }

	void TrussAdjointElement3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
		rSerializer.save("mIscompressed", this->mIsCompressed);
		rSerializer.save("LinerEle", this->mIsLinearElement);

	}
	void TrussAdjointElement3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
		rSerializer.load("mIscompressed", this->mIsCompressed);
		rSerializer.load("LinerEle", this->mIsLinearElement);
	}
} // namespace Kratos.



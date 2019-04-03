//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "q4ristab_element.hpp"
#include "multiscale_application_variables.h"

#include <string>
#include <iomanip>

#define Q4RI_STAB_Q_RLX 0.0

namespace Kratos
{

	namespace Utilities
	{

		template<class TVec>
		inline void ShapeFunc(double xi, double eta, TVec & N)
		{
			N(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
			N(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
			N(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
			N(3) = 0.25 * (1.0 - xi) * (1.0 + eta);
		}

		template<class TMat>
		inline void ShapeFunc_NaturalDerivatives(double xi, double eta, TMat & dN)
		{
			dN(0, 0) = -(1.0 - eta) * 0.25;
			dN(1, 0) =  (1.0 - eta) * 0.25;
			dN(2, 0) =  (1.0 + eta) * 0.25;
			dN(3, 0) = -(1.0 + eta) * 0.25;

			dN(0, 1) = -(1.0 - xi)  * 0.25;
			dN(1, 1) = -(1.0 + xi)  * 0.25;
			dN(2, 1) =  (1.0 + xi)  * 0.25;
			dN(3, 1) =  (1.0 - xi)  * 0.25;
		}

		template<class TMat>
		inline void Jacobian(const Element::GeometryType& geom,
							 double xi, double eta, const TMat& dN,
							 TMat& J, double& detJ)
		{
			const Element::NodeType& p1 = geom[0];
			const Element::NodeType& p2 = geom[1];
			const Element::NodeType& p3 = geom[2];
			const Element::NodeType& p4 = geom[3];

			J(0, 0) = dN(0, 0) * p1.X0() + dN(1, 0) * p2.X0() + dN(2, 0) * p3.X0() + dN(3, 0) * p4.X0();
			J(0, 1) = dN(0, 0) * p1.Y0() + dN(1, 0) * p2.Y0() + dN(2, 0) * p3.Y0() + dN(3, 0) * p4.Y0();
			J(1, 0) = dN(0, 1) * p1.X0() + dN(1, 1) * p2.X0() + dN(2, 1) * p3.X0() + dN(3, 1) * p4.X0();
			J(1, 1) = dN(0, 1) * p1.Y0() + dN(1, 1) * p2.Y0() + dN(2, 1) * p3.Y0() + dN(3, 1) * p4.Y0();

			detJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
		}

		template<class TMat>
		inline void Jacobian(const Element::GeometryType& geom,
							 double xi, double eta, const TMat& dN,
							 TMat& J, double& detJ, TMat& Jinv)
		{
			const Element::NodeType& p1 = geom[0];
			const Element::NodeType& p2 = geom[1];
			const Element::NodeType& p3 = geom[2];
			const Element::NodeType& p4 = geom[3];

			J(0, 0) = dN(0, 0) * p1.X0() + dN(1, 0) * p2.X0() + dN(2, 0) * p3.X0() + dN(3, 0) * p4.X0();
			J(0, 1) = dN(0, 0) * p1.Y0() + dN(1, 0) * p2.Y0() + dN(2, 0) * p3.Y0() + dN(3, 0) * p4.Y0();
			J(1, 0) = dN(0, 1) * p1.X0() + dN(1, 1) * p2.X0() + dN(2, 1) * p3.X0() + dN(3, 1) * p4.X0();
			J(1, 1) = dN(0, 1) * p1.Y0() + dN(1, 1) * p2.Y0() + dN(2, 1) * p3.Y0() + dN(3, 1) * p4.Y0();

			detJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
			double mult = 1.0 / detJ;

			Jinv(0, 0) =   J(1, 1) * mult;
			Jinv(0, 1) = - J(0, 1) * mult;
			Jinv(1, 0) = - J(1, 0) * mult;
			Jinv(1, 1) =   J(0, 0) * mult;
		}

		template<class TMat1, class TMat2>
		inline void BMatrix(const TMat1& dNxy, TMat2& B)
		{
			B(0,  0) =  dNxy(0, 0);   B(0, 2) =  dNxy(1, 0);   B(0, 4) =  dNxy(2, 0);   B(0, 6) =  dNxy(3, 0);
			B(1,  1) =  dNxy(0, 1);   B(1, 3) =  dNxy(1, 1);   B(1, 5) =  dNxy(2, 1);   B(1, 7) =  dNxy(3, 1);
			B(2,  0) =  dNxy(0, 1);   B(2, 2) =  dNxy(1, 1);   B(2, 4) =  dNxy(2, 1);   B(2, 6) =  dNxy(3, 1);
			B(2,  1) =  dNxy(0, 0);   B(2, 3) =  dNxy(1, 0);   B(2, 5) =  dNxy(2, 0);   B(2, 7) =  dNxy(3, 0);
		}


		inline double detQ(const Matrix& C, double a)
		{
			double sinA = std::sin(a);
			double cosA = std::cos(a);
			double a1 = sinA*sinA;
			double a2 = cosA*cosA;
			double a3 = C(1, 2);
			double a4 = C(2, 0);
			double a5 = C(2, 1);
			double a6 = C(0, 2);
			double a7 = C(1, 0);
			double a8 = C(0, 1);
			double a9 = C(1, 1);
			double a10 = C(0, 0);
			double a11 = C(2, 2);
			double a12 = a2*a2;
			double a13 = a1*a1;
			double detQ =
					(a1*a4*a8*cosA + a2*a3*a8*cosA - a1*a3*a10*cosA +
					a1*a6*a7*cosA + a2*a5*a7*cosA - a2*a4*a9*cosA -
					a1*a5*a10*cosA - a2*a6*a9*cosA)*sinA +
					a9*a11*a12 - a4*a6*a13 - a3*a5*a12 + a10*a11*a13 +
					a1*a2*a3*a4 + a1*a2*a5*a6 - a1*a2*a7*a8 -
					a1*a2*a7*a11 - a1*a2*a8*a11 + a1*a2*a9*a10;
			return detQ;
		}

		inline double min_detQ(const Matrix& C)
		{
			double dq = std::numeric_limits<double>::max();
			for(int i = 0; i <= 180.0; i++)
			{
				double a = double(i)/180.0*Globals::Pi;
				dq = std::min(dq, detQ(C,a));
			}
			return dq;
		}

	}

    // =====================================================================================
    //
    // Class Q4RIStabElement
    //
    // =====================================================================================

    Q4RIStabElement::Q4RIStabElement(IndexType NewId,
									 GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
		m_H = ZeroMatrix(2,2);
		m_H_calculated = false;
    }

    Q4RIStabElement::Q4RIStabElement(IndexType NewId,
									 GeometryType::Pointer pGeometry,
									 PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
		m_H = ZeroMatrix(2,2);
		m_H_calculated = false;
    }

    Q4RIStabElement::~Q4RIStabElement()
    {
    }

    Element::Pointer Q4RIStabElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new Q4RIStabElement(NewId, newGeom, pProperties) );
    }

    Q4RIStabElement::IntegrationMethod Q4RIStabElement::GetIntegrationMethod() const
    {
        return GeometryData::GI_GAUSS_1;
    }

    void Q4RIStabElement::Initialize()
    {
		KRATOS_TRY

		const GeometryType& geom = this->GetGeometry();
		const PropertiesType& props = this->GetProperties();

		if(mConstitutiveLaw) // check whether the constitutive law pointers have been already set up
		{
			mConstitutiveLaw->InitializeMaterial( props, geom, row( geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1 ), 0 ) );
		}
		else
		{
			const ConstitutiveLaw::Pointer& newConstitutiveLaw = props[CONSTITUTIVE_LAW];
			if(newConstitutiveLaw)
			{
				mConstitutiveLaw = newConstitutiveLaw->Clone();
				mConstitutiveLaw->InitializeMaterial( props, geom, row( geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1 ), 0 ) );
			}
			else
			{
				KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
			}
		}

		m_H_calculated = false;
		m_dQn = 1.0;
		m_dQ_trial = 1.0;

		KRATOS_CATCH( "" )
    }

    void Q4RIStabElement::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        const GeometryType & geom = GetGeometry();
        const Properties& props = GetProperties();
		mConstitutiveLaw->ResetMaterial( props, geom, row( geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1 ), 0 ) );

        KRATOS_CATCH("")
    }

    void Q4RIStabElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        if(rResult.size() != 8)
            rResult.resize(8, false);

        GeometryType & geom = this->GetGeometry();

        for(SizeType i = 0; i < geom.size(); i++)
        {
            int index = i * 2;
            NodeType & iNode = geom[i];

            rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void Q4RIStabElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        ElementalDofList.reserve(8);

        GeometryType & geom = this->GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
        }
    }

    int Q4RIStabElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered","");
        if(THICKNESS.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
        if(CONSTITUTIVE_LAW.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");

        // verify that the dofs exist
        for(unsigned int i=0; i<geom.size(); i++)
        {
            if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
            if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
            if(geom[i].GetBufferSize() < 2)
                KRATOS_THROW_ERROR(std::logic_error, "This Element needs at least a buffer size = 2", "");
        }

        // check properties
        if(this->pGetProperties() == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

        const PropertiesType & props = this->GetProperties();

        if(!props.Has(CONSTITUTIVE_LAW))
			KRATOS_THROW_ERROR(std::logic_error,"CONSTITUTIVE_LAW not provided for element ",this->Id());
		const ConstitutiveLaw::Pointer & claw = props[CONSTITUTIVE_LAW];
        if(claw == NULL)
            KRATOS_THROW_ERROR(std::logic_error,"CONSTITUTIVE_LAW not provided for element ",this->Id());
        claw->Check(props, geom, rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("")
    }

    void Q4RIStabElement::CleanMemory()
    {
    }

    void Q4RIStabElement::GetValuesVector(Vector& values, int Step)
    {
        if(values.size() != 8)
            values.resize(8, false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
			const NodeType & iNode = geom[i];
            int index = i*2;
            values[index]     = iNode.FastGetSolutionStepValue(DISPLACEMENT_X, Step);
            values[index + 1] = iNode.FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
        }
    }

    void Q4RIStabElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != 8)
            values.resize(8,false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            int index = i * 2;
            values[index]        = iNode.FastGetSolutionStepValue(VELOCITY_X, Step);
            values[index + 1]    = iNode.FastGetSolutionStepValue(VELOCITY_Y, Step);
        }
    }

    void Q4RIStabElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != 8)
            values.resize(8,false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            int index = i * 2;
            values[index]        = iNode.FastGetSolutionStepValue(ACCELERATION_X, Step);
            values[index + 1]    = iNode.FastGetSolutionStepValue(ACCELERATION_Y, Step);
        }
    }

    void Q4RIStabElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
	}

    void Q4RIStabElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void Q4RIStabElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        const PropertiesType& props = GetProperties();
        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
		mConstitutiveLaw->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, 0), CurrentProcessInfo);
		if(!m_H_calculated)
		{
			this->Precompute(CurrentProcessInfo);
			m_H_calculated = true;
		}
    }

    void Q4RIStabElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        const PropertiesType& props = GetProperties();
        const GeometryType& geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
		mConstitutiveLaw->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, 0), CurrentProcessInfo);

		m_dQn = m_dQ_trial;
    }

    void Q4RIStabElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rMassMatrix.size1() != 8) || (rMassMatrix.size2() != 8))
            rMassMatrix.resize(8, 8, false);
        noalias(rMassMatrix) = ZeroMatrix(8, 8);

        // Compute the total mass

		double total_mass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];
        double lumped_mass = total_mass / 4.0;

		// loop on nodes
        for(size_t i = 0; i < 4; i++)
        {
            size_t index = i * 4;
            rMassMatrix(index, index)            = lumped_mass;
            rMassMatrix(index + 1, index + 1)    = lumped_mass;
        }
    }

    void Q4RIStabElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rDampingMatrix.size1() != 8) || (rDampingMatrix.size2() != 8))
            rDampingMatrix.resize(8, 8, false);

        noalias( rDampingMatrix ) = ZeroMatrix(8, 8);
    }

    void Q4RIStabElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void Q4RIStabElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    // =====================================================================================
    //
    // Class Q4RIStabElement - Results on Gauss Points
    //
    // =====================================================================================

	void Q4RIStabElement::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
													  std::vector<double>& rValues,
													  const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == 1)
			mConstitutiveLaw->SetValue(rVariable, rValues[0], rCurrentProcessInfo);
	}

    void Q4RIStabElement::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
													  std::vector<Vector>& rValues,
													  const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == 1)
			mConstitutiveLaw->SetValue(rVariable, rValues[0], rCurrentProcessInfo);
	}

    void Q4RIStabElement::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
													  std::vector<Matrix>& rValues,
													  const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == 1)
			mConstitutiveLaw->SetValue(rVariable, rValues[0], rCurrentProcessInfo);
	}

	void Q4RIStabElement::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													  std::vector<ConstitutiveLaw::Pointer>& rValues,
													  const ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY
		if(rValues.size() != 1)
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", rValues.size() );
		mConstitutiveLaw = rValues[0];
		if(mConstitutiveLaw == NULL)
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law is NULL ", "" );
		KRATOS_CATCH("")
	}

    void Q4RIStabElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                      std::vector<double>& rValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
    {
        if(rValues.size() != 1) rValues.resize(1);
		if(rVariable == RI_STABILIZATION)
		{
			if(m_dQ0 > 0.0)
			{
				double alpha = m_dQn/m_dQ0;
				if(alpha > -0.1)
					alpha = 1.0; //alpha max
				else
					alpha = 0.001; // alpha min
				rValues[0] = alpha;
				return;
			}
		}
        mConstitutiveLaw->GetValue(rVariable, rValues[0]);
    }

    void Q4RIStabElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                      std::vector<Vector>& rValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
    {
		if(rValues.size() != 1) rValues.resize(1);
		mConstitutiveLaw->GetValue(rVariable, rValues[0]);
    }

    void Q4RIStabElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                      std::vector<Matrix>& rValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() != 1) rValues.resize(1);

		if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||
		   rVariable == PK2_STRESS_TENSOR ||
		   rVariable == CAUCHY_STRESS_TENSOR)
		{

			Matrix& output = rValues[0];
			if(output.size1() != 2 || output.size2() != 2)
				output.resize(2,2,false);

			const GeometryType&   geom  = GetGeometry();
			const PropertiesType& props = GetProperties();

			// geometric data

			double x1 = geom[0].X0();
			double x2 = geom[1].X0();
			double x3 = geom[2].X0();
			double x4 = geom[3].X0();
			double y1 = geom[0].Y0();
			double y2 = geom[1].Y0();
			double y3 = geom[2].Y0();
			double y4 = geom[3].Y0();

			double J0 = ((x1-x3)*(y2-y4) - (x2-x4)*(y1-y3))/ 8.0;

			// strain displacement matrices

			Matrix dNxy(4,2);
			dNxy(0,0) =  y2-y4; dNxy(0,1) = -x2+x4;
			dNxy(1,0) = -y1+y3; dNxy(1,1) =  x1-x3;
			dNxy(2,0) = -y2+y4; dNxy(2,1) =  x2-x4;
			dNxy(3,0) =  y1-y3; dNxy(3,1) = -x1+x3;
			dNxy /= (8.0*J0);

			Matrix B(3,8,0.0);
			for (unsigned int i = 0; i < 4; i++)
			{
				B(0,2*i)    = dNxy(i,0);
				B(1,2*i+1)  = dNxy(i,1);
				B(2,2*i)    = dNxy(i,1);
				B(2,2*i+1)  = dNxy(i,0);
			}

			// get the current displacement vector and calculate the strain vector

			Vector U(8);
			this->GetValuesVector(U);
			Vector E( prod( B, U ) );

			if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
			{
				noalias(output) = MathUtils<double>::StrainVectorToTensor(E);
			}
			else
			{
				// compute material response

				Vector S(3, 0.0);  // stress vector
				Matrix C(3,3,0.0); // material tangent

				ConstitutiveLaw::Parameters matpar;
				matpar.SetElementGeometry( geom );
				matpar.SetMaterialProperties( props );
				matpar.SetProcessInfo( rCurrentProcessInfo );
				matpar.SetStrainVector( E );
				matpar.SetStressVector( S );
				matpar.SetConstitutiveMatrix( C );
				matpar.SetShapeFunctionsDerivatives( dNxy );
				const Matrix & shapeFunctions = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
				Vector N(4);
				for(int nodeid = 0; nodeid < 4; nodeid++) N(nodeid) = shapeFunctions(0, nodeid);
				matpar.SetShapeFunctionsValues( N );
				Matrix F(IdentityMatrix(2,2));
				double detF(1.0);
				matpar.SetDeformationGradientF(F);
				matpar.SetDeterminantF(detF);
				Flags& options = matpar.GetOptions();
				options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
				options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

				mConstitutiveLaw->CalculateMaterialResponseCauchy( matpar );

				noalias(output) = MathUtils<double>::StressVectorToTensor(S);
			}
		}
		else
		{
			mConstitutiveLaw->GetValue(rVariable, rValues[0]);
		}
	}

	void Q4RIStabElement::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													  std::vector<ConstitutiveLaw::Pointer>& rValues,
													  const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() != 1) rValues.resize(1);
		rValues[0] = mConstitutiveLaw;
	}

	void Q4RIStabElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
													   std::vector<double>& rOutput,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void Q4RIStabElement::CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
													   std::vector< array_1d<double, 3 > >& rOutput,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		//GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void Q4RIStabElement::CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
													   std::vector< Vector >& rOutput,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void Q4RIStabElement::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
													   std::vector< Matrix >& rOutput,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

    // =====================================================================================
    //
    // Class Q4RIStabElement - Private methods
    //
    // =====================================================================================

	void Q4RIStabElement::Precompute(ProcessInfo& CurrentProcessInfo)
	{
		const GeometryType&   geom  = GetGeometry();
		const PropertiesType& props = GetProperties();

		double x1 = geom[0].X0();
		double x2 = geom[1].X0();
		double x3 = geom[2].X0();
		double x4 = geom[3].X0();
		double y1 = geom[0].Y0();
		double y2 = geom[1].Y0();
		double y3 = geom[2].Y0();
		double y4 = geom[3].Y0();

		double J0 = ((x1-x3)*(y2-y4) - (x2-x4)*(y1-y3))/ 8.0;

		Matrix Jac(2,2);
		Jac(0,0) = (x2-x1+x3-x4)/4.0; Jac(0,1) = (x3-x2-x1+x4)/4.0;
		Jac(1,0) = (y2-y1+y3-y4)/4.0; Jac(1,1) = (y3-y2-y1+y4)/4.0;

		Matrix invJ(2,2);
		invJ(0,0) = -(y1+y2-y3-y4)/(4.0*J0); invJ(0,1) =  (x1+x2-x3-x4)/(4.0*J0);
		invJ(1,0) =  (y1-y2-y3+y4)/(4.0*J0); invJ(1,1) = -(x1-x2-x3+x4)/(4.0*J0);

		Matrix dNxy(4,2);
		dNxy(0,0) =  y2-y4; dNxy(0,1) = -x2+x4;
		dNxy(1,0) = -y1+y3; dNxy(1,1) =  x1-x3;
		dNxy(2,0) = -y2+y4; dNxy(2,1) =  x2-x4;
		dNxy(3,0) =  y1-y3; dNxy(3,1) = -x1+x3;
		dNxy /= (8.0*J0);

		double hx = x1 - x2 + x3 - x4;
		double hy = y1 - y2 + y3 - y4;

		array_1d<double,4> gamma;
		gamma[0] = 0.25*( 1.0 - hx*dNxy(0,0) - hy*dNxy(0,1));
		gamma[1] = 0.25*(-1.0 - hx*dNxy(1,0) - hy*dNxy(1,1));
		gamma[2] = 0.25*( 1.0 - hx*dNxy(2,0) - hy*dNxy(2,1));
		gamma[3] = 0.25*(-1.0 - hx*dNxy(3,0) - hy*dNxy(3,1));

		array_1d<double, 2> g1;
		array_1d<double, 2> g2;
		g1[0] = Jac(0,0);
		g1[1] = Jac(1,0);
		g2[0] = Jac(0,1);
		g2[1] = Jac(1,1);
		g1 /= norm_2(g1);
		g2 /= norm_2(g2);

		double th = props[THICKNESS];
		Matrix I( outer_prod(g1,g1) );
		noalias(I) += outer_prod(g2,g2);
		I *=  4.0/3.0*th*J0;

		double Hss = (I(0,0)*invJ(1,0)*invJ(1,0) + I(0,1)*invJ(0,0)*invJ(1,0) + I(1,1)*invJ(0,0)*invJ(0,0))*0.25;
		double Htt = (I(0,0)*invJ(1,1)*invJ(1,1) + I(0,1)*invJ(0,1)*invJ(1,1) + I(1,1)*invJ(0,1)*invJ(0,1))*0.25;
		double Hst = (I(0,0)*invJ(1,1)*invJ(1,0) + I(0,1)*(invJ(1,0)*invJ(0,1) + invJ(1,1)*invJ(0,0)) + I(1,1)*invJ(0,1)*invJ(0,0))*0.25;

		Matrix C0(3,3);
		Vector E(3,0.0); // note: check if initial strain is set
		Vector S(3,0.0);
		ConstitutiveLaw::Parameters matpar;
		matpar.SetElementGeometry( geom );
		matpar.SetMaterialProperties( props );
		matpar.SetProcessInfo( CurrentProcessInfo );
		matpar.SetStrainVector( E );
		matpar.SetStressVector( S );
		matpar.SetConstitutiveMatrix( C0 );
		matpar.SetShapeFunctionsDerivatives( dNxy );
		const Matrix & shapeFunctions = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
		Vector N(4);
		for(int nodeid = 0; nodeid < 4; nodeid++) N(nodeid) = shapeFunctions(0, nodeid);
		matpar.SetShapeFunctionsValues( N );
		Matrix F(IdentityMatrix(2,2));
		double detF(1.0);
		matpar.SetDeformationGradientF(F);
		matpar.SetDeterminantF(detF);
		Flags& options = matpar.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
		mConstitutiveLaw->CalculateMaterialResponseCauchy( matpar );

		m_H(0,0) = (C0(0,0) - (C0(0,1) + C0(1,0)) + C0(1,1))*Hss;
		m_H(0,1) = (C0(0,1) - (C0(0,0) + C0(1,1)) + C0(1,0))*Hst;
		m_H(1,0) = (C0(1,0) - (C0(0,0) + C0(1,1)) + C0(0,1))*Hst;
		m_H(1,1) = (C0(1,1) - (C0(0,1) + C0(1,0)) + C0(0,0))*Htt;

		m_dQ0 = Utilities::min_detQ(C0);
		m_dQn = m_dQ0;

		// override the m_dQ0 with the shear modulus
		m_dQ0 = C0(2,2);
		double stab = 1.0e-3;
		if(props.Has(RI_STABILIZATION))
			stab = props[RI_STABILIZATION];
		m_dQ0 *= stab;
	}

	void Q4RIStabElement::AddBodyForces(double V, VectorType& rRightHandSideVector)
	{
		const GeometryType& geom = GetGeometry();

		// Get shape functions
		const Matrix & N = GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

		double rho = GetProperties()[DENSITY];

		// gauss loop to integrate the external force vector (1-gp)
		// interpolate nodal volume accelerations to this gauss point
		// and obtain the body force vector
		array_1d<double, 3> bf;
		for(unsigned int inode = 0; inode < 4; inode++)
		{
			if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
				bf += N(0,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
		}
		bf *= (rho * V);

		// add it to the RHS vector
		for(unsigned int inode = 0; inode < 4; inode++)
		{
			unsigned int index = inode*2;
			double iN = N(0,inode);
			rRightHandSideVector[index + 0] += iN * bf[0];
			rRightHandSideVector[index + 1] += iN * bf[1];
		}
	}

    void Q4RIStabElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       ProcessInfo& rCurrentProcessInfo,
                                       const bool LHSrequired,
                                       const bool RHSrequired)
    {

		const GeometryType&   geom  = GetGeometry();
		const PropertiesType& props = GetProperties();

        // Resize the Left Hand Side if necessary,
        // and initialize it to Zero

        if((rLeftHandSideMatrix.size1() != 8) || (rLeftHandSideMatrix.size2() != 8))
            rLeftHandSideMatrix.resize(8, 8, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(8, 8);

        // Resize the Right Hand Side if necessary,
        // and initialize it to Zero

        if(rRightHandSideVector.size() != 8)
            rRightHandSideVector.resize(8, false);
        noalias(rRightHandSideVector) = ZeroVector(8);

        // geometric data

		double x1 = geom[0].X0();
		double x2 = geom[1].X0();
		double x3 = geom[2].X0();
		double x4 = geom[3].X0();
		double y1 = geom[0].Y0();
		double y2 = geom[1].Y0();
		double y3 = geom[2].Y0();
		double y4 = geom[3].Y0();

		double J0 = ((x1-x3)*(y2-y4) - (x2-x4)*(y1-y3))/ 8.0;
		double th = props[THICKNESS];
		double V  = 4.0*J0*th;

		// strain displacement matrices

		Matrix dNxy(4,2);
		dNxy(0,0) =  y2-y4; dNxy(0,1) = -x2+x4;
		dNxy(1,0) = -y1+y3; dNxy(1,1) =  x1-x3;
		dNxy(2,0) = -y2+y4; dNxy(2,1) =  x2-x4;
		dNxy(3,0) =  y1-y3; dNxy(3,1) = -x1+x3;
		dNxy /= (8.0*J0);

		double hx = x1 - x2 + x3 - x4;
		double hy = y1 - y2 + y3 - y4;

		array_1d<double,4> gamma;
		gamma[0] = 0.25*( 1.0 - hx*dNxy(0,0) - hy*dNxy(0,1));
		gamma[1] = 0.25*(-1.0 - hx*dNxy(1,0) - hy*dNxy(1,1));
		gamma[2] = 0.25*( 1.0 - hx*dNxy(2,0) - hy*dNxy(2,1));
		gamma[3] = 0.25*(-1.0 - hx*dNxy(3,0) - hy*dNxy(3,1));

		Matrix B(3,8,0.0);
		Matrix Bh(2,8,0.0);

		for (unsigned int i = 0; i < 4; i++)
		{
			B(0,2*i)    = dNxy(i,0);
			B(1,2*i+1)  = dNxy(i,1);
			B(2,2*i)    = dNxy(i,1);
			B(2,2*i+1)  = dNxy(i,0);

			Bh(0,2*i)   = gamma[i];
			Bh(1,2*i+1) = gamma[i];
		}

		// get the current displacement vector and calculate the strain vector

		Vector U(8);
		this->GetValuesVector(U);
		Vector E( prod( B, U ) );

		// compute material response

		Vector S(3, 0.0);  // stress vector
		Matrix C(3,3,0.0); // material constitutive tensor

		ConstitutiveLaw::Parameters matpar;
		matpar.SetElementGeometry( geom );
		matpar.SetMaterialProperties( props );
		matpar.SetProcessInfo( rCurrentProcessInfo );
		matpar.SetStrainVector( E );
		matpar.SetStressVector( S );
		matpar.SetConstitutiveMatrix( C );
		matpar.SetShapeFunctionsDerivatives( dNxy );
		const Matrix & shapeFunctions = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
		Vector N(4);
		for(int nodeid = 0; nodeid < 4; nodeid++) N(nodeid) = shapeFunctions(0, nodeid);
		matpar.SetShapeFunctionsValues( N );
		Matrix F(IdentityMatrix(2,2));
		double detF(1.0);
		matpar.SetDeformationGradientF(F);
		matpar.SetDeterminantF(detF);
		Flags& options = matpar.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);
		mConstitutiveLaw->CalculateMaterialResponseCauchy( matpar );

		//if(LHSrequired)
		//{
		//  m_dQ_trial = Utilities::min_detQ(C);
		//	/*if(m_dQ_trial > m_dQn)
		//		m_dQ_trial = m_dQn;*/
		//	if(m_dQ_trial > m_dQ0)
		//		m_dQ_trial = m_dQ0;
		//	if(std::abs(m_dQ_trial-m_dQn) < 1.0E-5*m_dQ0) {
		//		m_dQ_trial = m_dQn;
		//	}
		//}
		//double alpha = m_dQn/m_dQ0;
		//if(alpha > -0.1)
		//	alpha = 1.0; //alpha max
		//else
		//	alpha = 0.001; // alpha min
		double alpha = m_dQ0;
		if(LHSrequired) {
			m_dQ_trial = Utilities::min_detQ(C);
			if(m_dQ_trial > m_dQn) m_dQ_trial = m_dQn;
		}
		if(m_dQn < 0.0 && props.Has(RI_STABILIZATION_RESIDUAL))
			alpha *= props[RI_STABILIZATION_RESIDUAL];

		// internal forces and tangent

		Matrix BhTH( prod(trans(Bh), m_H*alpha) );
		noalias( rLeftHandSideMatrix ) = prod( BhTH, Bh ); // now LHS contains Kh(stabilization)

		noalias( rRightHandSideVector ) = -prod( rLeftHandSideMatrix, U ); // stabilizing hourglass forces
		noalias( rRightHandSideVector ) -= prod( trans( B ), S * V ); // physical force vector

		if(LHSrequired)
		{
			Matrix BTD( prod( trans( B ), C * V ) );
			noalias( rLeftHandSideMatrix ) += prod( BTD, B );
		}

        // Add body forces contributions.

        AddBodyForces(V, rRightHandSideVector);
    }

    // =====================================================================================
    //
    // Class Q4RIStabElement - Serialization
    //
    // =====================================================================================

    void Q4RIStabElement::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLaw);
    }

    void Q4RIStabElement::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLaw);
    }

}

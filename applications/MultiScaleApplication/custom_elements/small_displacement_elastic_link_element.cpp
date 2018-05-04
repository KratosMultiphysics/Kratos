//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "small_displacement_elastic_link_element.hpp"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"
#include <string>
#include <iomanip>

namespace Kratos
{

    SmallDisplacementElasticLinkElement::SmallDisplacementElasticLinkElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    SmallDisplacementElasticLinkElement::SmallDisplacementElasticLinkElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    SmallDisplacementElasticLinkElement::~SmallDisplacementElasticLinkElement()
    {
    }

    Element::Pointer SmallDisplacementElasticLinkElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new SmallDisplacementElasticLinkElement(NewId, newGeom, pProperties) );
    }

    void SmallDisplacementElasticLinkElement::Initialize()
    {
    }

    void SmallDisplacementElasticLinkElement::ResetConstitutiveLaw()
    {
    }

    void SmallDisplacementElasticLinkElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        GeometryType & geom = this->GetGeometry();
        if(rResult.size() != 6)
            rResult.resize(6, false);

        for(SizeType i = 0; i < 2; i++)
        {
			unsigned int index = i*3;
            NodeType & iNode = geom[i];
            rResult[index  ] = iNode.GetDof(DISPLACEMENT_X).EquationId();
			rResult[index+1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index+2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void SmallDisplacementElasticLinkElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        ElementalDofList.reserve(6);
		GeometryType& geom = GetGeometry();
        for (SizeType i = 0; i < 2; i++)
        {
            NodeType & iNode = geom[i];
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));
        }
    }

    int SmallDisplacementElasticLinkElement::Check(const ProcessInfo& rCurrentProcessInfo)
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
        if(YOUNG_MODULUS.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero! (check if the application is correctly registered","");

        // verify that the dofs exist
        for(unsigned int i=0; i<geom.PointsNumber(); i++)
        {
            if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
            if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false || geom[i].HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
        }

        // check properties
        if(this->pGetProperties() == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

        const PropertiesType & props = this->GetProperties();

        if(!props.Has(YOUNG_MODULUS))
			KRATOS_THROW_ERROR(std::logic_error,"YOUNG_MODULUS not provided for element ",this->Id());

        return 0;

        KRATOS_CATCH("")
    }

    void SmallDisplacementElasticLinkElement::CleanMemory()
    {
    }

    void SmallDisplacementElasticLinkElement::GetValuesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();
        if(values.size() != 6)
            values.resize(6,false);
        for (SizeType i = 0; i < 2; i++)
        {
			unsigned int index = i*3;
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& val = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
            values[index  ] = val[0];
			values[index+1] = val[1];
			values[index+2] = val[2];
        }
    }

    void SmallDisplacementElasticLinkElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();
        if(values.size() != 6)
            values.resize(6,false);
        for (SizeType i = 0; i < 2; i++)
        {
			unsigned int index = i*3;
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& val = iNode.FastGetSolutionStepValue(VELOCITY, Step);
            values[index  ] = val[0];
			values[index+1] = val[1];
			values[index+2] = val[2];
        }
    }

    void SmallDisplacementElasticLinkElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();
        if(values.size() != 6)
            values.resize(6,false);
        for (SizeType i = 0; i < 2; i++)
        {
			unsigned int index = i*3;
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& val = iNode.FastGetSolutionStepValue(ACCELERATION, Step);
            values[index  ] = val[0];
			values[index+1] = val[1];
			values[index+2] = val[2];
        }
    }

    void SmallDisplacementElasticLinkElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rMassMatrix.size1() != 6) || (rMassMatrix.size2() != 6))
            rMassMatrix.resize(6, 6, false);
        noalias(rMassMatrix) = ZeroMatrix(6, 6);
    }

    void SmallDisplacementElasticLinkElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rDampingMatrix.size1() != 6) || (rDampingMatrix.size2() != 6))
            rDampingMatrix.resize(6, 6, false);
        noalias( rDampingMatrix ) = ZeroMatrix(6, 6);
    }

    void SmallDisplacementElasticLinkElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     ProcessInfo& rCurrentProcessInfo)
    {
		if((rLeftHandSideMatrix.size1() != 6) || (rLeftHandSideMatrix.size2() != 6))
			rLeftHandSideMatrix.resize(6, 6, false);

		if( rRightHandSideVector.size() != 6 )
			rRightHandSideVector.resize(6, false);

		double K = GetProperties()[YOUNG_MODULUS];
		rLeftHandSideMatrix.clear();
		rLeftHandSideMatrix(0,0) =  K;
		rLeftHandSideMatrix(3,3) =  K;
		rLeftHandSideMatrix(0,3) = -K;
		rLeftHandSideMatrix(3,0) = -K;

		array_1d<double,3> vx,vy,vz;
		double L;
		CalculateLocalAxes(vx,vy,vz,L);

		Matrix R(6,6,0.0);
		for(unsigned int i=0; i<3; i++)
		{
			R(0,i) = vx(i);
			R(1,i) = vy(i);
			R(2,i) = vz(i);

			R(3,i+3) = vx(i);
			R(4,i+3) = vy(i);
			R(5,i+3) = vz(i);
		}

		Matrix RTK( prod(trans(R), rLeftHandSideMatrix) );
		noalias(rLeftHandSideMatrix) = prod(RTK, R);

		Vector U(6);
		GetValuesVector(U);
		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, U);
    }

    void SmallDisplacementElasticLinkElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
    {
		Matrix rLeftHandSideMatrix;
		this->CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }

    // =====================================================================================
    //
    // Class SmallDisplacementElasticLinkElement - Results on Gauss Points
    //
    // =====================================================================================

	void SmallDisplacementElasticLinkElement::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
													   std::vector<double>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
	}

    void SmallDisplacementElasticLinkElement::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
													   std::vector<Vector>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
	}

    void SmallDisplacementElasticLinkElement::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
													   std::vector<Matrix>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void SmallDisplacementElasticLinkElement::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													 std::vector<ConstitutiveLaw::Pointer>& rValues,
													 const ProcessInfo& rCurrentProcessInfo )
	{
	}

    void SmallDisplacementElasticLinkElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                     std::vector<double>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                     std::vector<Vector>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void SmallDisplacementElasticLinkElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                     std::vector<Matrix>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void SmallDisplacementElasticLinkElement::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													 std::vector<ConstitutiveLaw::Pointer>& rValues,
													 const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void SmallDisplacementElasticLinkElement::CalculateLocalAxes(array_1d<double,3>& vx, array_1d<double,3>& vy, array_1d<double,3>& vz, double& L)
	{
		GeometryType& geom = GetGeometry();
		vx[0] = geom[1].X0() - geom[0].X0();
		vx[1] = geom[1].Y0() - geom[0].Y0();
		vx[2] = geom[1].Z0() - geom[0].Z0();
		L = norm_2(vx);
		if(L > 0.0)
		{
			vx /= L;
			if(std::abs(vx[2]) > 1.0-1.0e-7)
			{
				// aligned with global axes
				vx.clear();
				vx[2] = 1.0; // local x = global Z
				vy.clear();
				vy[1] = 1.0; // local y = global Y
				vz.clear();
				vz[1] = -1.0; // local z = - global X
			}
			else
			{
				array_1d<double, 3> global_z;
				global_z.clear();
				global_z[2] = 1.0;
				MathUtils<double>::CrossProduct(vy,  global_z, vx);
				MathUtils<double>::CrossProduct(vz,   vx, vy);
				double ly = norm_2(vy);
				double lz = norm_2(vz);
				vy /= ly;
				vz /= lz;
			}
		}
		else
		{
			// coinciding nodes. assume local axes alligned with global axes
			vx.clear();
			vx[0] = 1.0;
			vy.clear();
			vy[1] = 1.0;
			vz.clear();
			vz[2] = 1.0;
			L = 1.0;
		}
	}

	// =====================================================================================
    //
    // Class SmallDisplacementElasticLinkElement - Serialization
    //
    // =====================================================================================

	void SmallDisplacementElasticLinkElement::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    void SmallDisplacementElasticLinkElement::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

}

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "opt_triangle_element.hpp"
#include "multiscale_application_variables.h"

#include "geometries/triangle_2d_3.h"

#include <string>
#include <iomanip>

#define OPT_NUM_NODES 3
#define OPT_STRAIN_SIZE 3
#define OPT_NUM_DOFS 9

//----------------------------------------
// preprocessors for the integration
// method used by this element.

//#define OPT_1_POINT_INTEGRATION

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_1
#define OPT_NUM_GP 1
#else
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_2
#define OPT_NUM_GP 3
#endif // OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

//#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) Utilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

#define OPT_AVARAGE_RESULTS

namespace Kratos
{

	namespace Utilities
	{
		inline void InterpToStandardGaussPoints(double& v1, double& v2, double& v3)
		{
			double vg1 = v1;
			double vg2 = v2;
			double vg3 = v3;
#ifdef OPT_AVARAGE_RESULTS
			v1 = (vg1+vg2+vg3)/3.0;
			v2 = (vg1+vg2+vg3)/3.0;
			v3 = (vg1+vg2+vg3)/3.0;
#else
			v1 = (2.0*vg1)/3.0 - vg2/3.0       + (2.0*vg3)/3.0;
			v2 = (2.0*vg1)/3.0 + (2.0*vg2)/3.0 - vg3/3.0;
			v3 = (2.0*vg2)/3.0 - vg1/3.0       + (2.0*vg3)/3.0;
#endif // OPT_AVARAGE_RESULTS
		}

		inline void InterpToStandardGaussPoints(std::vector< double >& v)
		{
			if(v.size() != 3) return;
			InterpToStandardGaussPoints(v[0], v[1], v[2]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double,3> >& v)
		{
			if(v.size() != 3) return;
			for(size_t i = 0; i < 3; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double,6> >& v)
		{
			if(v.size() != 3) return;
			for(size_t i = 0; i < 6; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Vector >& v)
		{
			if(v.size() != 3) return;
			size_t ncomp = v[0].size();
			for(int i = 1; i < 3; i++)
				if(v[i].size() != ncomp)
					return;
			for(size_t i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Matrix >& v)
		{
			if(v.size() != 3) return;
			size_t nrows = v[0].size1();
			size_t ncols = v[0].size2();
			for(int i = 1; i < 3; i++)
				if(v[i].size1() != nrows || v[i].size2() != ncols)
					return;
			for(size_t i = 0; i < nrows; i++)
				for(size_t j = 0; j < ncols; j++)
					InterpToStandardGaussPoints(v[0](i,j), v[1](i,j), v[2](i,j));
		}

	}

    // =====================================================================================
    //
    // CalculationData
    //
    // =====================================================================================

    OptTriangleElement::CalculationData::CalculationData(const ProcessInfo& rCurrentProcessInfo)
	        : CurrentProcessInfo(rCurrentProcessInfo)

	{
	}

    // =====================================================================================
    //
    // Class OptTriangleElement
    //
    // =====================================================================================

    OptTriangleElement::OptTriangleElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
        mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
    }

    OptTriangleElement::OptTriangleElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry,
                                               PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
    }

    OptTriangleElement::~OptTriangleElement()
    {
    }

    Element::Pointer OptTriangleElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new OptTriangleElement(NewId, newGeom, pProperties) );
    }

    OptTriangleElement::IntegrationMethod OptTriangleElement::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    void OptTriangleElement::Initialize()
    {
        KRATOS_TRY

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		//Constitutive Law initialization
		if ( mConstitutiveLawVector.size() != integration_points.size() )
		{
			mConstitutiveLawVector.resize( integration_points.size() );
		}
		else
		{
			// check whether the constitutive law pointers have been already set up
			bool already_set_up = true;
			for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				if(mConstitutiveLawVector[i] == NULL)
					already_set_up = false;
			}
			if(already_set_up)
			{
				for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
				{
					mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
							row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
				}

				return; // if so, we are done here!
			}
		}

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
		{
			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
						row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
			}
		}
		else
		{
			KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
		}

		KRATOS_CATCH( "" )
    }

    void OptTriangleElement::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        const Properties& props = GetProperties();
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->ResetMaterial(props, geom, row(shapeFunctionsValues, i));

        KRATOS_CATCH("")
    }

    void OptTriangleElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        if(rResult.size() != OPT_NUM_DOFS)
            rResult.resize(OPT_NUM_DOFS, false);

        GeometryType & geom = this->GetGeometry();

        for(SizeType i = 0; i < geom.size(); i++)
        {
            int index = i * 3;
            NodeType & iNode = geom[i];

            rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = iNode.GetDof(ROTATION_Z).EquationId();
        }
    }

    void OptTriangleElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        ElementalDofList.reserve(OPT_NUM_DOFS);

        GeometryType & geom = this->GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
        }
    }

    int OptTriangleElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
        if(ROTATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");
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
            if(geom[i].SolutionStepsDataHas(ROTATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ",geom[i].Id());
            if(geom[i].HasDofFor(ROTATION_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ",geom[i].Id());

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

    void OptTriangleElement::CleanMemory()
    {
    }

    void OptTriangleElement::GetValuesVector(Vector& values, int Step)
    {
        if(values.size() != OPT_NUM_DOFS)
            values.resize(OPT_NUM_DOFS, false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
			const NodeType & iNode = geom[i];
            int index = i*3;
            values[index]     = iNode.FastGetSolutionStepValue(DISPLACEMENT_X, Step);
            values[index + 1] = iNode.FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
            values[index + 2] = iNode.FastGetSolutionStepValue(ROTATION_Z, Step);
        }
    }

    void OptTriangleElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != OPT_NUM_DOFS)
            values.resize(OPT_NUM_DOFS,false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            int index = i * 3;
            values[index]        = iNode.FastGetSolutionStepValue(VELOCITY_X, Step);
            values[index + 1]    = iNode.FastGetSolutionStepValue(VELOCITY_Y, Step);
            values[index + 2]    = 0.0;
        }
    }

    void OptTriangleElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != OPT_NUM_DOFS)
            values.resize(OPT_NUM_DOFS,false);

        const GeometryType & geom = GetGeometry();

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            int index = i * 3;
            values[index]        = iNode.FastGetSolutionStepValue(ACCELERATION_X, Step);
            values[index + 1]    = iNode.FastGetSolutionStepValue(ACCELERATION_Y, Step);
            values[index + 2]    = 0.0;
        }
    }

    void OptTriangleElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        /*const GeometryType & geom = this->GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->InitializeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);*/
	}

    void OptTriangleElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        /*const GeometryType & geom = this->GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->FinalizeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);*/
    }

    void OptTriangleElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        const PropertiesType& props = GetProperties();
        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
    }

    void OptTriangleElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        const PropertiesType& props = GetProperties();
        const GeometryType& geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
    }

    void OptTriangleElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rMassMatrix.size1() != OPT_NUM_DOFS) || (rMassMatrix.size2() != OPT_NUM_DOFS))
            rMassMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
        noalias(rMassMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

        // Compute the total mass

		double total_mass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];
        double lumped_mass = total_mass / 3.0;

		// loop on nodes
        for(size_t i = 0; i < 3; i++)
        {
            size_t index = i * 3;

            // translational mass
            rMassMatrix(index, index)            = lumped_mass;
            rMassMatrix(index + 1, index + 1)    = lumped_mass;

            // rotational mass - neglected for the moment...
        }
    }

    void OptTriangleElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rDampingMatrix.size1() != OPT_NUM_DOFS) || (rDampingMatrix.size2() != OPT_NUM_DOFS))
            rDampingMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);

        noalias( rDampingMatrix ) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);
    }

    void OptTriangleElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void OptTriangleElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, false, true);
    }

    // =====================================================================================
    //
    // Class OptTriangleElement - Results on Gauss Points
    //
    // =====================================================================================

	void OptTriangleElement::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
														 std::vector<ConstitutiveLaw::Pointer>& rValues,
														 const ProcessInfo& rCurrentProcessInfo )
	{
		if ( mConstitutiveLawVector.size() != rValues.size() )
		{
			mConstitutiveLawVector.resize(rValues.size());
			if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod  ) )
				KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
		}
		for(unsigned int i=0; i<rValues.size(); i++)
			mConstitutiveLawVector[i] = rValues[i];
	}

    void OptTriangleElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                         std::vector<double>& rValues,
                                                         const ProcessInfo& rCurrentProcessInfo)
    {
        if(rValues.size() != OPT_NUM_GP)
            rValues.resize(OPT_NUM_GP);

        for(int i = 0; i < OPT_NUM_GP; i++) {
			rValues[i] = 0.0;
            mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}

		OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
    }

    void OptTriangleElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                         std::vector<Vector>& rValues,
                                                         const ProcessInfo& rCurrentProcessInfo)
    {
		if(rValues.size() != OPT_NUM_GP)
            rValues.resize(OPT_NUM_GP);

        for(int i = 0; i < OPT_NUM_GP; i++) {
            mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}

		OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
    }

    void OptTriangleElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                         std::vector<Matrix>& rValues,
                                                         const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() != OPT_NUM_GP)
			rValues.resize(OPT_NUM_GP);

		if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||
		   rVariable == PK2_STRESS_TENSOR ||
		   rVariable == CAUCHY_STRESS_TENSOR)
		{
			CalculationData data(rCurrentProcessInfo);
			data.CalculateLHS = false;
			data.CalculateRHS = true;
			InitializeCalculationData(data);

			for(size_t i = 0; i < OPT_NUM_GP; i++)
			{
				Matrix& output = rValues[i];
				if(output.size1() != 2 || output.size2() != 2)
					output.resize(2,2,false);

				data.gpIndex = i;

				// calculate beta0
				CalculateBeta0( data );

				// calculate the total strain displ. matrix
				CalculateBMatrix( data );

				// compute strain vector
				noalias( data.E ) = prod( data.B, data.U );

				if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
				{
					noalias(output) = MathUtils<double>::StrainVectorToTensor(data.E);
				}
				else
				{
					// calculate material response
					CalculateConstitutiveLawResponse( data );
					noalias(output) = MathUtils<double>::StressVectorToTensor(data.S);
				}
			}
		}
		else
		{
			for(size_t i = 0; i < OPT_NUM_GP; i++)
				mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}

		OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
	}

	void OptTriangleElement::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
														 std::vector<ConstitutiveLaw::Pointer>& rValues,
														 const ProcessInfo& rCurrentProcessInfo)
	{
		if ( rValues.size() != mConstitutiveLawVector.size() )
			rValues.resize(mConstitutiveLawVector.size());

		for(unsigned int i=0; i<rValues.size(); i++)
			rValues[i] = mConstitutiveLawVector[i];
	}


    // =====================================================================================
    //
    // Class OptTriangleElement - Private methods
    //
    // =====================================================================================

    void OptTriangleElement::DecimalCorrection(Vector& a)
    {
        double norm = norm_2(a);
        double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
        for(SizeType i = 0; i < a.size(); i++)
            if(std::abs(a(i)) < tolerance)
                a(i) = 0.0;
    }

	void OptTriangleElement::InitializeCalculationData(CalculationData& data)
	{
		//-------------------------------------
		// Computation of all stuff that remain
		// constant throughout the calculations

		//-------------------------------------
		// geometry data

		const GeometryType& geom = GetGeometry();
		const NodeType& p1 = geom[0];
		const NodeType& p2 = geom[1];
		const NodeType& p3 = geom[2];

		const double x12 = p1.X0() - p2.X0();
		const double x23 = p2.X0() - p3.X0();
		const double x31 = p3.X0() - p1.X0();
		const double x21 = -x12;
		const double x32 = -x23;
		const double x13 = -x31;

		const double y12 = p1.Y0() - p2.Y0();
		const double y23 = p2.Y0() - p3.Y0();
		const double y31 = p3.Y0() - p1.Y0();
		const double y21 = -y12;
		const double y32 = -y23;
		const double y13 = -y31;

		const double A  = 0.5*(y21*x13 - x21*y13);
		const double A2 = 2.0*A;
		const double A4 = 4.0*A;
		const double AA4 = A * A4;

		const double LL21 = x21*x21 + y21*y21;
		const double LL32 = x32*x32 + y32*y32;
		const double LL13 = x13*x13 + y13*y13;

		// Note: here we compute the avarage thickness,
		// since L is constant over the element.
		// Now it is not necessary to compute the avarage
		// because the current implementation of the cross section
		// doesn't have a variable thickness
		// (for example as a function of the spatial coordinates...).
		// This is just a place-holder for future
		// implementation of a variable thickness

		double h = GetProperties()[THICKNESS];
		double volume = A * h;

		// this is the integration weight
		// used during the gauss loop.

		data.dV = volume / (double)OPT_NUM_GP;

		// crete the integration point locations
		if(data.gpLocations.size() != 0) data.gpLocations.clear();
		data.gpLocations.resize( OPT_NUM_GP );
#ifdef OPT_1_POINT_INTEGRATION
		array_1d<double,3>& gp0 = data.gpLocations[0];
		gp0[0] = 1.0/3.0; gp0[1] = 1.0/3.0; gp0[2] = 1.0/3.0;
#else
		array_1d<double,3>& gp0 = data.gpLocations[0];
		array_1d<double,3>& gp1 = data.gpLocations[1];
		array_1d<double,3>& gp2 = data.gpLocations[2];
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
		gp0[0] = 1.0/6.0; gp0[1] = 1.0/6.0; gp0[2] = 2.0/3.0;
		gp1[0] = 2.0/3.0; gp1[1] = 1.0/6.0; gp1[2] = 1.0/6.0;
		gp2[0] = 1.0/6.0; gp2[1] = 2.0/3.0; gp2[2] = 1.0/6.0;
#else
		gp0[0] = 0.5; gp0[1] = 0.5; gp0[2] = 0.0;
		gp1[0] = 0.0; gp1[1] = 0.5; gp1[2] = 0.5;
		gp2[0] = 0.5; gp2[1] = 0.0; gp2[2] = 0.5;
#endif // OPT_USES_INTERIOR_GAUSS_POINTS

#endif // OPT_1_POINT_INTEGRATION

		// cartesian derivatives
		data.dNxy.resize(3, 2, false);
		data.dNxy(0, 0) = (y13 - y12)/A2; data.dNxy(0, 1) = (x12 - x13)/A2;
		data.dNxy(1, 0) =        -y13/A2; data.dNxy(1, 1) =         x13/A2;
		data.dNxy(2, 0) =         y12/A2; data.dNxy(2, 1) =        -x12/A2;

		//-------------------------------------
		// template parameters

		const double b1    =  1.0;
		const double b2    =  2.0;
		const double b3    =  1.0;
		const double b4    =  0.0;
		const double b5    =  1.0;
		const double b6    = -1.0;
		const double b7    = -1.0;
		const double b8    = -1.0;
		const double b9    = -2.0;

		const double alpha   =  1.5;
		const double alpha_6 =  alpha / 6.0;

		//--------------------------------------
		// calculate L - Lumping matrix
		// for the construction of the basic
		// stiffness

		double L_mult = 0.5/A;

		data.L.resize(3, 9, false);

		data.L(0, 0) = L_mult * y23;
		data.L(1, 0) = 0.00;
		data.L(2, 0) = L_mult * x32;
		data.L(0, 1) = 0.00;
		data.L(1, 1) = L_mult * x32;
		data.L(2, 1) = L_mult * y23;
		data.L(0, 2) = L_mult * y23*(y13-y21)*alpha_6;
		data.L(1, 2) = L_mult * x32*(x31-x12)*alpha_6;
		data.L(2, 2) = L_mult * 2.00*(x31*y13-x12*y21)*alpha_6;
		data.L(0, 3) = L_mult * y31;
		data.L(1, 3) = 0.00;
		data.L(2, 3) = L_mult * x13;
		data.L(0, 4) = 0.00;
		data.L(1, 4) = L_mult * x13;
		data.L(2, 4) = L_mult * y31;
		data.L(0, 5) = L_mult * y31*(y21-y32)*alpha_6;
		data.L(1, 5) = L_mult * x13*(x12-x23)*alpha_6;
		data.L(2, 5) = L_mult * 2.00*(x12*y21-x23*y32)*alpha_6;
		data.L(0, 6) = L_mult * y12;
		data.L(1, 6) = 0.00;
		data.L(2, 6) = L_mult * x21;
		data.L(0, 7) = 0.00;
		data.L(1, 7) = L_mult * x21;
		data.L(2, 7) = L_mult * y12;
		data.L(0, 8) = L_mult * y12*(y32-y13)*alpha_6;
		data.L(1, 8) = L_mult * x21*(x23-x31)*alpha_6;
		data.L(2, 8) = L_mult * 2.00*(x23*y32-x31*y13)*alpha_6;

		/*double a6 = alpha/6.0;
		double a3 = alpha/3.0;
		double y21_2 = y21*y21;
		double y31_2 = y31*y31;
		double x21_2 = x21*x21;
		double x31_2 = x31*x31;
		double y12_2 = y12*y12;
		double y32_2 = y32*y32;
		double x12_2 = x12*x12;
		double x32_2 = x32*x32;
		double y13_2 = y13*y13;
		double y23_2 = y23*y23;
		double x13_2 = x13*x13;
		double x23_2 = x23*x23;
		data.L(0, 0) = L_mult*y23;
		data.L(1, 0) = 0;
		data.L(2, 0) = -L_mult*x23;
		data.L(0, 1) = 0;
		data.L(1, 1) = -L_mult*x23;
		data.L(2, 1) = L_mult*y23;
		data.L(0, 2) = -L_mult*a6*(y21_2 - y31_2);
		data.L(1, 2) = -L_mult*a6*(x21_2 - x31_2);
		data.L(2, 2) = L_mult*a3*(x21*y21 - x31*y31);
		data.L(0, 3) = L_mult*y31;
		data.L(1, 3) = 0;
		data.L(2, 3) = -L_mult*x31;
		data.L(0, 4) = 0;
		data.L(1, 4) = -L_mult*x31;
		data.L(2, 4) = L_mult*y31;
		data.L(0, 5) = L_mult*a6*(y12_2 - y32_2);
		data.L(1, 5) = L_mult*a6*(x12_2 - x32_2);
		data.L(2, 5) = -L_mult*a3*(x12*y12 - x32*y32);
		data.L(0, 6) = L_mult*y12;
		data.L(1, 6) = 0;
		data.L(2, 6) = -L_mult*x12;
		data.L(0, 7) = 0;
		data.L(1, 7) = -L_mult*x12;
		data.L(2, 7) = L_mult*y12;
		data.L(0, 8) = -L_mult*a6*(y13_2 - y23_2);
		data.L(1, 8) = -L_mult*a6*(x13_2 - x23_2);
		data.L(2, 8) = L_mult*a3*(x13*y13 - x23*y23);*/

		//--------------------------------------
		// calculate Q1,Q2,Q3 - matrices
		// for the construction of the
		// higher order stiffness

		data.Q1.resize(3, 3, false);

		data.Q1(0,0) = b1*A2/(LL21*3.00);
		data.Q1(0,1) = b2*A2/(LL21*3.00);
		data.Q1(0,2) = b3*A2/(LL21*3.00);
		data.Q1(1,0) = b4*A2/(LL32*3.00);
		data.Q1(1,1) = b5*A2/(LL32*3.00);
		data.Q1(1,2) = b6*A2/(LL32*3.00);
		data.Q1(2,0) = b7*A2/(LL13*3.00);
		data.Q1(2,1) = b8*A2/(LL13*3.00);
		data.Q1(2,2) = b9*A2/(LL13*3.00);

		data.Q2.resize(3, 3, false);

		data.Q2(0,0) = b9*A2/(LL21*3.00);
		data.Q2(0,1) = b7*A2/(LL21*3.00);
		data.Q2(0,2) = b8*A2/(LL21*3.00);
		data.Q2(1,0) = b3*A2/(LL32*3.00);
		data.Q2(1,1) = b1*A2/(LL32*3.00);
		data.Q2(1,2) = b2*A2/(LL32*3.00);
		data.Q2(2,0) = b6*A2/(LL13*3.00);
		data.Q2(2,1) = b4*A2/(LL13*3.00);
		data.Q2(2,2) = b5*A2/(LL13*3.00);

		data.Q3.resize(3, 3, false);

		data.Q3(0,0) = b5*A2/(LL21*3.00);
		data.Q3(0,1) = b6*A2/(LL21*3.00);
		data.Q3(0,2) = b4*A2/(LL21*3.00);
		data.Q3(1,0) = b8*A2/(LL32*3.00);
		data.Q3(1,1) = b9*A2/(LL32*3.00);
		data.Q3(1,2) = b7*A2/(LL32*3.00);
		data.Q3(2,0) = b2*A2/(LL13*3.00);
		data.Q3(2,1) = b3*A2/(LL13*3.00);
		data.Q3(2,2) = b1*A2/(LL13*3.00);

		//--------------------------------------
		// calculate Te, TTu -
		// transformation matrices
		// for the construction of the
		// higher order stiffness

		data.Te.resize(3, 3, false);

		data.Te(0,0) = 1.0/AA4 * y23*y13*LL21;
		data.Te(0,1) = 1.0/AA4 * y31*y21*LL32;
		data.Te(0,2) = 1.0/AA4 * y12*y32*LL13;
		data.Te(1,0) = 1.0/AA4 * x23*x13*LL21;
		data.Te(1,1) = 1.0/AA4 * x31*x21*LL32;
		data.Te(1,2) = 1.0/AA4 * x12*x32*LL13;
		data.Te(2,0) = 1.0/AA4 * (y23*x31+x32*y13)*LL21;
		data.Te(2,1) = 1.0/AA4 * (y31*x12+x13*y21)*LL32;
		data.Te(2,2) = 1.0/AA4 * (y12*x23+x21*y32)*LL13;

		data.TTu.resize(3, 9, false);

		for(unsigned int i=0; i<3; i++)
		{
			data.TTu(i, 0) = 1.0/A4 * x32;
			data.TTu(i, 1) = 1.0/A4 * y32;
			data.TTu(i, 2) = 0.0;
			data.TTu(i, 3) = 1.0/A4 * x13;
			data.TTu(i, 4) = 1.0/A4 * y13;
			data.TTu(i, 5) = 0.0;
			data.TTu(i, 6) = 1.0/A4 * x21;
			data.TTu(i, 7) = 1.0/A4 * y21;
			data.TTu(i, 8) = 0.0;
		}
		data.TTu(0, 2) = 1.0;
		data.TTu(1, 5) = 1.0;
		data.TTu(2, 8) = 1.0;

		//--------------------------------------
		// calculate the displacement vector
		// in global and local coordinate systems

		data.U.resize(OPT_NUM_DOFS, false);
		GetValuesVector( data.U );

		//--------------------------------------
		// Finally allocate all auxiliary
		// matrices to be used later on
		// during the element integration.
		// Just to avoid re-allocations

		data.B.resize(OPT_STRAIN_SIZE, OPT_NUM_DOFS, false);
		data.D.resize(OPT_STRAIN_SIZE, OPT_STRAIN_SIZE, false);
		data.BTD.resize(OPT_NUM_DOFS, OPT_STRAIN_SIZE, false);

		data.E.resize(OPT_STRAIN_SIZE, false);
		data.S.resize(OPT_STRAIN_SIZE, false);

		data.N.resize(3, false);

		data.Q.resize(3, 3, false);
		data.Qh.resize(3, 9, false);
		data.TeQ.resize(3, 3, false);

		//--------------------------------------
		// Initialize the section parameters

		data.MaterialParameters.SetElementGeometry( GetGeometry() );
		data.MaterialParameters.SetMaterialProperties( GetProperties() );
		data.MaterialParameters.SetProcessInfo( data.CurrentProcessInfo );

		data.MaterialParameters.SetStrainVector( data.E );
		data.MaterialParameters.SetStressVector( data.S );
		data.MaterialParameters.SetConstitutiveMatrix( data.D );

		data.MaterialParameters.SetShapeFunctionsDerivatives( data.dNxy );

		Flags& options = data.MaterialParameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, data.CalculateLHS);

		data.detF = 1.0;
		data.F  = IdentityMatrix(2,2);
		data.MaterialParameters.SetDeterminantF(data.detF);
		data.MaterialParameters.SetDeformationGradientF(data.F);
	}

	void OptTriangleElement::CalculateBMatrix(CalculationData& data)
	{
		//---------------------------------------------
		// geom data
		array_1d<double, 3>& gpLoc = data.gpLocations[data.gpIndex];
		double loc1 = gpLoc[0];
		double loc2 = gpLoc[1];
		double loc3 = gpLoc[2];

		const GeometryType& geom = GetGeometry();
		const NodeType& p1 = geom[0];
		const NodeType& p2 = geom[1];
		const NodeType& p3 = geom[2];

		const double x12 = p1.X0() - p2.X0();
		const double x23 = p2.X0() - p3.X0();
		const double x31 = p3.X0() - p1.X0();
		const double y12 = p1.Y0() - p2.Y0();
		const double y23 = p2.Y0() - p3.Y0();
		const double y31 = p3.Y0() - p1.Y0();

		//---------------------------------------------
		// membrane basic part L
		// already computed. it is constant over the
		// element

		//---------------------------------------------
		// membrane higher order part Qh
		noalias( data.Q )  = loc1 * data.Q1;
		noalias( data.Q ) += loc2 * data.Q2;
		noalias( data.Q ) += loc3 * data.Q3;
		noalias( data.TeQ ) = prod( data.Te, data.Q );
		//noalias( data.Qh  ) = 1.5*std::sqrt(data.beta0) * prod( data.TeQ, data.TTu );
		noalias( data.Qh  ) = prod( data.TeQ, data.TTu );

		//---------------------------------------------
		// membrane strain displacement matrix
		//noalias(data.B) = data.L + 1.5*std::sqrt(data.beta0)*data.Qh;
		noalias(data.B) = data.L + std::sqrt(3.0/4.0*data.beta0)*data.Qh;
	}

	void OptTriangleElement::CalculateBeta0(CalculationData& data)
	{
		data.beta0 = 1.0; // to be changed!
	}

	void OptTriangleElement::CalculateConstitutiveLawResponse(CalculationData& data)
	{
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
		const Matrix & shapeFunctions = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
		for(int nodeid = 0; nodeid < OPT_NUM_NODES; nodeid++)
			data.N(nodeid) = shapeFunctions(data.gpIndex, nodeid);
#else
		const array_1d<double,3>& loc = data.gpLocations[data.gpIndex];
		data.N(0) = 1.0 - loc[1] - loc[2];
		data.N(1) = loc[1];
		data.N(2) = loc[2];
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

		ConstitutiveLaw::Pointer& claw = mConstitutiveLawVector[data.gpIndex];
		data.MaterialParameters.SetShapeFunctionsValues( data.N );
		claw->CalculateMaterialResponseCauchy( data.MaterialParameters );
	}

	void OptTriangleElement::CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS)
	{
		// calculate beta0
		CalculateBeta0( data );

		// calculate the total strain displ. matrix
		CalculateBMatrix( data );

		// compute strain vector
		noalias( data.E ) = prod( data.B, data.U );

		if(data.CalculateRHS || data.CalculateLHS)
		{
			// calculate material response
			CalculateConstitutiveLawResponse( data );

			if(data.CalculateLHS)
			{
				// Add all contributions to the Stiffness Matrix
				noalias( data.BTD ) = prod( trans( data.B ), data.D * data.dV );
				noalias( LHS ) += prod( data.BTD, data.B );
			}

			// Add all contributions to the residual vector
			noalias( RHS ) -= prod( trans( data.B ), data.S * data.dV );
		}
	}

	void OptTriangleElement::AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector)
	{
		const GeometryType& geom = GetGeometry();

		// Get shape functions
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
		const Matrix & N = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
#else
		Matrix N(3,3);
		for(unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
		{
			const array_1d<double,3>& loc = data.gpLocations[igauss];
			N(igauss,0) = 1.0 - loc[1] - loc[2];
			N(igauss,1) = loc[1];
			N(igauss,2) = loc[2];
		}
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

		double rho = GetProperties()[DENSITY];

		// auxiliary
		array_1d<double, 3> bf;

		// gauss loop to integrate the external force vector
		for(unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
		{
			// interpolate nodal volume accelerations to this gauss point
			// and obtain the body force vector
			bf.clear();
			for(unsigned int inode = 0; inode < 3; inode++)
			{
				if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
					bf += N(igauss,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
			}
			bf *= (rho * data.dV);

			// add it to the RHS vector
			for(unsigned int inode = 0; inode < 3; inode++)
			{
				unsigned int index = inode*3;
				double iN = N(igauss,inode);
				rRightHandSideVector[index + 0] += iN * bf[0];
				rRightHandSideVector[index + 1] += iN * bf[1];
				rRightHandSideVector[index + 2]  = 0.0;
			}
		}
	}

    void OptTriangleElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo,
                                            const bool LHSrequired,
                                            const bool RHSrequired)
    {
        // Resize the Left Hand Side if necessary,
        // and initialize it to Zero

        if((rLeftHandSideMatrix.size1() != OPT_NUM_DOFS) || (rLeftHandSideMatrix.size2() != OPT_NUM_DOFS))
            rLeftHandSideMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

        // Resize the Right Hand Side if necessary,
        // and initialize it to Zero

        if(rRightHandSideVector.size() != OPT_NUM_DOFS)
            rRightHandSideVector.resize(OPT_NUM_DOFS, false);
        noalias(rRightHandSideVector) = ZeroVector(OPT_NUM_DOFS);

        // Initialize common calculation variables

        CalculationData data(rCurrentProcessInfo);
		data.CalculateLHS = LHSrequired;
		data.CalculateRHS = RHSrequired;
		InitializeCalculationData(data);

        // Gauss Loop.

		for(size_t i = 0; i < OPT_NUM_GP; i++)
		{
			data.gpIndex = i;
			CalculateGaussPointContribution(data, rLeftHandSideMatrix, rRightHandSideVector);
		}

        // Add body forces contributions. This doesn't depend on the coordinate system

        if(data.CalculateRHS)
			AddBodyForces(data, rRightHandSideVector);
    }

    // =====================================================================================
    //
    // Class OptTriangleElement - Serialization
    //
    // =====================================================================================

    void OptTriangleElement::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
        rSerializer.save("IntM", (int)mThisIntegrationMethod);
    }

    void OptTriangleElement::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
		int temp;
        rSerializer.load("IntM", temp);
		mThisIntegrationMethod = (IntegrationMethod)temp;
    }

}

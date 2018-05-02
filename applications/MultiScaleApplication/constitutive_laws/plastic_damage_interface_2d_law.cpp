/*
==============================================================================
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

/* *********************************************************
 *
 *   Last Modified by:    $Author:   Massimo Petracca$
 *   Date:                $Date:     30-10-2013$
 *   Revision:            $Revision: 1.0$
 *
 * ***********************************************************/

#include "plastic_damage_interface_2d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

//#define PDINTERF_2D_TEST_ELASTIC

namespace Kratos
{


	namespace PDIUtils
	{

		// B3 Function

		inline double b3_eval_area(double x1,double x2,double x3,
								   double y1,double y2,double y3)
		{
			return (x2*y1/3.0 + x3*y1/6.0 - x2*y3/3.0 + x3*y2/3.0 + x3*y3/2.0
					-x1*(y1/2.0 + y2/3.0 + y3/6.0));
		}
		inline double b3_eval_bezier(double xi,double x1,double x2,double x3,
									 double y1,double y2,double y3)
		{
			double A = x1-2.0*x2+x3;
			double B = 2.0*(x2-x1);
			double C = x1-xi;
			if(std::abs(A)<1.0E-12) {
				x2 = x2+1.0E-6*(x3-x1);
				A = x1-2.0*x2+x3;
				B = 2.0*(x2-x1);
				C = x1-xi;
			}
			double D = B*B-4.0*A*C;
			double t = (-B+std::sqrt(D))/(2.0*A);
			return ((y1-2.0*y2+y3)*t*t+2.0*(y2-y1)*t+y1);
		}
		inline void b3_eval_G(double sp,double sj,double sk,double sr,double su,
							  double ep,double ej,double ek,double er,double eu,
							  double& G1, double& G2, double& G3)
		{
			G1 = sp*ep/2.0;
			G2 = b3_eval_area(ep,ej,ek,sp,sj,sk);
			G3 = b3_eval_area(ek,er,eu,sk,sr,su);
		}
		inline double b3_hardening_function(double x,double e0,double ei,double ep,
											double ej,double ek,double er,double eu,
											double s0,double si,double sp,
											double sj,double sk,double sr,double su)
		{
			if(x <= ep)
				return b3_eval_bezier(x,e0,ei,ep,s0,si,sp);
			else if(x <= ek)
				return b3_eval_bezier(x,ep,ej,ek,sp,sj,sk);
			else if(x < eu)
				return b3_eval_bezier(x,ek,er,eu,sk,sr,su);
			else
				return su;
		}
		inline void b3_stretch(double S,double ep,
							   double& ej,double& ek,double& er,double& eu)
		{
			ej = ej+(ej-ep)*S;
			ek = ek+(ek-ep)*S;
			er = er+(er-ep)*S;
			eu = eu+(eu-ep)*S;
		}

	}


    PlasticDamageInterface2DLaw::PlasticDamageInterface2DLaw()
        : ConstitutiveLaw()
		, m_initialized(false)
		, m_error_code(0.0)
#ifdef PDI_USE_SUBINCR
		, m_strain_n()
		, m_strain_n_converged()
#endif // PDI_USE_SUBINCR
	{
    }

    ConstitutiveLaw::Pointer PlasticDamageInterface2DLaw::Clone() const
    {
        return ConstitutiveLaw::Pointer( new PlasticDamageInterface2DLaw() );
    }

    PlasticDamageInterface2DLaw::SizeType PlasticDamageInterface2DLaw::WorkingSpaceDimension()
    {
        return 2;
    }

    PlasticDamageInterface2DLaw::SizeType PlasticDamageInterface2DLaw::GetStrainSize()
    {
        return 2;
    }

    bool PlasticDamageInterface2DLaw::Has(const Variable<double>& rThisVariable)
    {
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE) return true;
		if(rThisVariable == DAMAGE_T) return true;
		if(rThisVariable == DAMAGE_C) return true;
		if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_I) return true;
		if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_II) return true;
		if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_III) return true;
        return false;
    }

    bool PlasticDamageInterface2DLaw::Has(const Variable<Vector>& rThisVariable)
    {
		if(rThisVariable == INTERFACE_PLASTIC_DISPLACEMENT_JUMP) return true;
        return false;
    }

    bool PlasticDamageInterface2DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool PlasticDamageInterface2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }

    bool PlasticDamageInterface2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& PlasticDamageInterface2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
		rValue = 0.0;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		else if(rThisVariable == DAMAGE_T)
			rValue = m_damage_T;
		else if(rThisVariable == DAMAGE_C)
			rValue = m_damage_C;
		else if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_I)
			rValue = (m_lambda(0) + m_alpha*m_lambda(1));
		else if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_II)
			rValue = (m_lambda(0)/m_alpha + m_lambda(1));
		else if(rThisVariable == EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_III)
			rValue = m_lambda(2);
		return rValue;
    }

    Vector& PlasticDamageInterface2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
		if(rThisVariable == INTERFACE_PLASTIC_DISPLACEMENT_JUMP) {
			if(rValue.size() != 3)
				rValue.resize(3,false);
			rValue(0) = m_up(1);
			rValue(1) = m_up(0);
			rValue(2) = 0.0;
		}
        return rValue;
    }

    Matrix& PlasticDamageInterface2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & PlasticDamageInterface2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable,
																 array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & PlasticDamageInterface2DLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable,
																 array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

    void PlasticDamageInterface2DLaw::SetValue(const Variable<double>& rVariable,
                                              const double& rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
		if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			m_lch_multiplier = rValue;
    }

    void PlasticDamageInterface2DLaw::SetValue(const Variable<Vector >& rVariable,
                                              const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::SetValue(const Variable<Matrix >& rVariable,
                                              const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                              const array_1d<double, 3 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                              const array_1d<double, 6 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool PlasticDamageInterface2DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) ) return false;
		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) ) return false;
#ifndef PDINTERF_2D_TEST_ELASTIC
		if( !rMaterialProperties.Has(INTERFACE_TENSILE_LAW_S0) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) ) return false;
		if( !rMaterialProperties.Has(INITIAL_COHESION) ) return false;
		if( !rMaterialProperties.Has(INITIAL_FRICTION_ANGLE) ) return false;
		if( !rMaterialProperties.Has(INITIAL_DILATANCY_ANGLE) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) ) return false;
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_S0) ) return false;
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_SP) ) return false;
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_EP) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_III) ) return false;
#endif // !PDINTERF_2D_TEST_ELASTIC
		return true;
    }

    PlasticDamageInterface2DLaw::StrainMeasure PlasticDamageInterface2DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }

    PlasticDamageInterface2DLaw::StressMeasure PlasticDamageInterface2DLaw::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool PlasticDamageInterface2DLaw::IsIncremental()
    {
        return true;
    }

    void PlasticDamageInterface2DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const Vector& rShapeFunctionsValues)
    {
		if(!m_initialized)
		{
			m_error_code = 0.0;
			m_lch_multiplier = 1.0;
			m_lambda.clear();
			m_up.clear();
			m_up_aux.clear();
			m_lambda_converged.clear();
			m_up_converged.clear();
			m_up_aux_converged.clear();
			m_damage_T = 0.0;
			m_damage_C = 0.0;
			m_alpha = 1.0;
			m_initialized = true;
#ifdef PDI_USE_SUBINCR
			m_strain_n = ZeroVector(this->GetStrainSize());
			m_strain_n_converged = ZeroVector(this->GetStrainSize());
#endif // PDI_USE_SUBINCR
		}
    }

    void PlasticDamageInterface2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
		// save converged values
		m_lambda_converged = m_lambda;
		m_up_converged = m_up;
		m_up_aux_converged = m_up_aux;

#ifdef PDI_USE_SUBINCR
		m_strain_n_converged = m_strain_n;
#endif // PDI_USE_SUBINCR

    }

    void PlasticDamageInterface2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                                  const GeometryType& rElementGeometry,
                                                                  const Vector& rShapeFunctionsValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                                const GeometryType& rElementGeometry,
                                                                const Vector& rShapeFunctionsValues,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void PlasticDamageInterface2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
#ifdef PDI_USE_SUBINCR

		array_1d<double,3> m_lambda_converged__save = m_lambda_converged;
		array_1d<double,2> m_up_converged__save     = m_up_converged;
		array_1d<double,2> m_up_aux_converged__save = m_up_aux_converged;
		Vector m_strain_n__save                     = m_strain_n;
		Vector m_strain_n_converged__save           = m_strain_n_converged;
		double DT_save                              = m_damage_T;
		double DC_save                              = m_damage_C;
		double alpha_save                           = m_alpha;

		unsigned int num_sub_incr = 20;

		Vector dE = rValues.GetStrainVector() - m_strain_n_converged;
		Vector dEsub = dE/double(num_sub_incr);

		noalias(rValues.GetStrainVector()) = m_strain_n_converged;
		for(unsigned int isub = 0; isub < num_sub_incr; isub++)
		{
			rValues.GetStrainVector() += dEsub;
			CalculateAll(rValues);
			if(m_error_code != 0.0)
			{
				m_lambda_converged   = m_lambda_converged__save;
				m_up_converged       = m_up_converged__save;
				m_up_aux_converged   = m_up_aux_converged__save;
				m_strain_n_converged = m_strain_n_converged__save;
				m_damage_T           = DT_save;
				m_damage_C           = DC_save;
				m_alpha              = alpha_save;
				return;
			}
			FinalizeSolutionStep(rValues.GetMaterialProperties(),rValues.GetElementGeometry(),
				rValues.GetShapeFunctionsValues(),rValues.GetProcessInfo());
		}


		if(m_error_code != 0.0) return;
		double hh = 1.0e-9;

		array_1d<double, 3> copy_lambda	= m_lambda;
		array_1d<double, 2> copy_up		= m_up;
		array_1d<double, 2> copy_up_aux	= m_up_aux;
		double copy_alpha	 = m_alpha;
		double copy_damage_T = m_damage_T;
		double copy_damage_C = m_damage_C;

		Vector E = rValues.GetStrainVector();
		Vector S = rValues.GetStressVector();
		Matrix& C = rValues.GetConstitutiveMatrix();

		for(size_t i = 0; i < 2; i++)
		{
			m_lambda_converged   = m_lambda_converged__save;
			m_up_converged       = m_up_converged__save;
			m_up_aux_converged   = m_up_aux_converged__save;
			m_strain_n_converged = m_strain_n_converged__save;
			m_damage_T           = DT_save;
			m_damage_C           = DC_save;
			m_alpha              = alpha_save;

			hh = std::max(1.0e-5, 1.0e-4*std::abs(E(i)));
			noalias(rValues.GetStrainVector()) = E;
			rValues.GetStrainVector()(i) += hh;

			noalias(dE) = rValues.GetStrainVector() - m_strain_n_converged;
			noalias(dEsub) = dE/double(num_sub_incr);

			noalias(rValues.GetStrainVector()) = m_strain_n_converged;
			for(unsigned int isub = 0; isub < num_sub_incr; isub++)
			{
				rValues.GetStrainVector() += dEsub;
				CalculateAll(rValues);
				if(m_error_code != 0.0)
				{
					m_lambda_converged   = m_lambda_converged__save;
					m_up_converged       = m_up_converged__save;
					m_up_aux_converged   = m_up_aux_converged__save;
					m_strain_n_converged = m_strain_n_converged__save;
					m_damage_T           = DT_save;
					m_damage_C           = DC_save;
					m_alpha              = alpha_save;
					return;
				}
				FinalizeSolutionStep(rValues.GetMaterialProperties(),rValues.GetElementGeometry(),
					rValues.GetShapeFunctionsValues(),rValues.GetProcessInfo());
			}

			for(size_t j = 0; j < 2; j++)
			{
				C(j,i) = (rValues.GetStressVector()(j)-S(j))/(hh);
			}
		}

		m_lambda	  = copy_lambda	 ;
		m_up		  = copy_up		 ;
		m_up_aux	  = copy_up_aux	 ;
		m_alpha	      = copy_alpha	 ;
		m_damage_T    = copy_damage_T;
		m_damage_C    = copy_damage_C;

		noalias(rValues.GetStressVector()) = S;
		noalias(rValues.GetStrainVector()) = E;
#else

		CalculateAll(rValues);
		if(m_error_code != 0.0) return;

		double hh = 1.0e-9;

		Vector E = rValues.GetStrainVector();
		Vector S = rValues.GetStressVector();

		array_1d<double, 3> copy_lambda	= m_lambda;
		array_1d<double, 2> copy_up		= m_up;
		array_1d<double, 2> copy_up_aux	= m_up_aux;
		double copy_alpha	 = m_alpha;
		double copy_damage_T = m_damage_T;
		double copy_damage_C = m_damage_C;

		Vector S2(2);
		Vector S1(2);

		Matrix& C = rValues.GetConstitutiveMatrix();
		for(size_t i = 0; i < 2; i++)
		{
			hh = std::max(1.0e-8, 1.0e-5*std::abs(E(i)));
			noalias(rValues.GetStrainVector()) = E;
			rValues.GetStrainVector()(i) -= hh;
			CalculateAll(rValues);
			noalias(S1) = rValues.GetStressVector();
			rValues.GetStrainVector()(i) += 2.0*hh;
			CalculateAll(rValues);
			noalias(S2) = rValues.GetStressVector();
			for(size_t j = 0; j < 2; j++)
			{
				C(j,i) = (S2(j)-S1(j))/(2.0*hh);
			}
		}

		noalias(rValues.GetStressVector()) = S;
		noalias(rValues.GetStrainVector()) = E;

		m_lambda	  = copy_lambda	 ;
		m_up		  = copy_up		 ;
		m_up_aux	  = copy_up_aux	 ;
		m_alpha	      = copy_alpha	 ;
		m_damage_T    = copy_damage_T;
		m_damage_C    = copy_damage_C;

#endif // PDI_USE_SUBINCR
    }

    void PlasticDamageInterface2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void PlasticDamageInterface2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {
    }

    void PlasticDamageInterface2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector& rShapeFunctionsValues)
    {
		m_initialized = false;
		m_lch_multiplier = 1.0;
		m_error_code = 0.0;

		m_lambda.clear();
		m_up.clear();
		m_up_aux.clear();

		m_lambda_converged.clear();
		m_up_converged.clear();
		m_up_aux_converged.clear();
    }

    void PlasticDamageInterface2DLaw::GetLawFeatures(Features& rFeatures)
    {
        //Set the type of law
		rFeatures.mOptions.Set( PLANE_STRESS_LAW ); // TODO: INTERFACE 2D LAW
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }

    int PlasticDamageInterface2DLaw::Check(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: NORMAL_STIFFNESS", "");
		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TANGENTIAL_STIFFNESS", "");
#ifndef PDINTERF_2D_TEST_ELASTIC
		if( !rMaterialProperties.Has(INTERFACE_TENSILE_LAW_S0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERFACE_TENSILE_LAW_S0", "");
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_I", "");
		if( !rMaterialProperties.Has(INITIAL_COHESION) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_COHESION", "");
		if( !rMaterialProperties.Has(INITIAL_FRICTION_ANGLE) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_FRICTION_ANGLE", "");
		if( !rMaterialProperties.Has(INITIAL_DILATANCY_ANGLE) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_DILATANCY_ANGLE", "");
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_II", "");
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_S0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERFACE_COMPRESSIVE_LAW_S0", "");
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_SP) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERFACE_COMPRESSIVE_LAW_SP", "");
		if( !rMaterialProperties.Has(INTERFACE_COMPRESSIVE_LAW_EP) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERFACE_COMPRESSIVE_LAW_EP", "");
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_III) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_III", "");
#endif // !PDINTERF_2D_TEST_ELASTIC

		return 0;

		KRATOS_CATCH("");
	}

	void PlasticDamageInterface2DLaw::InitializeData(Parameters& rValues, CalculationData& d)
	{
		// properties
		const Properties& props = rValues.GetMaterialProperties();
		const ProcessInfo& pinfo = rValues.GetProcessInfo();

		// elasticity
		d.kn = props[NORMAL_STIFFNESS];
		d.kt = props[TANGENTIAL_STIFFNESS];
		// tensile data (mode I)
		d.ft  = props[INTERFACE_TENSILE_LAW_S0];
		d.GfI = props[FRACTURE_ENERGY_MODE_I];
		// shear data (mode II)
		d.c = props[INITIAL_COHESION];
		d.tanphi0 = std::tan( Globals::Pi/180.0 * props[INITIAL_FRICTION_ANGLE] );
		d.tanphiu = d.tanphi0;
		if(props.Has(RESIDUAL_FRICTION_ANGLE))
			d.tanphiu = std::tan( Globals::Pi/180.0 * props[RESIDUAL_FRICTION_ANGLE] );
		d.tanpsi0 = std::tan( Globals::Pi/180.0 * props[INITIAL_DILATANCY_ANGLE] );
		d.tanpsiu = d.tanpsi0;
		if(props.Has(RESIDUAL_DILATANCY_ANGLE))
			d.tanpsiu = std::tan( Globals::Pi/180.0 * props[RESIDUAL_DILATANCY_ANGLE] );
		d.GfII = props[FRACTURE_ENERGY_MODE_II];
		// cap data (mode III)
		d.fc0 = props[INTERFACE_COMPRESSIVE_LAW_S0];
		d.fcp = props[INTERFACE_COMPRESSIVE_LAW_SP];
		d.fcr = 0.0;
		if(props.Has(INTERFACE_COMPRESSIVE_LAW_SR)) d.fcr = props[INTERFACE_COMPRESSIVE_LAW_SR];
		d.ecp = props[INTERFACE_COMPRESSIVE_LAW_EP];
		d.Css = 9.0;
		if(props.Has(INTERFACE_CAP_VALUE)) d.Css = props[INTERFACE_CAP_VALUE];
		d.GfIII = props[FRACTURE_ENERGY_MODE_III];
		d.c1 = 0.6;
		d.c2 = 0.5;
		d.c3 = 1.5;
		d.c4 = 0.0;
		if(props.Has(INTERFACE_COMPRESSIVE_LAW_C1)) d.c1 = props[INTERFACE_COMPRESSIVE_LAW_C1];
		if(props.Has(INTERFACE_COMPRESSIVE_LAW_C2)) d.c2 = props[INTERFACE_COMPRESSIVE_LAW_C2];
		if(props.Has(INTERFACE_COMPRESSIVE_LAW_C3)) d.c3 = props[INTERFACE_COMPRESSIVE_LAW_C3];
		if(props.Has(INTERFACE_COMPRESSIVE_LAW_C4)) d.c4 = props[INTERFACE_COMPRESSIVE_LAW_C4];
		// regularization (interfaces don't need it! but it can be used with the lch multiplier...)
		d.lch = 1.0 * m_lch_multiplier;
		d.GfI /= d.lch;
		d.GfII /= d.lch;
		d.GfIII /= d.lch;
		// viscosity
		d.eta_over_dt = 0.0;
		if(props.Has(VISCOSITY)) {
			double eta = props[VISCOSITY];
			if(eta > 0.0) {
				double dt = pinfo[DELTA_TIME];
				if(dt > 0.0) {
					d.eta_over_dt = eta/dt;
				}
			}
		}
		// plastic damage factors
		d.Hpt = 0.0;
		if(props.Has(INTERFACE_PLASTIC_DAMAGE_FACTOR_T))
			d.Hpt = std::max(0.0, std::min(1.0, props[INTERFACE_PLASTIC_DAMAGE_FACTOR_T]));
		d.Hpc = 0.0;
		if(props.Has(INTERFACE_PLASTIC_DAMAGE_FACTOR_C))
			d.Hpc = std::max(0.0, std::min(1.0, props[INTERFACE_PLASTIC_DAMAGE_FACTOR_C]));
		// misc
		d.alpha = d.GfI / d.GfII * d.c / d.ft;
		m_alpha = d.alpha;
		d.use_damage = true; // TODO: optional in properties!!
	}

	void PlasticDamageInterface2DLaw::HardeningModeI(CalculationData& d, double lambdaI, double& q, double& dq)
	{
		double A = std::exp(-d.ft/d.GfI*lambdaI);
		q = d.ft*A;
		dq = -d.ft*d.ft/d.GfI*A;
	}
	void PlasticDamageInterface2DLaw::HardeningModeII(CalculationData& d, double lambdaII, double& q, double& dq, double& tanphi)
	{
		double A = std::exp(-d.c/d.GfII*lambdaII);
		q = d.c*A;
		tanphi = d.tanphi0+(d.tanphiu-d.tanphi0)*(d.c-q)/d.c;
		dq = -d.c*d.c/d.GfII*A;
	}
	void PlasticDamageInterface2DLaw::HardeningModeIII(CalculationData& d, double lambdaIII, double& q, double& dq)
	{
		double kn    = d.kn;
		double s0    = d.fc0;
		double sp    = d.fcp;
		double sr    = d.fcr;
		double ep    = d.ecp;
		double GfIII = d.GfIII;
		double c1    = d.c1;
		double c2    = d.c2;
		double c3    = d.c3;
		double c4    = d.c4;
		// fracture energy
		double lch   = d.lch;
		double Gdis = GfIII; // /lch; // note: already regularized in initialize data
		// auto-computed stress parameters
		double si = sp;
		double sj = sp;
		double sk = sr + (sp-sr)*c1;
		double su = sr;
		// auto-computed strain parameters
		double e0 = s0/kn;
		double ei = si/kn;
		ei = ei+(ep-ei)*c4;
		double gamma = 2.0*(ep-sp/kn);
		double ej = ep + gamma*c2;
		double ek = ej + gamma*(1.0-c2);
		double er = (ek-ej)/(sp-sk)*(sp-sr)+ej;
		double eu = er*c3;
		// this is the initial area under the curve
		// (i.e. the area given by the curve inputed by the user)
		double G1,G2,G3;
		PDIUtils::b3_eval_G(sp,sj,sk,sr,su,ep,ej,ek,er,eu,G1,G2,G3);
		// modify G1
		double Hpc     = d.Hpc;
		double spp     = Hpc*kn*ep + (1.0-Hpc)*sp;
		double damage  = 1.0-sp/spp;
		double KD      = (1.0-damage)*kn;
		double eei     = sp/KD;
		G1      = sp*eei/2.0;
		// sum
		double Gtot = G1+G2+G3;
		// this computes the stretch factor required to obtain
		// a curve whose area is equal to G/lch
		double S = (Gdis-G1)/(Gtot-G1)-1.0;
		if(S<=-1.0) {
			std::stringstream ss;
			ss << "WARNING - Interface Hardening Compression. G is too low\n" << std::endl;
			S = -0.99999;
		}
		PDIUtils::b3_stretch(S,ep,ej,ek,er,eu);
		// move to "plastic_strain - stress" space
		e0=e0-s0/kn;
		ei=ei-si/kn;
		ep=ep-sp/kn;
		ej=ej-sj/kn;
		ek=ek-sk/kn;
		er=er-sr/kn;
		eu=eu-su/kn;
		q = PDIUtils::b3_hardening_function(lambdaIII,e0,ei,ep,ej,ek,er,eu,s0,si,sp,sj,sk,sr,su);
		double hh = 1.0e-10;
		double q1 = PDIUtils::b3_hardening_function(lambdaIII+hh,e0,ei,ep,ej,ek,er,eu,s0,si,sp,sj,sk,sr,su);
		dq = (q1-q)/hh;
	}

	double PlasticDamageInterface2DLaw::GetBetaTension(CalculationData& d,double l1,double l2,double up1)
	{
		double alpha = d.alpha;
		double upeq = l1 + alpha * l2;
		double q,dq;
		HardeningModeI(d,upeq,q,dq);
		double kn  = d.kn;
		double ft  = d.ft;
		double et = up1+q/kn;
		double Hpt = d.Hpt;
		double qeff;
		if(Hpt < 0.5) {
			double H = Hpt*2.0;
			qeff = H*ft + (1-H)*q;
		}
		else {
			double H = (Hpt-0.5)*2.0;
			qeff = ft + H*kn*(et-ft/kn);
		}
		double upeff = et-qeff/kn;
		double beta_t = (up1-upeff)/up1;
		if(beta_t < 1.0e-10) beta_t=0.0;
		return beta_t;
	}
	double PlasticDamageInterface2DLaw::GetBetaCompression(CalculationData& d,double l3,double up3)
	{
		// up3 = abs(up3); // <<<<<<<<<<<< WARNING
		// l3 = up3;
		up3 = l3;
		// real hardening variable
		double q,dq;
		HardeningModeIII(d,l3,q,dq);
		// material parameters
		double kn    = d.kn;
		double s0    = d.fc0;
		double sp    = d.fcp;
		double sr    = d.fcr;
		double ep    = d.ecp;
		double GfIII = d.GfIII;
		double c1 = d.c1;
		double c2 = d.c2;
		double c3 = d.c3;
		double c4 = d.c4;
		// fracture energy
		double lch   = d.lch;
		double Gdis = GfIII; // /lch; // note: already regularized in initialize data
		// auto-computed stress parameters
		double si = sp;
		double sj = sp;
		double sk = sr + (sp-sr)*c1;
		double su = sr;
		// auto-computed strain parameters
		double e0 = s0/kn;
		double ei = si/kn;
		ei = ei+(ep-ei)*c4;
		double gamma = 2.0*(ep-sp/kn);
		double ej = ep + gamma*c2;
		double ek = ej + gamma*(1.0-c2);
		double er = (ek-ej)/(sp-sk)*(sp-sr)+ej;
		double eu = er*c3;
		// this is the initial area under the curve
		// (i.e. the area given by the curve inputed by the user)
		double G1,G2,G3;
		PDIUtils::b3_eval_G(sp,sj,sk,sr,su,ep,ej,ek,er,eu,G1,G2,G3);
		// modify G1
		double Hpc     = d.Hpc;
		double spp     = Hpc*kn*ep + (1.0-Hpc)*sp;
		double damage  = 1.0-sp/spp;
		double KD      = (1.0-damage)*kn;
		double eei     = sp/KD;
		G1      = sp*eei/2.0;
		// sum
		double Gtot = G1+G2+G3;
		// this computes the stretch factor required to obtain
		// a curve whose area is equal to G/lch
		double S = (Gdis-G1)/(Gtot-G1)-1.0;
		if(S<=-1.0) {
			std::stringstream ss;
			ss << "WARNING - Interface Hardening Compression. G is too low\n" << std::endl;
			S = -0.99999;
		}
		PDIUtils::b3_stretch(S,ep,ej,ek,er,eu);
		// modify the control points for the effective hardening variable
		sj = Hpc*kn*ej + (1-Hpc)*sj;
		sk = Hpc*kn*ek + (1-Hpc)*sk;
		sr = Hpc*kn*er + (1-Hpc)*sr;
		double HKR = (sr-sk)/(er-ek);
		if(HKR<0.0)
			su = sr;
		else
			su = sr+HKR*(eu-er);
		// total strain
		double et = up3+q/kn;
		// evaluate the effective hardening variable
		double qeff;
		if(Hpc==1.0)
			qeff = kn*et;
		else if(Hpc == 0.0)
			qeff = q;
		else
			qeff = PDIUtils::b3_hardening_function(
				et,e0,si/kn,ep,ej,ek,er,eu,s0,si,spp,sj,sk,sr,su);
		// compute beta
		double upeff = et-qeff/kn;
		double beta_c = (up3-upeff)/up3;
		if(beta_c < 1.0e-10) beta_c=0.0;
		return beta_c;
	}

	void PlasticDamageInterface2DLaw::CalculateDamage(CalculationData& d,double l1,double l2,double l3,double up1,double up3,
													  double& dt, double& dc)
	{
		double q,dq,upeq,et;
		// mode 1
		upeq = l1 + d.alpha * l2;
		HardeningModeI(d,upeq,q,dq);
		et = up1+q/d.kn;
		double kdt = q/et;
		dt = 1.0-kdt/d.kn;
		// mode 3
		upeq = l3;
		HardeningModeIII(d,upeq,q,dq);
		et = std::abs(up3)+q/d.kn;
		double kdc = q/et;
		dc = 1.0-kdc/d.kn;
	}

	bool PlasticDamageInterface2DLaw::RMapModeI(CalculationData& d, const Matrix& tule, const Vector& trat,
												double tolerance, double sign_tau,
												Vector& tra, Matrix& tunl, double& ft,
												double& upeq, array_1d<double,2>& dup, double& dup1, double& dl, double& ovs1)
	{
		// constant parameters
		int miter = 100;

		// initialize
		dl = 0.0;
		double dyiedk = 0.0;
		tra.clear();

		// shear is constant in this return mapping
		tra(1) = trat(1);

		// Newton loop to solve for the plastic multiplier
		bool converged = false;
		double f=1.0;
		for (int iter=1; iter<=miter; iter++)
		{
			HardeningModeI(d, upeq+dl, ft, dyiedk);
			tra(0) = trat(0)-dl*tule(0,0);
			f = tra(0)-ft;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			f = f -d.eta_over_dt*dl;
			ovs1 = d.eta_over_dt*dl;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			if(std::abs(f) < tolerance) {
				converged = true;
				break;
			}

			double dfdl = -tule(0,0)-dyiedk;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			dfdl = dfdl -d.eta_over_dt;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			dl = dl-f/dfdl;
		}

		// check convergence
		if(!converged) {
			std::stringstream ss;
			ss << "Return mapping - Mode I : maximum iterations reached. F = " << f << std::endl;
			std::cout << ss.str();
		}

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		dup(0)=dl; dup(1)=0.0;
		dup1 = dl;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// updated plastic multiplier
		upeq = upeq+dl;

		// correct the sign of tau
		tra(1) = tra(1)*sign_tau;

		// consistent tangent operator
		tunl(0,0) = dyiedk*tule(0,0)/(dyiedk+tule(0,0));
		tunl(1,1) = tule(1,1);
		tunl(0,1) = 0.0;
		tunl(1,0) = 0.0;

		// return
		return converged;
	}

	bool PlasticDamageInterface2DLaw::RMapModeII(CalculationData& d, const Matrix& tule, const Vector& trat,
												 double tolerance, double sign_tau,
												 Vector& tra, Matrix& tunl, double& c, double& tanphi,
												 double& upeq, array_1d<double,2>& dup, double& dl, double& ovs2)
	{
		// constant parameters
		int miter = 100;

		// get material paramters
		double c0      = d.c;
		double tanphi0 = d.tanphi0;
		double tanphiu = d.tanphiu;
		double tanpsi0 = d.tanpsi0;
		double tanpsiu = d.tanpsiu;
		double unconf  = -1.0e20 * c0;

		// initialize
		dl     = 0.0;
		double dyiedk = 0.0;
		tra.clear();

		double tanpsi = 0.0;
		double rdum = 0.0;

		// Newton loop to solve for the plastic multiplier
		bool converged = false;
		double f=1.0;
		for (int iter = 1; iter <= miter; iter++)
		{
			HardeningModeII(d, upeq+dl, c, dyiedk, tanphi);
			tanpsi = tanpsi0 + (tanpsiu-tanpsi0) * (c0-c)/c0;
			rdum   = 1.0/(1.0 - dl*tule(0,0)*tanpsi/unconf);
			tra(0) = rdum*(trat(0)-dl*tule(0,0)*tanpsi);

			// adjust the dilatancy angle to zero for tra(0) higher than unconf
			if(tra(0) <= unconf) {
				rdum = 0.0;
				tra(0) = trat(0);
			}

			tra(1) = trat(1) - dl*tule(1,1);
			f      = tra(1) + tra(0)*tanphi - c;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			f      = f -d.eta_over_dt*dl;
			ovs2   = d.eta_over_dt*dl;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			if(std::abs(f) < tolerance) {
				converged = true;
				break;
			}

			double dfdl = - tule(1,1) + (tra(0)/unconf-1.0) *
							tule(0,0) * (tanpsi-dl*(tanpsiu-tanpsi0)/c0*dyiedk ) *
							rdum * tanphi -
							tra(0) * (tanphiu-tanphi0)/c0*dyiedk - dyiedk;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			dfdl = dfdl -d.eta_over_dt;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			dl = dl - f/dfdl;

		}

		// check convergence
		if(!converged) {
			std::stringstream ss;
			ss << "Return mapping - Mode II : maximum iterations reached. F = " << f << std::endl;
			std::cout << ss.str();
		}

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		dup(0)=dl*tanpsi; dup(1)=dl*sign_tau;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// updated plastic multiplier
		upeq = upeq+dl;

		// correct the sign of tau
		tra(1) = tra(1)*sign_tau;

		// ------------------  consistent tangent operator -------------------------

		rdum  = (1.0+tra(0)*(tanphiu-tanphi0)/c0)*dyiedk;

		// adjust the dilatancy angle to zero for tra(0) higher than unconf
		if(tra(1) <= unconf) {
			tanpsi0 = 0.0;
			tanpsiu = 0.0;
		}

		tanpsi = tanpsi0 + ( tanpsiu - tanpsi0 ) * ( c0 - c ) / c0;
		tanpsi = tanpsi * ( 1.0 - tra(0) / unconf );
		array_1d<double,2> gradf, gradg;
		gradf(0)=tanphi; gradf(1)=sign_tau;
		gradg(0)=tanpsi; gradg(1)=sign_tau;

		gradg(0) = gradg(0) - dl*(tanpsiu-tanpsi0)/c0*dyiedk*(1.0-tra(0)/unconf);

		tunl(0,0) = 1.0/(1.0/tule(0,0)+dl*tanpsi/unconf);
		tunl(1,1) = tule(1,1);
		tunl(0,1) = 0.0;
		tunl(1,0) = 0.0;

		array_1d<double,2> dgradg;
		noalias(dgradg) = prod(tunl,gradg);
		double beta = inner_prod(gradf,dgradg);

		array_1d<double,2> dgradf;
		noalias(dgradf) = prod(trans(tunl),gradf);

		tunl -= 1.0/(beta+rdum)*outer_prod(dgradg,dgradf);

		// return
		return converged;
	}

	bool PlasticDamageInterface2DLaw::RMapModeIII(CalculationData& d, const Matrix& tule, const Vector& trat,
												  double tolerance, double sign_tau,
												  Vector& tra, Matrix& tunl, double& fc,
												  double& upeq, array_1d<double,2>& dup, double& dup3, double& dl, double& ovs3)
	{
		// constant parameters
		int miter = 100;

		// some compressive data
		double Css = d.Css;
		Matrix p(2,2,0.0);
		p(0,0)=2.0;
		p(1,1)=2.0*Css;

		// initialize
		dl     = 0.0;
		double dyiedk = 0.0;
		tra.clear();

		// inverse eleasticity tensor (diagonal)
		Matrix invd(2,2,0.0);
		invd(0,0)=1.0/tule(0,0);
		invd(1,1)=1.0/tule(1,1);

		array_1d<double,2> uela;
		noalias(uela) = prod(invd,trat);

		// Newton loop to solve for the
		// plastic multiplier
		bool converged = false;
		double f=1.0;
		Matrix a(2,2,0.0);
		array_1d<double,2> gradf;
		double norm_gradf(0.0);
		double dupeq(0.0);
		double dfdk(0.0);
		array_1d<double,2> dkdsig;
		array_1d<double,2> dsigdl;
		double dkdl(0.0);
		for (int iter=1; iter<=miter; iter++)
		{
			//a = inv(invd+dl*p); // diagonal
			a(0,0) = 1.0/(invd(0,0)+dl*p(0,0));
			a(1,1) = 1.0/(invd(1,1)+dl*p(1,1));
			noalias(tra) = prod(a,uela);

			// Calculate equivalent plastic strain
			noalias(gradf) = prod(p,tra);
			norm_gradf = norm_2(gradf);
			dupeq = dl*norm_gradf;

			// hardening function
			HardeningModeIII(d, upeq+dupeq, fc, dyiedk);

			// yield function
			f = 0.5 * inner_prod(tra,gradf) - fc*fc;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			f = f -d.eta_over_dt*dupeq;
			ovs3 = d.eta_over_dt*dupeq;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			if(std::abs(f) < tolerance) {
				converged = true;
				break;
			}

			dfdk   = -2.0 * fc * dyiedk;
			noalias(dkdsig) = dl/norm_gradf * prod(p,gradf);
			noalias(dsigdl) = -prod(a,gradf);
			dkdl   = norm_gradf;

			double dfdl = inner_prod(gradf,dsigdl) +
						  dfdk*inner_prod(dkdsig,dsigdl) + dfdk*dkdl;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			dfdl = dfdl -d.eta_over_dt*norm_gradf;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			dl = dl-f/dfdl;

		}

		// check convergence
		if(!converged) {
			std::stringstream ss;
			ss << "Return mapping - Mode III : maximum iterations reached. F = " << f << std::endl;
			std::cout << ss.str();
		}

		// updated plastic multiplier
		upeq = upeq+dupeq;

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		// dup = dl*gradf; // or dupeq*gradf/norm(gradf)
		// Note: gradf3 is the normal to the cap surface defined on the positive
		// part of tau. so the second component of gradf3 should be
		// modified with the sign of tau!
		dup(0) = dl*gradf(0); dup(1)=dl*gradf(1)*sign_tau;
		dup3 = dup(0);
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// correct the sign of tau
		tra(1) = tra(1)*sign_tau;

		// consistent tangent operator
		noalias(gradf) = prod(p,tra);
		dkdl   = norm_2(gradf);
		noalias(dkdsig) = dl/norm_gradf * prod(p,gradf);
		dfdk   = -2.0 * fc * dyiedk;

		// correct [d]
		// tunl = invd + dl*p;
		// tunl = inv(tunl); // note: diagonal
		tunl.clear();
		tunl(0,0) = 1.0/(invd(0,0)+dl*p(0,0));
		tunl(1,1) = 1.0/(invd(1,1)+dl*p(1,1));

		// make gradg
		array_1d<double,2> gradg = gradf;
		// correct gradf
		gradf += dfdk*dkdsig;

		array_1d<double,2> dgradg(prod(tunl,gradg));
		double beta = inner_prod(gradf,dgradg);
		array_1d<double,2> dgradf(prod(trans(tunl),gradf));

		tunl += -1.0/(beta-dfdk*dkdl) * outer_prod(dgradg,dgradf);

		// return
		return converged;
	}

	bool PlasticDamageInterface2DLaw::RMapCorner_I_II(CalculationData& d, const Matrix& tule, const Vector& trat,
													  double tolerance, double sign_tau,
													  Vector& tra, Matrix& tunl, double& ft, double& c, double& tanphi,
													  double& upeq1, double& upeq2, array_1d<double,2>& dup, double& dup1,
													  double& dl1, double& dl2, double& ovs1, double& ovs2)
	{
		// constant parameters
		int miter = 100;

		// get material paramters
		double c0      = d.c;
		double tanphi0 = d.tanphi0;
		double tanphiu = d.tanphiu;
		double tanpsi0 = d.tanpsi0;
		double tanpsiu = d.tanpsiu;
		double alpha   = d.alpha;

		// initialize
		dl1 = 0.0;
		dl2 = 0.0;
		array_1d<double,2> x;
		x(0)=0.0;
		x(1)=0.0;
		tra.clear();
		Matrix jacob(2,2,0.0);
		Matrix invmat(2,2,0.0);
		double dummy_det=0.0;
		double dyi1dk=0.0;
		double dyi2dk=0.0;
		double dk1=0.0;
		double dk2=0.0;

		double dk1dl1, dk1dl2, dk2dl1, dk2dl2;
		double tanpsi;

		// Newton loop to solve for the plastic multiplier
		bool converged = false;
		double f = 1.0;
		array_1d<double,2> func;
		for (int iter = 1; iter <=miter; iter++)
		{
			dl1 = x(0);
			dl2 = x(1);

			dk1 = std::sqrt( dl1 * dl1 + dl2 * dl2 * alpha * alpha );
			dk2 = std::sqrt( dl1 * dl1 / alpha / alpha + dl2 * dl2 );

			double rdum1 = upeq1 + alpha * upeq2 + dk1;
			double rdum2 = upeq1 / alpha + upeq2 + dk2;
			HardeningModeI(d, rdum1, ft, dyi1dk);
			HardeningModeII(d, rdum2, c, dyi2dk, tanphi);
			tanpsi = tanpsi0 + ( tanpsiu - tanpsi0 ) * ( c0 - c ) / c0;

			tra(0) = trat(0) - dl1 * tule(0,0) - dl2 * tule(0,0) * tanpsi;
			tra(1) = trat(1) - dl2 * tule(1,1);
			double f1 = tra(0) - ft;
			double f2 = tra(1) + tra(0) * tanphi - c;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			f1 = f1 -d.eta_over_dt*dl1;
			f2 = f2 -d.eta_over_dt*dl2;
			ovs1 = d.eta_over_dt*dl1;
			ovs2 = d.eta_over_dt*dl2;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			f = std::sqrt(f1*f1+f2*f2);
			if(std::abs(f) < tolerance) {
				converged=true;
				break;
			}

			func(0)=f1;
			func(1)=f2;

			if (dk1 > 1.0e-14) {
				dk1dl1 = dl1 / dk1;
				dk1dl2 = alpha * alpha * dl2 / dk1;
			}
			else {
				dk1dl1 = 1.0;
				dk1dl2 = 0.0;
			}
			if(dk2 > 1.0e-14) {
				dk2dl1 = dl1 / alpha / alpha / dk2;
				dk2dl2 = dl2 / dk2;
			}
			else {
				dk2dl1 = 0.0;
				dk2dl2 = 1.0;
			}

			double psidk2 = - ( tanpsiu - tanpsi0 ) / c0 * dyi2dk;
			double phidk2 = - ( tanphiu - tanphi0 ) / c0 * dyi2dk ;
			double sigdl1 = - tule(0,0) - dl2 * tule(0,0) * psidk2 * dk2dl1;
			double sigdl2 = - tule(0,0) * tanpsi - dl2 * tule(0,0) * psidk2 * dk2dl2;
			double taudl2 = - tule(1,1);

			jacob(0,0) = sigdl1 - dyi1dk * dk1dl1;
			jacob(0,1) = sigdl2 - dyi1dk * dk1dl2;
			jacob(1,0) =          sigdl1 * tanphi + tra(0) * phidk2 * dk2dl1 - dyi2dk * dk2dl1;
			jacob(1,1) = taudl2 + sigdl2 * tanphi + tra(0) * phidk2 * dk2dl2 - dyi2dk * dk2dl2;

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			jacob(0,0) -= d.eta_over_dt;
			jacob(1,1) -= d.eta_over_dt;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			MathUtils<double>::InvertMatrix2(jacob,invmat,dummy_det);
			x -= prod(invmat,func);

		}

		// check convergence
		if(!converged) {
			std::stringstream ss;
			ss << "Return mapping - Corner I & II : maximum iterations reached. F = " << f << std::endl;
			std::cout << ss.str();
		}

		// Do not add the entire dk1 and dk2. A division directly proportional
		// to dl1  and dl2 is assumed
		upeq1 = upeq1 + dk1 * dl1 / ( dl1 + dl2 );
		upeq2 = upeq2 + dk2 * dl2 / ( dl1 + dl2 );
		tra(1) = tra(1) * sign_tau;

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		dup(0)=dl1 + dl2*tanpsi;
		dup(1)=dl2*sign_tau;
		dup1 = dl1;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// Consistent Tangent
		double rdum = ( 1.0 + tra(0) * ( tanphiu - tanphi0 ) / c0 ) * dyi2dk;
		array_1d<double,2> gradf1, gradf2, gradg1, gradg2;
		gradf1(0)=1.0; gradf1(1)=0.0;
		gradg1(0)=1.0; gradg1(1)=0.0;
		gradf2(0)=tanphi; gradf2(1)=sign_tau;
		gradg2(0)=tanpsi; gradg2(1)=sign_tau;

		if (dk1 > 1.0e-14) {
			dk1dl1 = dl1 / dk1;
			dk1dl2 = alpha * alpha * dl2 / dk1;
		}
		else {
			dk1dl1 = 1.0;
			dk1dl2 = 0.0;
		}
		if(dk2 > 1.0e-14) {
			dk2dl1 = dl1 / alpha / alpha / dk2;
			dk2dl2 = dl2 / dk2;
		}
		else {
			dk2dl1 = 0.0;
			dk2dl2 = 1.0;
		}

		// correct gradg
		gradg1(0) = gradg1(0) - dl2 * ( tanpsiu-tanpsi0 ) / c0 * dyi2dk * dk2dl1;
		gradg2(0) = gradg2(0) - dl2 * ( tanpsiu-tanpsi0 ) / c0 * dyi2dk * dk2dl2;

		Matrix u(2,2);
		u(0,0)=gradg1(0); u(0,1)=gradg2(0);
		u(1,0)=gradg1(1); u(1,1)=gradg2(1);
		Matrix v(2,2);
		v(0,0)=gradf1(0); v(0,1)=gradf2(0);
		v(1,0)=gradf1(1); v(1,1)=gradf2(1);
		Matrix e(2,2);
		e(0,0) = - dyi1dk * dk1dl1;
		e(0,1) = - dyi1dk * dk1dl2;
		e(1,0) = - rdum * dk2dl1;
		e(1,1) = - rdum * dk2dl2;

		Matrix vtd(prod(trans(v),tule));
		Matrix mdum(prod(vtd,u));
		mdum-=e;
		MathUtils<double>::InvertMatrix2(mdum,invmat,dummy_det);
		noalias(mdum) = prod(invmat,vtd);
		Matrix du (prod(tule,u));
		noalias(tunl)= tule - prod(du,mdum);

		// return
		return converged;
	}

	bool PlasticDamageInterface2DLaw::RMapCorner_II_III(CalculationData& d, const Matrix& tule, const Vector& trat,
														double tolerance, double sign_tau,
														Vector& tra, Matrix& tunl, double& fc, double& c, double& tanphi,
														double& upeq2, double& upeq3, array_1d<double,2>& dup, double& dup3,
														double& dl2, double& dl3, double& ovs2, double& ovs3)
	{
		// constant parameters
		int miter = 1000;

		// get material paramters
		double c0      = d.c;
		double tanphi0 = d.tanphi0;
		double tanphiu = d.tanphiu;
		double tanpsi0 = d.tanpsi0;
		double tanpsiu = d.tanpsiu;
		double unconf  = -1.0e20 * c0;
		double Css = d.Css;
		Matrix p(2,2,0.0);
		p(0,0)=2.0;
		p(1,1)=2.0*Css;

		// initialize
		dl2    = 0.0;
		dl3    = 0.0;
		array_1d<double,2> x;
		x(0)=0.0;
		x(1)=0.0;
		tra.clear();
		Matrix jacob(2,2,0.0);
		Matrix invmat(2,2,0.0);
		double dummy_det=0.0;
		double dyi2dk = 0.0;
		double dyi3dk = 0.0;

		double rdum2, rdum3;
		double tanpsi;
		array_1d<double,2> gradf3;
		array_1d<double,2> func;
		double norm_gradf3;
		double dk3dl3;

		// Newton loop to solve for the plastic multiplier
		bool converged = false;
		double f = 1.0;
		for (int iter = 1; iter <= miter; iter++)
		{
			dl2 = x(0);
			dl3 = x(1);
			dl3 = std::max(dl3,0.0);

			// equivalent plastic strain for shear
			rdum2 = upeq2 + dl2;
			HardeningModeII(d,rdum2,c,dyi2dk,tanphi);
			tanpsi = tanpsi0 + (tanpsiu-tanpsi0) * (c0-c)/c0;
			double rdum   = 1.0/(1.0 - dl2*tule(0,0)*tanpsi/unconf + 2.0*dl3*tule(0,0));
			double rdum1  = 1.0/(1.0 + 2.0*Css*dl3*tule(1,1));
			tra(0) = rdum * (trat(0)-dl2*tule(0,0)*tanpsi);

			// adjust the dilatancy angle to zero for tra(1) higher than unconf
			if (tra(0) <= unconf) {
				tanpsi0 = 0.0;
				tanpsiu = 0.0;
				tanpsi = 0.0;
				rdum   = 1.0/(2.0*dl3*tule(0,0));
				tra(0) = rdum * trat(0);
			}
			tra(1) = rdum1 * (trat(1)-dl2*tule(1,1));

			// calculate equivalent plastic strain for cap
			noalias(gradf3) = prod(p,tra);
			norm_gradf3 = norm_2(gradf3);
			rdum3 = upeq3 + dl3*norm_gradf3;
			HardeningModeIII(d,rdum3,fc,dyi3dk);

			double f2 = tra(1) + tra(0) * tanphi - c;
			double f3 = 0.5 * inner_prod(tra,gradf3) - fc*fc;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			f2 = f2 -d.eta_over_dt*dl2;
			f3 = f3 -d.eta_over_dt*dl3*norm_gradf3;
			ovs2 = d.eta_over_dt*dl2;
			ovs3 = d.eta_over_dt*dl3*norm_gradf3;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			f = std::sqrt(f2*f2+f3*f3);
			if(std::abs(f) < tolerance) {
				converged=true;
				break;
			}

			func(0)=f2;
			func(1)=f3;

			double sigma2 = std::sqrt( tra(0)*tra(0) + (Css*tra(1))*(Css*tra(1)) );
			double taudl2 = - tule(1,1) * rdum1;
			double psidk2 = - (tanpsiu-tanpsi0)/c0 * dyi2dk;
			double sigdl2 = tule(0,0) * (dl2*psidk2+tanpsi) * (tra(0)/unconf-1.0) * rdum;
			double phidk2 = - ( tanphiu - tanphi0 ) / c0 * dyi2dk;
			double taudl3 = -2.0 * Css * tule(1,1) * tra(1) * rdum1;
			double sigdl3 = -2.0 * tule(0,0) * tra(0) * rdum;
			double k3dsig = 2.0 * dl3 * tra(0) / sigma2;
			double k3dtau = 2.0 * dl3 * tra(1) * Css*Css / sigma2;
			dk3dl3 = 2.0 * sigma2;

			jacob(0,0) = taudl2 + sigdl2 * tanphi + tra(0) * phidk2 - dyi2dk;
			jacob(0,1) = taudl3 + sigdl3 * tanphi;
			jacob(1,0) = 2.0*tra(0)*sigdl2 + 2.0*Css*tra(1)*taudl2 -
						 2.0*fc*dyi3dk * (k3dsig*sigdl2 + k3dtau*taudl2);
 			jacob(1,1) = 2.0*tra(0)*sigdl3 + 2.0*Css*tra(1)*taudl3 -
						 2.0*fc*dyi3dk * (dk3dl3 + k3dsig*sigdl3 + k3dtau*taudl3);

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			jacob(0,0) -= d.eta_over_dt;
			jacob(1,1) -= d.eta_over_dt*norm_gradf3;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			MathUtils<double>::InvertMatrix2(jacob,invmat,dummy_det);
			x -= prod(invmat,func);
			if(x(1) <= 0.0) {
				x(1)=0.0;
				x(0)=dl2-f2/jacob(0,0);
			}

		}

		// check convergence
		if(!converged) {
			std::stringstream ss;
			ss << "Return mapping - Corner II & III : maximum iterations reached. F = " << f << std::endl;
			std::cout << ss.str();
		}

		// updated equivalent plastic strains
		upeq2 = rdum2;
		upeq3 = rdum3;

		// correct the sign of tau
		tra(1) = tra(1)*sign_tau;

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		// dup = dl2*[tanpsi;sign_tau] + dl3*gradf3;
		// Note: gradf3 is the normal to the cap surface defined on the positive
		// part of tau. so the second component of gradf3 should be
		// modified with the sign of tau!
		dup(0) = dl2*tanpsi + dl3*gradf3(0);
		dup(1) = dl2*sign_tau + dl3*gradf3(1)*sign_tau;
		dup3 = dl3*gradf3(0);
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// Consistent Tangent
		// Adjust the dilatancy angle to zero for tra(1) higher than unconf
		if(tra(0) <= unconf) {
			tanpsi0 = 0.0;
			tanpsiu = 0.0;
		}
		double tanps1 = tanpsi0 + ( tanpsiu - tanpsi0 ) * ( c0 - c ) / c0;
		tanpsi = tanps1 * ( 1.0 - tra(0) / unconf );
		array_1d<double,2> gradf2;
		gradf2(0)=tanphi; gradf2(1)=sign_tau;
		array_1d<double,2> gradg2;
		gradg2(0)=tanpsi; gradg2(1)=sign_tau;
		noalias(gradf3) = prod(p,tra);
		array_1d<double,2> gradg3 = gradf3;

		// correct gradg2
		gradg2(0) = gradg2(0) - dl2*(tanpsiu-tanpsi0)/c0*dyi2dk*(1.0-tra(0)/unconf);

		// correct gradf3
		double df3dk3 = -2.0 * fc * dyi3dk;
		dk3dl3 = norm_2(gradf3);
		array_1d<double,2> dk3dsg = dl3/dk3dl3*gradf3;
		gradf3 += df3dk3*dk3dsg;

		// correct [d] for [h]
		Matrix h(2,2,0.0);
		h(0,0) = 1.0 / (1.0/tule(0,0) + dl2*tanps1/unconf + dl3*p(0,0));
		h(1,1) = 1.0 / (1.0/tule(1,1) + dl3*p(0,0));

		// fill matrices [u], [v], [e]
		Matrix u(2,2);
		u(0,0) = gradg2(0);
		u(1,0) = gradg2(1);
		u(0,1) = gradg3(0);
		u(1,1) = gradg3(1);
		Matrix v(2,2);
		v(0,0) = gradf2(0);
		v(1,0) = gradf2(1);
		v(0,1) = gradf3(0);
		v(1,1) = gradf3(1);
		Matrix e(2,2);
		e(0,0) = - ( 1.0 + tra(0) * ( tanphiu - tanphi0 ) / c0 ) * dyi2dk;
		e(0,1) = 0.0;
		e(1,0) = 0.0;
		e(1,1) = df3dk3 * dk3dl3;

		// tangent matrix
		// tunl = h-h*u*inv(v'*h*u-e)*v'*h;

		Matrix aux1( prod(trans(v),h) ); // v'*h
		Matrix aux2( prod(aux1,u) ); // v'*h*u
		aux2 -= e; // // v'*h*u-e
		MathUtils<double>::InvertMatrix2(aux2,invmat,dummy_det); // invmat=inv(v'*h*u-e)
		noalias(aux2) = prod(invmat,aux1);// aux2 = inv(v'*h*u-e)*v'*h;
		noalias(aux1) = prod(u,aux2); // aux1 = u*inv(v'*h*u-e)*v'*h;
		noalias(tunl) = h; // tunl = h;
		noalias(tunl) -= prod(h,aux1); // tunl = h-h*u*inv(v'*h*u-e)*v'*h;

		// return
		return converged;
	}

//#define CHULAFUN(X)( X>0.0 ? 0.0 : 1.0 )
#define CHULAFUN(X)(1.0)
#define RMAP_TOL 1.0e-8
#define RMAP_ABS_TOL 1.0e-6

	void PlasticDamageInterface2DLaw::CalculateAll(Parameters& rValues)
	{
		// get some references
		const ProcessInfo&  pinfo = rValues.GetProcessInfo();
		const GeometryType& geom  = rValues.GetElementGeometry();
		const Properties&   props = rValues.GetMaterialProperties();

#ifdef PDINTERF_2D_TEST_ELASTIC
		Matrix& KKK = rValues.GetConstitutiveMatrix();
		if(KKK.size1() != 2 || KKK.size2() != 2)
			KKK.resize(2,2,false);
		KKK.clear();
		KKK(0,0) = props[TANGENTIAL_STIFFNESS];
		KKK(1,1) = props[NORMAL_STIFFNESS];
		rValues.GetStressVector() = prod(KKK, rValues.GetStrainVector());
		return;
#endif // PDINTERF_2D_TEST_ELASTIC

		// input and output
		Vector strain = rValues.GetStrainVector();
		double temp = strain(0); strain(0)=strain(1); strain(1)=temp;

		Vector       tra    = Vector(2);
		Matrix       tunl   = Matrix(2,2);

		if(tra.size() != 2) tra.resize(2,false);
		if(tunl.size1() != 2 || tunl.size2() != 2) tunl.resize(2,2,false);
		tra.clear();
		tunl.clear();

		// initialize calculation data
		CalculationData data;
		this->InitializeData(rValues, data);

		// recover converged internal variables
		double             lambda1 = m_lambda_converged(0);
		double             lambda2 = m_lambda_converged(1);
		double             lambda3 = m_lambda_converged(2);
		array_1d<double,2> up      = m_up_converged;
		double             up1     = m_up_aux_converged(0);
		double             up3     = m_up_aux_converged(1);

		// reset the error code
		m_error_code = 0.0;

		// some auxiliary variables
		double rdum;
		Matrix p(2,2,0.0);
		double Css = data.Css;
		p(0,0) = 2.0; p(1,1) = 2.0*Css;
		array_1d<double,2> pt;

		// elastic stress and tangent
		Matrix tule(2,2,0.0);
		tule(0,0) = data.kn;
		tule(1,1) = data.kt;
		array_1d<double,2> uela;
		uela(0) = strain(0)-up(0);
		uela(1) = strain(1)-up(1);
		Vector trat( prod(tule,uela) );

		// first damaged trial
		double damage_T = 0.0;
		double damage_C = 0.0;
		double T1; // to recover damaged normal stress prediction
		double u1trial=0.0;
		if(data.use_damage) {
			CalculateDamage(data,lambda1,lambda2,lambda3,up1,up3,damage_T,damage_C);
			u1trial = strain(0)-up(0)+up1+up3;
			if( u1trial > 0.0 )
				trat(0)=(1.0-damage_T)*data.kn*u1trial;
			else
				trat(0)=(1.0-damage_C)*data.kn*u1trial;
			T1 = trat(0);
		}

		// save the sign of the shear stress, and make it positive for the integration algorithm.
		double sign_tau = 1.0;
		if(trat(1)<0.0)
			sign_tau = -1.0;
		trat(1)  = std::abs(trat(1));

		// plastic state at trial stress
		double upeq, ft, c, tanphi, fc;
		upeq = lambda1 + data.alpha * lambda2;
		HardeningModeI(data,upeq,ft,rdum);
		upeq = lambda1 / data.alpha + lambda2;
		HardeningModeII(data,upeq,c,rdum,tanphi);
		upeq = lambda3;
		HardeningModeIII(data,upeq,fc,rdum);
		double ft0    = ft;
		double c0     = c;
		double tanph0 = tanphi;
		double fc0    = fc;

		// trial yield functions
		double f1 = trat(0) - ft;
		double f2 = trat(1) + trat(0)*tanphi - c;
		noalias(pt) = prod(p,trat);
		double f3 = 0.5 * inner_prod(trat,pt) - fc*fc;
		f3 *= CHULAFUN(trat(0));
		double f30 = f3;

		// compute a tolerance parameter for the evaluation of yielding
		double tolarance_0 = RMAP_TOL;
		double abs_tolerance = RMAP_ABS_TOL;
		double favg        = std::sqrt( std::pow(std::max(f1,0.0),2) +
										std::pow(std::max(f2,0.0),2) +
										std::pow(std::max(f3,0.0),2) );
		double tolerance   = std::max(tolarance_0*favg, abs_tolerance);
		double tolerance3  = tolerance;

		// trial active surfaces
		bool tensile = (f1 > tolerance);
		bool shear   = (f2 > tolerance);
		bool cap     = (f3 > tolerance3);

		// ================ return mapping for multisurface plasticity ================

		int yield_mode = 0;
		double dl1 = 0.0;
		double dl2 = 0.0;
		double dl3 = 0.0;
		array_1d<double,2> dup;
		dup.clear();
		double dup1 = 0.0;
		double dup3 = 0.0;
		double ovs1 = 0.0;
		double ovs2 = 0.0;
		double ovs3 = 0.0;

		// notes on yield_mode:
		// 0 = elastic
		// 1 = tensile
		// 2 = shear
		// 3 = compression
		// 4 = corner tensile - shear
		// 5 = corner compression - shear

		if(tensile || shear || cap)
		{
			// test 01: assume that only the tensile surface is active
			if(tensile) {
				if(data.use_damage) {
					// for the return mapping we need to include the whole plastic strain
					trat(0) = data.kn*(strain(0)-up(0)+up3);
				}
				// return mapping on tensile surface (mode I)
				lambda1 = m_lambda_converged(0);
				lambda2 = m_lambda_converged(1);
				lambda3 = m_lambda_converged(2);
				upeq = lambda1 + data.alpha * lambda2;
				bool rmres = RMapModeI(data,tule,trat,tolerance,sign_tau,tra,tunl,
									   ft,upeq,dup,dup1,dl1,ovs1);
				lambda1 = upeq - data.alpha * lambda2;
				// check other modes
				f2 = std::abs(tra(1)) + tra(0)*tanph0 - c0;
				noalias(pt) = prod(p,tra);
				f3 = 0.5 * inner_prod(tra,pt) - fc0*fc0;
				f3 *= CHULAFUN(tra(0));
				if(f2 < tolerance && f3 < tolerance3 && rmres) {
					yield_mode = 1;
					// updated plastic displacement vector
					up += dup;
					double BETA_T = GetBetaTension(data,lambda1,lambda2,up1+dup1);
					up1 += dup1*BETA_T;
				}
				else {
					if(data.use_damage) {
						// restore with the damaged prediction
						trat(0) = T1;
					}
				}
			}

			// test 02: assume that only the shear surface is active
			if(shear && yield_mode==0) {
				// return mapping on shear surface (mode II)
				lambda1 = m_lambda_converged(0);
				lambda2 = m_lambda_converged(1);
				lambda3 = m_lambda_converged(2);
				upeq = lambda1 / data.alpha + lambda2;
				bool rmres = RMapModeII(data,tule,trat,tolerance,sign_tau,tra,tunl,
										c,tanphi,upeq,dup,dl2,ovs2);
				lambda2 = upeq - lambda1 / data.alpha;
				// check other modes
				f1 = tra(0)-ft0;
				noalias(pt) = prod(p,tra);
				f3 = 0.5 * inner_prod(tra,pt) - fc0*fc0;
				f3 *= CHULAFUN(tra(0));
				if(f1 < tolerance && f3 < tolerance3 && rmres) {
					yield_mode = 2;
					// updated plastic displacement vector
					up += dup;
				}
			}

			// test 03: assume that only the cap surface is active
			if(cap && yield_mode==0) {
				if(data.use_damage) {
					// for the return mapping we need to include the whole plastic strain
					trat(0) = data.kn*(strain(0)-up(0)+up1);
				}
				// return mapping on cap surface (mode III)
				lambda1 = m_lambda_converged(0);
				lambda2 = m_lambda_converged(1);
				lambda3 = m_lambda_converged(2);
				upeq = lambda3;
				bool rmres = RMapModeIII(data,tule,trat,tolerance3,sign_tau,tra,tunl,
										 fc,upeq,dup,dup3,dl3,ovs3);
				lambda3 = upeq;
				// check other modes
				f1 = tra(0)-ft0;
				f2 = std::abs(tra(1)) + tra(0)*tanph0 - c0;
				if(f1 < tolerance && f2 < tolerance3 && rmres) {
					yield_mode = 3;
					// updated plastic displacement vector
					up += dup;
					double BETA_C = GetBetaCompression(data,lambda3,up3+dup3);
					up3 += dup3*BETA_C;
				}
				else {
					if(data.use_damage) {
						// restore with the damaged prediction
						trat(0) = T1;
					}
				}
			}

			// test 04:  both tensile and shear surfaces are active
			if((tensile || shear) && yield_mode==0) {
				if(data.use_damage) {
					// for the return mapping we need to include the whole plastic strain
					trat(0) = data.kn*(strain(0)-up(0)+up3);
				}
				// return mapping to corner (mode I and II)
				lambda1 = m_lambda_converged(0);
				lambda2 = m_lambda_converged(1);
				lambda3 = m_lambda_converged(2);
				bool rmres = RMapCorner_I_II(data,tule,trat,tolerance,sign_tau,tra,tunl,
											 ft,c,tanphi,lambda1,lambda2,dup,dup1,dl1,dl2,ovs1,ovs2);
				// check other modes
				noalias(pt) = prod(p,tra);
				f3 = 0.5 * inner_prod(tra,pt) - fc0*fc0;
				f3 *= CHULAFUN(tra(0));
				if(f3 < tolerance3 && dl1>0.0 && dl2>0.0 && rmres) {
					yield_mode = 4;
					// updated plastic displacement vector
					up += dup;;
					double BETA_T = GetBetaTension(data,lambda1,lambda2,up1+dup1);
					up1 += dup1*BETA_T;
				}
				else {
					if(data.use_damage) {
						// restore with the damaged prediction
						trat(0) = T1;
					}
				}
			}

			// test 05: both shear and cap surfaces are active
			if((shear || cap) && yield_mode==0) {
				if(data.use_damage) {
					// for the return mapping we need to include the whole plastic strain
					trat(0) = data.kn*(strain(0)-up(0)+up1);
				}
				lambda1 = m_lambda_converged(0);
				lambda2 = m_lambda_converged(1);
				lambda3 = m_lambda_converged(2);
				double upeq2 = lambda1 / data.alpha + lambda2;
				double upeq3 = lambda3;
				bool rmres = RMapCorner_II_III(data,tule,trat,tolerance3,sign_tau,tra,tunl,
											   fc,c,tanphi,upeq2,upeq3,dup,dup3,dl2,dl3,ovs2,ovs3);
				lambda2 = upeq2 - lambda1 / data.alpha;
				lambda3 = upeq3;
				// check other modes
				f1 = tra(0) - ft0;
				if(f1 < tolerance && dl2>0.0 && dl3>0.0 && rmres) {
					yield_mode = 5;
					// updated plastic displacement vector
					up += dup;
					double BETA_C = GetBetaCompression(data,lambda3,up3+dup3);
					up3 += dup3*BETA_C;
				}
				else {
					if(data.use_damage) {
						// restore with the damaged prediction
						trat(0) = T1;
					}
				}
			}

			// failed
			if(yield_mode==0) {
				std::stringstream ss;
				ss << "Interface return mapping failed: [" << std::boolalpha << tensile
					<< ", " << std::boolalpha << shear
					<< ", " << std::boolalpha << cap << "]" << std::endl;
				if(cap && !tensile && !shear)
				{
					ss << "F1: " << f1 << ", F2: " << f2 << std::endl;
				}
				std::cout << ss.str();
				m_error_code = -1.0;
			}
		}
		else
		{
			noalias(tra)   = trat;
			tra(1)         = tra(1)*sign_tau;
			noalias(tunl)  = tule;
		}

		// apply damage
		if(data.use_damage) {
			if(yield_mode==0) {
				CalculateDamage(data,lambda1,lambda2,lambda3,up1,up3,damage_T,damage_C);
				double kdt = (1-damage_T)*data.kn;
				double kdc = (1-damage_C)*data.kn;
				double epp = strain(0)-up(0)+up1+up3;
				double snd = kdt*epp;
				if(snd>0.0) {
					tra(0)    = snd;
					tunl(0,0) = kdt;
				}
				else {
					snd       = kdc*epp;
					tra(0)    = snd;
					tunl(0,0) = kdc;
				}
			}
			else {
				if(dl1>0.0 || dl2>0.0) {
					CalculateDamage(data,lambda1,lambda2,lambda3,up1,up3,damage_T,damage_C);
					if(dl1>0.0 && (yield_mode==1 || yield_mode==4)) {
						double ovs = dl1*data.eta_over_dt;
						ovs = std::min(tra(0),ovs);
						tra(0)=tra(0)-ovs+(1.0-damage_T)*ovs;
					}
					if(dl2>0.0 && (yield_mode==2 || yield_mode==4)) {
						double ovs = dl2*data.eta_over_dt;
						ovs = std::min(ovs,std::abs(tra(1)));
						double damage_S = 1.0-c/data.c;
						tra(1)=sign_tau*(std::abs(tra(1))-ovs+(1.0-damage_S)*ovs);
					}
				}
			}
		}

		// store trial internal variables
		if(m_error_code == 0.0)
		{
			m_lambda(0) = lambda1;
			m_lambda(1) = lambda2;
			m_lambda(2) = lambda3;
			m_up		= up;
			m_up_aux(0) = up1;
			m_up_aux(1) = up3;

			// store damage
			m_damage_T = damage_T;
			m_damage_C = damage_C;
		}

		// save results
		Vector& SS = rValues.GetStressVector();
		if(SS.size() != 2) SS.resize(2,false);
		SS(0) = tra(1); SS(1) = tra(0);
		Matrix& CC = rValues.GetConstitutiveMatrix();
		if(CC.size1() != 2 || CC.size2() != 2)
			CC.resize(2,2,false);
		CC(0,0) = tunl(1,1); CC(0,1) = tunl(1,0);
		CC(1,0) = tunl(0,1); CC(1,1) = tunl(0,0);
	}

} /* namespace Kratos.*/

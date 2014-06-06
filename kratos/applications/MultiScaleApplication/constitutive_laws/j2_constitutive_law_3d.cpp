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

#include "j2_constitutive_law_3d.h"
#include "multiscale_application.h"

namespace Kratos
{

    J2ConstitutiveLaw3D::J2ConstitutiveLaw3D() 
        : ConstitutiveLaw()
		, m_initialized(false)
    {
    }
 
    ConstitutiveLaw::Pointer J2ConstitutiveLaw3D::Clone() const
    {
        return ConstitutiveLaw::Pointer(new J2ConstitutiveLaw3D());
    }

    J2ConstitutiveLaw3D::SizeType J2ConstitutiveLaw3D::WorkingSpaceDimension()
    {
        return 3;
    }

    J2ConstitutiveLaw3D::SizeType J2ConstitutiveLaw3D::GetStrainSize()
    {
        return 6;
    }

    bool J2ConstitutiveLaw3D::Has(const Variable<double>& rThisVariable)
    {
		if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN) return true;
        return false;
    }
    
    bool J2ConstitutiveLaw3D::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    bool J2ConstitutiveLaw3D::Has(const Variable<Matrix>& rThisVariable)
    {
		if(rThisVariable == PLASTIC_STRAIN_TENSOR) return true;
        return false;
    }

    bool J2ConstitutiveLaw3D::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }
    
    bool J2ConstitutiveLaw3D::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& J2ConstitutiveLaw3D::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
		rValue = 0.0;
		if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN) 
			rValue = m_xi_n1;
        return rValue;
    }

    Vector& J2ConstitutiveLaw3D::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
        return rValue;
    }

    Matrix& J2ConstitutiveLaw3D::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
		if(rThisVariable == PLASTIC_STRAIN_TENSOR) {
			if(rValue.size1() != 1 || rValue.size2() != 6) rValue.resize(1, 6, false);
			for(int i = 0; i < 6; i++)
				rValue(0,i) = m_eps_p_n1(i);
		}
        return rValue;
    }

    array_1d<double, 3 > & J2ConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & J2ConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

     void J2ConstitutiveLaw3D::SetValue(const Variable<double>& rVariable,
                                         const double& rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::SetValue(const Variable<Vector >& rVariable,
                                        const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::SetValue(const Variable<Matrix >& rVariable,
                                        const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                         const array_1d<double, 3 > & rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                        const array_1d<double, 6 > & rValue,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool J2ConstitutiveLaw3D::ValidateInput(const Properties& rMaterialProperties)
    {
		if(!rMaterialProperties.Has(YOUNG_MODULUS)) return false;
		if(!rMaterialProperties.Has(POISSON_RATIO)) return false;
		if(!rMaterialProperties.Has(DENSITY)) return false;
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING)) return false;
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING_EXPONENT)) return false;
		if(!rMaterialProperties.Has(KINEMATIC_HARDENING)) return false;
		if(!rMaterialProperties.Has(YIELD_STRESS)) return false;
		if(!rMaterialProperties.Has(YIELD_STRESS_INFINITY)) return false;
        return true;
    }

    J2ConstitutiveLaw3D::StrainMeasure J2ConstitutiveLaw3D::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }
    
    J2ConstitutiveLaw3D::StressMeasure J2ConstitutiveLaw3D::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool J2ConstitutiveLaw3D::IsIncremental()
    {
        return false;
    }

    void J2ConstitutiveLaw3D::InitializeMaterial(const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues)
    {
		if(!m_initialized)
		{
			m_xi_n = 0.0;
			m_xi_n1 = 0.0;
			m_eps_p_n.clear();
			m_eps_p_n1.clear();
			m_initialized = true;
		}
    }

    void J2ConstitutiveLaw3D::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                      const GeometryType& rElementGeometry,
                                                      const Vector& rShapeFunctionsValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
    {
		m_xi_n1 = m_xi_n;
		m_eps_p_n1 = m_eps_p_n;
    }

    void J2ConstitutiveLaw3D::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                    const GeometryType& rElementGeometry,
                                                    const Vector& rShapeFunctionsValues,
                                                    const ProcessInfo& rCurrentProcessInfo)
    {
		m_xi_n = m_xi_n1;
		m_eps_p_n = m_eps_p_n1;
    }

    void J2ConstitutiveLaw3D::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void J2ConstitutiveLaw3D::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
		/*CalculateNumericalVersion(rValues);
		return;*/

		// get some references
		const Properties& props = rValues.GetMaterialProperties();
		const Vector& eps = rValues.GetStrainVector();

		// get the options
		Flags& options = rValues.GetOptions();
		bool compute_tangent = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
		compute_tangent = true;
		// get a reference of the output stress vector and tangent operator
		// and initialize them
		Vector& stressVector = rValues.GetStressVector();
		Matrix& tangentMatrix = rValues.GetConstitutiveMatrix();
		
		if(stressVector.size() != 6)
			stressVector.resize(6, false);
		noalias(stressVector) = ZeroVector(6);

		if(compute_tangent) {
			if(tangentMatrix.size1() != 6 || tangentMatrix.size2() != 6)
				tangentMatrix.resize(6, 6, false);
			noalias(tangentMatrix) = ZeroMatrix(6, 6);
		}

        // get material parameters
		double E  = props[YOUNG_MODULUS];
		double ni = props[POISSON_RATIO];
		double Hi = props[ISOTROPIC_HARDENING];
		double Hk = props[KINEMATIC_HARDENING];
		double y0 = props[YIELD_STRESS];

		// elastic constants
		double G = E/(2.0*(1.0+ni));
		double K = E/(3.0*(1.0-2.0*ni));

		// trial state
		m_eps_p_n1 = m_eps_p_n;
		m_xi_n1 = m_xi_n;

		// volumetric strain
		double eps_vol = eps(0) + eps(1) + eps(2);

		// deviatoric strains
		array_1d<double, 6> eps_dev;
		eps_dev(0) = eps(0) - eps_vol / 3.0;
		eps_dev(1) = eps(1) - eps_vol / 3.0;
		eps_dev(2) = eps(2) - eps_vol / 3.0;
		eps_dev(3) = eps(3)/2.0;
		eps_dev(4) = eps(4)/2.0;
		eps_dev(5) = eps(5)/2.0;

		// trial deviatoric stress
		array_1d<double, 6> sig_dev;
		noalias( sig_dev ) = 2.0 * G * (eps_dev  - m_eps_p_n1);

		// trial relative stress and its norm
		array_1d<double, 6> sig_rel;
		noalias( sig_rel ) = sig_dev - 2.0/3.0*Hk * m_eps_p_n1;
		double norm_sig_rel = std::sqrt(
			                  sig_rel(0)*sig_rel(0) + 
			                  sig_rel(1)*sig_rel(1) + 
							  sig_rel(2)*sig_rel(2) + 
							  2.0*sig_rel(3)*sig_rel(3) + 
							  2.0*sig_rel(4)*sig_rel(4) + 
							  2.0*sig_rel(5)*sig_rel(5)
							  );

		// isotropic hardening variable
		double q = -Hi * m_xi_n1;

		// initialize a parameter for the tangent operator here to 1
		double delta1 = 1.0;

		// Yield function
		double F = norm_sig_rel - std::sqrt(2.0/3.0)*(y0 - q);

		// Plastic corrector
		if( F > 0.0 )
		{
			double dgamma = F / (2.0*G + 2.0/3.0*(Hi + Hk));
			m_xi_n1 += dgamma * std::sqrt(2.0/3.0);

			array_1d<double, 6> N;
			noalias(N) = sig_rel / norm_sig_rel;

			noalias(m_eps_p_n1) += dgamma * N;
			noalias(sig_dev) -= dgamma*2.0*G*N;

			if(compute_tangent) {
				delta1 = 1.0-2.0*G*dgamma/norm_sig_rel;
				double delta2 = 1.0/(1.0+1.0/(3.0*G)*(Hi+Hk))-(1.0-delta1);
				noalias( tangentMatrix ) -= 2.0*G*delta2*outer_prod(N,N);
			}
		}

		// compute stresses
		double pressure = K * eps_vol;
		stressVector(0) = sig_dev(0) + pressure;
		stressVector(1) = sig_dev(1) + pressure;
		stressVector(2) = sig_dev(2) + pressure;
		stressVector(3) = sig_dev(3);
		stressVector(4) = sig_dev(4);
		stressVector(5) = sig_dev(5);

		// compute the tangent operator
		if(compute_tangent)
		{
			double D = 2.0*G*delta1;
			double K1 = K + 2.0/3.0*D;
			double K2 = K - D/3.0;
			double K3 = D/2.0;
			tangentMatrix(0,0) += K1;
			tangentMatrix(1,1) += K1;
			tangentMatrix(2,2) += K1;
			tangentMatrix(3,3) += K3;
			tangentMatrix(4,4) += K3;
			tangentMatrix(5,5) += K3;
			tangentMatrix(0,1) += K2; tangentMatrix(0,2) += K2; tangentMatrix(1,2) += K2;
			tangentMatrix(1,0) += K2; tangentMatrix(2,0) += K2; tangentMatrix(2,1) += K2;
		}
    }

    void J2ConstitutiveLaw3D::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy( rValues );
    }

    void J2ConstitutiveLaw3D::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {
    }

    void J2ConstitutiveLaw3D::ResetMaterial(const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues)
    {
        if(m_initialized)
		{
			m_xi_n = 0.0;
			m_xi_n1 = 0.0;
			m_eps_p_n.clear();
			m_eps_p_n1.clear();
			m_initialized = false;
		}
    }

    void J2ConstitutiveLaw3D::GetLawFeatures(Features& rFeatures)
    {
        rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ISOTROPIC );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
        rFeatures.mSpaceDimension = WorkingSpaceDimension();
		rFeatures.mStrainSize = GetStrainSize();
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    }

    int J2ConstitutiveLaw3D::Check(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
		if(!rMaterialProperties.Has(YOUNG_MODULUS)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YOUNG_MODULUS", "");
		}
		if(!rMaterialProperties.Has(POISSON_RATIO)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing POISSON_RATIO", "");
		}
		if(!rMaterialProperties.Has(DENSITY)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing DENSITY", "");
		}
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing ISOTROPIC_HARDENING", "");
		}
		if(!rMaterialProperties.Has(KINEMATIC_HARDENING)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing KINEMATIC_HARDENING", "");
		}
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING_EXPONENT)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing ISOTROPIC_HARDENING_EXPONENT", "");
		}
		if(!rMaterialProperties.Has(YIELD_STRESS)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YIELD_STRESS", "");
		}
		if(!rMaterialProperties.Has(YIELD_STRESS_INFINITY)) {
			KRATOS_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YIELD_STRESS_INFINITY", "");
		}
		return 0;
        KRATOS_CATCH("");
    }


	void J2ConstitutiveLaw3D::CalculateNumericalVersion(Parameters& rValues)
    {
		const Vector& eps = rValues.GetStrainVector();

		Flags& options = rValues.GetOptions();
		bool compute_tangent = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
		compute_tangent = true;

		Vector& stressVector = rValues.GetStressVector();
		Matrix& tangentMatrix = rValues.GetConstitutiveMatrix();
		
		if(stressVector.size() != 6)
			stressVector.resize(6, false);
		noalias(stressVector) = ZeroVector(6);

		if(compute_tangent) {
			if(tangentMatrix.size1() != 6 || tangentMatrix.size2() != 6)
				tangentMatrix.resize(6, 6, false);
			noalias(tangentMatrix) = ZeroMatrix(6, 6);
		}

        
		CalculateStressVector(rValues, eps, stressVector);
		
		if(compute_tangent)
		{
			double save_xi = m_xi_n1;
			array_1d<double,6> save_ep = m_eps_p_n1;

			CalculateNumericalTangent(rValues);

			m_xi_n1 = save_xi;
			m_eps_p_n1 = save_ep;
		}
	}

	void J2ConstitutiveLaw3D::CalculateStressVector(Parameters& rValues, const Vector& eps, Vector& stressVector)
	{
		// get some references
		const Properties& props = rValues.GetMaterialProperties();

        // get material parameters
		double E  = props[YOUNG_MODULUS];
		double ni = props[POISSON_RATIO];
		double Hi = props[ISOTROPIC_HARDENING];
		double Hk = props[KINEMATIC_HARDENING];
		double y0 = props[YIELD_STRESS];

		// elastic constants
		double G = E/(2.0*(1.0+ni));
		double K = E/(3.0*(1.0-2.0*ni));

		// trial state
		m_eps_p_n1 = m_eps_p_n;
		m_xi_n1 = m_xi_n;

		// volumetric strain
		double eps_vol = eps(0) + eps(1) + eps(2);

		// deviatoric strains
		array_1d<double, 6> eps_dev;
		eps_dev(0) = eps(0) - eps_vol / 3.0;
		eps_dev(1) = eps(1) - eps_vol / 3.0;
		eps_dev(2) = eps(2) - eps_vol / 3.0;
		eps_dev(3) = eps(3)/2.0;
		eps_dev(4) = eps(4)/2.0;
		eps_dev(5) = eps(5)/2.0;

		// trial deviatoric stress
		array_1d<double, 6> sig_dev;
		noalias( sig_dev ) = 2.0 * G * (eps_dev  - m_eps_p_n1);

		// trial relative stress and its norm
		array_1d<double, 6> sig_rel;
		noalias( sig_rel ) = sig_dev - 2.0/3.0*Hk * m_eps_p_n1;
		double norm_sig_rel = std::sqrt(
			                  sig_rel(0)*sig_rel(0) + 
			                  sig_rel(1)*sig_rel(1) + 
							  sig_rel(2)*sig_rel(2) + 
							  2.0*sig_rel(3)*sig_rel(3) + 
							  2.0*sig_rel(4)*sig_rel(4) + 
							  2.0*sig_rel(5)*sig_rel(5)
							  );

		// isotropic hardening variable
		double q = -Hi * m_xi_n1;

		// Yield function
		double F = norm_sig_rel - std::sqrt(2.0/3.0)*(y0 - q);

		// Plastic corrector
		if( F > 0.0 )
		{
			double dgamma = F / (2.0*G + 2.0/3.0*(Hi + Hk));
			m_xi_n1 += dgamma * std::sqrt(2.0/3.0);

			array_1d<double, 6> N;
			noalias(N) = sig_rel / norm_sig_rel;

			noalias(m_eps_p_n1) += dgamma * N;
			noalias(sig_dev) -= dgamma*2.0*G*N;
		}

		// compute stresses
		double pressure = K * eps_vol;
		stressVector(0) = sig_dev(0) + pressure;
		stressVector(1) = sig_dev(1) + pressure;
		stressVector(2) = sig_dev(2) + pressure;
		stressVector(3) = sig_dev(3);
		stressVector(4) = sig_dev(4);
		stressVector(5) = sig_dev(5);
	}

	void J2ConstitutiveLaw3D::CalculateNumericalTangent(Parameters& rValues)
	{
		const Vector& strainVector = rValues.GetStrainVector();
		const Vector& stressVector = rValues.GetStressVector();

		Matrix& tangentMatrix = rValues.GetConstitutiveMatrix();
		if(tangentMatrix.size1() != 6 || tangentMatrix.size2() != 6)
			tangentMatrix.resize(6, 6, false);
		noalias(tangentMatrix) = ZeroMatrix(6, 6);

		Vector perturbedStrainVector(6);
		Vector stressPerturbation(6);
		double perturbation = 1.0E-10;

		for(unsigned int j = 0; j < 6; j++)
		{
			noalias(perturbedStrainVector) = strainVector;

			double delta_strain = perturbation;
			if(perturbedStrainVector(j) < 0.0)
				delta_strain = -delta_strain;
			perturbedStrainVector(j) += delta_strain;

			CalculateStressVector(rValues, perturbedStrainVector, stressPerturbation);

			stressPerturbation -= stressVector;

			for(unsigned int i = 0; i < 6; i++)
				tangentMatrix(i, j) = stressPerturbation(i) / delta_strain;
		}
	}

} /* namespace Kratos.*/

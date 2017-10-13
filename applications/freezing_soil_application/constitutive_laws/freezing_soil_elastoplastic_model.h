/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
 *   Last Modified by:    $Author: janosch $
 *   Date:                $Date: 2008-01-25 08:37:31 $
 *   Revision:            $Revision: 1.10 $
 *
 * ***********************************************************/

#if !defined(KRATOS_FREEZING_SOIL_ELASTOPLASTIC_MODEL_H_INCLUDED )
#define  KRATOS_FREEZING_SOIL_ELASTOPLASTIC_MODEL_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
/**
 * Elasto-plastic consititutive model for freezing soils as integration of CASM and enhanced BBM
*/

class FreezingSoilElastoplasticModel : public ConstitutiveLaw
{

    public:

        typedef ConstitutiveLaw BaseType;

        /**
         * Counted pointer of FreezingSoilElastoplasticModel
         */
        KRATOS_CLASS_POINTER_DEFINITION ( FreezingSoilElastoplasticModel );

        /**
         * Life Cycle
         */
        FreezingSoilElastoplasticModel();

        /**
         * Destructor.
         */
        virtual ~FreezingSoilElastoplasticModel();

        /**
         * Clone function
         * @return a pointer to a new instance of this constitutive law
         */
        virtual boost::shared_ptr<BaseType> Clone() const
        {
            boost::shared_ptr<BaseType> p_clone ( new FreezingSoilElastoplasticModel() );
            return p_clone;
        }

        /**
         * @return the working space dimension of the current constitutive law
         */
        virtual SizeType WorkingSpaceDimension()
        {
            return ( 3 );
        }

        /**
         * returns the size of the strain vector of the current constitutive law
         */
        virtual SizeType GetStrainSize()
        {
            return ( 6 );
        }

        bool Has ( const Variable<double>& rThisVariable );
        bool Has ( const Variable<Vector>& rThisVariable );
        bool Has ( const Variable<Matrix>& rThisVariable );
        bool Has ( const Variable<array_1d<double, 3 > >& rThisVariable );
        bool Has ( const Variable<array_1d<double, 6 > >& rThisVariable );

        virtual double& GetValue ( const Variable<double>& rThisVariable, double& rValue );
        virtual Vector& GetValue ( const Variable<Vector>& rThisVariable, Vector& rValue );
        virtual Matrix& GetValue ( const Variable<Matrix>& rThisVariable, Matrix& rValue );
        virtual array_1d<double, 3 > & GetValue ( const Variable<array_1d<double, 3 > >& rVariable,
                array_1d<double, 3 > & rValue );
        virtual array_1d<double, 6 > & GetValue ( const Variable<array_1d<double, 6 > >& rVariable,
                array_1d<double, 6 > & rValue );
        virtual void SetValue ( const Variable<double>& rVariable,
                                const double& Value,
                                const ProcessInfo& rCurrentProcessInfo );
        virtual void SetValue ( const Variable<Vector >& rVariable,
                                const Vector& Value, const ProcessInfo& rCurrentProcessInfo );
        virtual void SetValue ( const Variable<Matrix >& rVariable,
                                const Matrix& Value, const ProcessInfo& rCurrentProcessInfo );
        virtual void SetValue ( const Variable<array_1d<double, 3 > >& rVariable,
                                const array_1d<double, 3 > & Value,
                                const ProcessInfo& rCurrentProcessInfo );
        virtual void SetValue ( const Variable<array_1d<double, 6 > >& rVariable,
                                const array_1d<double, 6 > & Value,
                                const ProcessInfo& rCurrentProcessInfo );
        virtual bool ValidateInput ( const Properties& props );
        virtual StrainMeasure GetStrainMeasure()
        {
            return StrainMeasure_Infinitesimal;
        }

        virtual StressMeasure GetStressMeasure()
        {
            return StressMeasure_PK1;
        }

        virtual bool IsIncremental()
        {
            return false;
        }

        virtual void InitializeMaterial ( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues );
					  
        virtual void InitializeSolutionStep ( const Properties& props,
					const GeometryType& geom, //this is just to give the array of nodes
					const Vector& ShapeFunctionsValues,
					const ProcessInfo& CurrentProcessInfo );
					      
        virtual void FinalizeSolutionStep ( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues,
					const ProcessInfo& CurrentProcessInfo );
					    
        virtual void InitializeNonLinearIteration ( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues,
					const ProcessInfo& CurrentProcessInfo );
		
        virtual void CalculateMaterialResponse ( const Vector& StrainVector,
					const Matrix& DeformationGradient,
					Vector& StressVector,
					Matrix& AlgorithmicTangent,
					const ProcessInfo& CurrentProcessInfo,
					const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues,
					bool CalculateStresses = true,
					int CalculateTangent = true,
					bool ShowIterativeResults /*SaveInternalVariables*/ = true );
		
        virtual void CalculateVolumetricResponse ( const double VolumetricStrain,
					const Matrix& DeformationGradient,
					double& VolumetricStress,
					double& AlgorithmicBulk,
					const ProcessInfo& CurrentProcessInfo,
					const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues,
					bool CalculateStresses = true,
					int CalculateTangent = true,
					bool SaveInternalVariables = true )
        {
            KRATOS_THROW_ERROR ( std::logic_error, "Volumetric response is not implemented in FreezingSoilElastoplasticModel", "" );
        }

        virtual void CalculateDeviatoricResponse ( const Vector& InputVector,
					const Matrix& InternalVariables,
					Vector& StressVector,
					Matrix& AlgorithmicTangent,
					const ProcessInfo& CurrentProcessInfo,
					const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues,
					bool CalculateStresses = true,
					int CalculateTangent = true,
					bool SaveInternalVariables = true )
        {
            KRATOS_THROW_ERROR ( std::logic_error, "Deviatoric response is not implemented in FreezingSoilElastoplasticModel", "" );
        }

        virtual void ResetMaterial ( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues );

        virtual void CalculateCauchyStresses ( Vector& Cauchy_StressVector,
					const Matrix& F,
					const Vector& PK2_StressVector,
					const Vector& GreenLagrangeStrainVector )
        {
            KRATOS_THROW_ERROR ( std::logic_error, "Cauchy stresses are not implemented in FreezingSoilElastoplasticModel", "" );
        }

        /**
         * This function is designed to be called once to perform all the checks needed
         * on the input provided. Checks can be "expensive" as the function is designed
         * to catch user's errors.
         * @param props
         * @param geom
         * @param CurrentProcessInfo
         * @return
         */
        virtual int Check ( const Properties& props,
                            const GeometryType& geom,
                            const ProcessInfo& CurrentProcessInfo );

    protected:

        /**
         * Member Variables
         */
	
	double mTol, mpmin, mMth, mK, mG, mMSL, mMSC, mS, mtheta, mSf, mE, mpi, mpL, mt;
	double mnu, mMmax, mphiCS, mlambda, mkappa, me, mnn, mrr;
	double mc, mtstar, mm, mn0, mUnitRatio;
	Vector mMaterialParameters;

        Vector mOldStress;
        Vector mOldStrain;
        double mOldPreconsolidation, mPreHydrostaticPressure;
	
        Vector mCurrentStress;
        Vector mCurrentStrain;	
        double mCurrentPreconsolidation;

    private:

        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw );
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw );
        }
        
	void UpdateMS ();
	double GetLiquidSaturation ( const double t );
	void UpdateKG ( const double p );
	Matrix GetElasticTangent();
	Matrix GetElastoPlasticTangent ( double pn, double q, double pn0, Vector& stressTr, double DEpsPp, double DEpsPq );
	Vector GetResidual ( double pn, double q, double pn0, double DEpsPp, double DEpsPq );
	Matrix GetStiffness ( double pn, double q, double pn0, double DEpsPp, double DEpsPq );

	double GetYieldFunctionAndDerivatives ( double pn, double q, double pn0, int index );
	double GetPotentialFunctionDerivatives ( double pn, double q, int index );

	double GetMth ( Vector& stress );
	double Getp ( Vector& stress );
	double Getq ( Vector& stress );
	Vector GetDevStress ( Vector& stress ); 
}; // Class FreezingSoilElastoplasticModel
}  // namespace Kratos.
#endif // KRATOS_FREEZING_SOIL_ELASTOPLASTIC_MODEL_H_INCLUDED  defined

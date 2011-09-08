#if !defined(KRATOS_UMAT_H_INCLUDED )
#define  KRATOS_UMAT_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

    class Umat : public ConstitutiveLaw
    {

        public:

            // Type Definitions
            typedef ConstitutiveLaw BaseType;
            typedef array_1d<double, 81> MaterialTensorType;
            typedef boost::shared_ptr<Umat> Pointer;
            typedef array_1d<double, 3 > PlaneArrayType;
            typedef array_1d<double, 6 > SpaceArrayType;

            //zero constructor
            Umat();
            //destructor
            virtual ~Umat();

            //clone
            virtual boost::shared_ptr<BaseType> Clone() const
            {
                boost::shared_ptr<BaseType> p_clone ( new Umat() );
                return p_clone;
            }

            size_t WorkingSpaceDimension()
            {
                return 3;
            }

            size_t GetStrainSize()
            {
                return 6;
            }

            //operators
            virtual bool Has ( const Variable<double>& rThisVariable )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::Has(double)", "" );
            }

            virtual bool Has ( const Variable<Vector>& rThisVariable )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::Has(Vector)", "" );
            }

            virtual bool Has ( const Variable<Matrix>& rThisVariable )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::Has(Matrix)", "" );
            }

            virtual bool Has ( const Variable<PlaneArrayType>& rThisVariable )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::Has(2DVector)", "" );
            }

            virtual bool Has ( const Variable<SpaceArrayType>& rThisVariable )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::Has(3DVector)", "" );
            }

            virtual double& GetValue ( const Variable<double>& rThisVariable,
                                       double& rValue )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::GetValue(double)", "" );
            }

            virtual Vector& GetValue ( const Variable<Vector>& rThisVariable,
                                       Vector& rValue )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::GetValue(Vector)", "" );
            }

            virtual Matrix& GetValue ( const Variable<Matrix>& rThisVariable,
                                       Matrix& rValue )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::GetValue(Matrix)", "" );
            }

            virtual PlaneArrayType& GetValue ( const Variable<PlaneArrayType>& rVariable,
                                               PlaneArrayType& rValue )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::GetValue(2DVector)", "" );
            }

            virtual SpaceArrayType & GetValue ( const Variable<SpaceArrayType>& rVariable,
                                                SpaceArrayType& rValue )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::GetValue(3DVector)", "" );
            }


            virtual void SetValue ( const Variable<double>& rThisVariable,
                                    const double& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::SetValue(Double)", "" );
            }

            virtual void SetValue ( const Variable<Vector>& rThisVariable,
                                    const Vector& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::SetValue(Vector)", "" );
            }

            virtual void SetValue ( const Variable<Matrix>& rThisVariable,
                                    const Matrix& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::SetValue(Matrix)", "" );
            }

            virtual void SetValue ( const Variable<PlaneArrayType>& rThisVariable,
                                    const PlaneArrayType& rValue,
                                    const ProcessInfo& rCurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::SetValue(2DVector)", "" );
            }

            virtual void SetValue ( const Variable<SpaceArrayType>& rVariable,
                                    const SpaceArrayType& Value,
                                    const ProcessInfo& rCurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::SetValue(3DVector)", "" );
            }

            virtual bool ValidateInput ( const Properties& props )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::ValidateInput called", "" );
            }

            StrainMeasure GetStrainMeasure()
            {
                return StrainMeasure_Linear;
            }

            StressMeasure GetStressMeasure()
            {
                return StressMeasure_PK1;
            }

            bool IsIncremental()
            {
                return false;
            }



            //Material parameters inizialization
            virtual void InitializeMaterial ( const Properties& props,
                                              const GeometryType& geom,
                                              const Vector& ShapeFunctionsValues );


            virtual void InitializeSolutionStep ( const Properties& props,
                                                  const GeometryType& geom,
                                                  const Vector& ShapeFunctionsValues ,
                                                  const ProcessInfo& CurrentProcessInfo );

            virtual void FinalizeSolutionStep ( const Properties& props,
                                                const GeometryType& geom,
                                                const Vector& ShapeFunctionsValues ,
                                                const ProcessInfo& CurrentProcessInfo );

            virtual void InitializeNonLinearIteration ( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues,
                    const ProcessInfo& CurrentProcessInfo )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::InitializeNonLinearIteration called", "" );
            }

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
                    bool SaveInternalVariables = true );


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
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::CalculateVolumetricResponse called", "" );
            }

            virtual void CalculateDeviatoricResponse ( const Vector& StrainVector,
                    const Matrix& DeformationGradient,
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
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::CalculateDeviatoricResponse called", "" );
            }

            virtual void ResetMaterial ( const Properties& props,
                                         const GeometryType& geom,
                                         const Vector& ShapeFunctionsValues )
            {
                KRATOS_ERROR ( std::logic_error, "virtual function Umat::ResetMaterial called", "" );
            }

            virtual void CalculateCauchyStresses ( Vector& Cauchy_StressVector,
                                                   const Matrix& F,
                                                   const Vector& PK2_StressVector,
                                                   const Vector& GreenLagrangeStrainVector );

            virtual int Check ( const Properties& props,
                                const GeometryType& geom,
                                const ProcessInfo& CurrentProcessInfo );

        protected:

        private:

            //member variables for umat
            double* STRESS;
            double* STATEV; 
            double* PROPS;
            int* NSTATV;            
            int* NPROPS; 
            
            //static member variables for umat
            double* STRAN;
            double* DSTRAN;
            int* NDI;
            int* NSHR;
            int* NTENS;
                       
            
            //variable for wrapper
            int* MaterialNumber;

            //serialization

            friend class Serializer;

            virtual void save ( Serializer& rSerializer ) const
            {
                rSerializer.save ( "name", "Umat" );
                KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, ConstitutiveLaw )
            }

            virtual void load ( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, ConstitutiveLaw )
            }







    }; // Class Umat
}  // namespace Kratos.

#endif // KRATOS_UMAT_H_INCLUDED  defined 


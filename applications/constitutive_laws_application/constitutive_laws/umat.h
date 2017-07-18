#if !defined(KRATOS_UMAT_H_INCLUDED )
#define  KRATOS_UMAT_H_INCLUDED

// System includes

// External includes
//#include "boost/smart_ptr.hpp"

// Project includes
//#include "includes/define.h"
//#include "includes/variables.h"
//#include "includes/constitutive_law.h"
//#include "includes/serializer.h"
//#include "includes/ublas_interface.h"

#include "custom_constitutive/hyperelastic_3D_law.hpp"


namespace Kratos
{

class Umat : public HyperElastic3DLaw
{

public:

    // Type Definitions
    typedef array_1d<double, 81> MaterialTensorType;
    KRATOS_CLASS_POINTER_DEFINITION( Umat );
    typedef array_1d<double, 3 > PlaneArrayType;
    typedef array_1d<double, 6 > SpaceArrayType;

    //zero constructor
    Umat();
    //destructor
    virtual ~Umat();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    virtual ConstitutiveLaw::Pointer Clone() const override
    {
       Umat::Pointer p_clone(new Umat() );
       return p_clone;
    };

    size_t WorkingSpaceDimension()
    {
        return 3;
    }

    size_t GetStrainSize()
    {
        return 6;
    }

     /* This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);

    StrainMeasure GetStrainMeasure()
    {
        return StrainMeasure_Infinitesimal;
    }

    StressMeasure GetStressMeasure()
    {
        return StressMeasure_Cauchy;
    }


    //Material parameters inizialization
    virtual void InitializeMaterial ( const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues ) override;


    virtual void InitializeSolutionStep ( const Properties& props,
                                          const GeometryType& geom,
                                          const Vector& ShapeFunctionsValues ,
                                          const ProcessInfo& CurrentProcessInfo ) override;

    virtual void FinalizeSolutionStep ( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues ,
                                        const ProcessInfo& CurrentProcessInfo ) override;



    virtual void CalculateMaterialResponseCauchy( Parameters & rValues) override;

    virtual void  FinalizeMaterialResponseCauchy( Parameters & rValues) override;

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;


protected:

    void LoadPreviousInformation();
private:

    //member variables for umat
    double* STRESS;
    double* STATEV;

    double* STRESS_FINALIZED;
    double* STATEV_FINALIZED;

    double* STRAN_FINALIZED; 

    int* NSTATV;

    double* PROPS;
    int* NPROPS;

    //static member variables for umat
    //double* STRAN;
    //double* DSTRAN;


    //variable for wrapper
    int* MaterialNumber;

    //serialization

    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const
    {
        rSerializer.save ( "name", "Umat" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, HyperElastic3DLaw )
    }

    virtual void load ( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, HyperElastic3DLaw )
    }







}; // Class Umat
}  // namespace Kratos.

#endif // KRATOS_UMAT_H_INCLUDED  defined 


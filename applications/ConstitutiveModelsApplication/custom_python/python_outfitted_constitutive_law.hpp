//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_PYTHON_OUTFITTED_CONSTITUTIVE_LAW_H_INCLUDED)
#define  KRATOS_PYTHON_OUTFITTED_CONSTITUTIVE_LAW_H_INCLUDED

// System includes
#include <iostream>

// External includes
#include <cmath>
#include <boost/python.h>

// Project includes
#include "includes/constitutive_law.h"

#include "constitutive_models_application_variables.h"

namespace Kratos
{
/**
 * Defines a law supplied by a python module
 */

class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) PythonOutfittedConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /**
     * Counted pointer of PythonOutfittedConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( PythonOutfittedConstitutiveLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    PythonOutfittedConstitutiveLaw();

    /**
     * Default constructor suppling python constitutive law
     */
    PythonOutfittedConstitutiveLaw(PyObject* pPyConstitutiveLaw);


    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    PythonOutfittedConstitutiveLaw (const PythonOutfittedConstitutiveLaw& rOther);


    /**
     * Assignment operator.
     */

    //PythonOutfittedConstitutiveLaw& operator=(const PythonOutfittedConstitutiveLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~PythonOutfittedConstitutiveLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
      return boost::python::call_method<int>(mpPyConstitutiveLaw->ptr(),"WorkingSpaceDimension");
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
      return boost::python::call_method<int>(mpPyConstitutiveLaw->ptr(),"GetStrainSize");
    };


    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );


    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rProperties,
                             const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues );


    void InitializeSolutionStep( const Properties& rProperties,
                                 const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                 const Vector& rShapeFunctionsValues ,
                                 const ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep( const Properties& rProperties,
                               const GeometryType& rElementGeometry, //this is just to give the array of nodes
                               const Vector& rShapeFunctionsValues ,
                               const ProcessInfo& rCurrentProcessInfo);

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK1 (Parameters & rValues);

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters & rValues);


    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponsePK1 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponsePK2 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponseKirchhoff (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponseCauchy (Parameters & rValues);


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    Matrix mInverseTotalDeformationMatrix;  //F0

    double mTotalDeformationDet;  //detF0

    double mStrainEnergy;

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
      * Updates the material response:
      * Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void UpdateInternalVariables (Parameters & rValues);

    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    //PyObject* mpPyConstitutiveLaw;
    boost::shared_ptr<boost::python::object> mpPyConstitutiveLaw;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     * if the matrix passed is 3D is does nothing
     * if the matrix passed is bigger or smaller throws an error
     * @param rMatrix : usually the DeformationGradientF
     */
    Matrix& Transform2DTo3D (Matrix& rMatrix);

    ///@}
    ///@name Private  Access
    ///@{
    ///@}


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
	rSerializer.save("mInverseTotalDeformationMatrix",mInverseTotalDeformationMatrix);
	rSerializer.save("mTotalDeformationDet",mTotalDeformationDet);
	rSerializer.save("mStrainEnergy",mStrainEnergy);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
	rSerializer.load("mInverseTotalDeformationMatrix",mInverseTotalDeformationMatrix);
	rSerializer.load("mTotalDeformationDet",mTotalDeformationDet);
	rSerializer.load("mStrainEnergy",mStrainEnergy);
    }


    ///@}

}; // Class PythonOutfittedConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_PYTHON_OUTFITTED_CONSTITUTIVE_LAW_H_INCLUDED  defined

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined (KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
/**
 * Defines a bingham non-newtonian constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */

class Newtonian2DLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of Newtonian3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(Newtonian2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    Newtonian2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    Newtonian2DLaw (const Newtonian2DLaw& rOther);

    /**
     * Destructor.
     */
    ~Newtonian2DLaw() override;

    /**
     * Operations needed by the base class:
     */

    /**
     * @return Working space dimension constitutive law
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @return Size of the strain vector (in Voigt notation) for the constitutive law
     */
    SizeType GetStrainSize() override;


    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;
    int& GetValue(const Variable<int>& rThisVariable, int& rValue) override;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;
    array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rVariable,array_1d<double, 3 > & rValue) override;
    array_1d<double, 6 > & GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue) override;

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
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}



private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    double mViscosity;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }


}; // Class Newtonian2DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED  defined 

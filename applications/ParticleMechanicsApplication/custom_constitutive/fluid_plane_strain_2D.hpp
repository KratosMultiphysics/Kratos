//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined (KRATOS_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes

#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "includes/checks.h"


namespace Kratos
{
/**
 * The Johnson Cook strain-rate and temperature senstive plastic 3D material law.
 * Requires a strain vector to be provided by the element, which
 * should ideally be objective to enable large displacements.
 * Only suitable for explicit time integration because calculate
 * constitutive tensor is not implemented.
 * Thermal softening may be disabled by setting TAYLOR_QUINNEY_COEFFICIENT=0.0
 * in materials.json
 * Reference:   Lu Ming, Olivier Pantale. An efficient and robust VUMAT implementation of elastoplastic constitutive
 *              laws in Abaqus / Explicit finite element code.Mechanics & Industry, EDP Sciences, 2018, 19 (3),
 *              pp.308.10.1051 / meca / 2018021.hal - 01905414
 */
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) FluidPlaneStrain2DLaw : public HyperElastic3DLaw
{
public:

    /// Type Definitions
    typedef ProcessInfo      ProcessInfoType;
    typedef HyperElastic3DLaw         BaseType;
    typedef std::size_t             SizeType;
    typedef Properties::Pointer            PropertiesPointer;

    /// Counted pointer of FluidPlaneStrain2DLaw
    KRATOS_CLASS_POINTER_DEFINITION(FluidPlaneStrain2DLaw);

    /**
     * Default constructor.
     */
    FluidPlaneStrain2DLaw();

    /**
     * Copy constructor.
     */
    FluidPlaneStrain2DLaw(const FluidPlaneStrain2DLaw & rOther);

    /**
     * Assignment operator.
     */
    FluidPlaneStrain2DLaw& operator=(const FluidPlaneStrain2DLaw & rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~FluidPlaneStrain2DLaw() override;

    /**
     * Operators
     */
     /**
  * Dimension of the law:
  */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
    };

    /**
     * Operations needed by the base class:
     */
    void GetLawFeatures(Features& rFeatures) override;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    void SetValue(const Variable<double>& rVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;

    bool Has(const Variable<double>&rThisVariable) override
    {
        if (rThisVariable == PRESSURE)
        {
            return true;
        }
        else return false;
    }

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    double mPressure = 0.0;
    Vector mStrainOld;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    bool CheckParameters(Parameters& rValues) override;

    void CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo);

private:





    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.save("mPressure", mPressure);
        rSerializer.save("mStrainOld", mStrainOld);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.load("mPressure", mPressure);
        rSerializer.load("mStrainOld", mStrainOld);
    }



}; // Class FluidPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED defined

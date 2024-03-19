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

#if !defined (KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED)
#define  KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_3D_law.hpp"

namespace Kratos
{
	/**
	 * The Johnson Cook strain-rate senstive plastic 2D plane strain material law
	 * derived from the Johnson Cook 3D material law.
	 * Requires a strain vector to be provided by the element, which
	 * should ideally be objective to enable large displacements.
	 * Only suitable for explicit time integration because calculate
	 * constitutive tensor is not implemented.
	 */

	class KRATOS_API(MPM_APPLICATION) JohnsonCookThermalPlastic2DPlaneStrainLaw : public JohnsonCookThermalPlastic3DLaw
	{
	public:

		/// Type Definitions
		typedef ProcessInfo          ProcessInfoType;
		typedef JohnsonCookThermalPlastic3DLaw         BaseType;
		typedef std::size_t             SizeType;
		typedef Properties::Pointer            PropertiesPointer;

		/// Counted pointer of JohnsonCookThermalPlastic2DPlaneStrainLaw
		KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalPlastic2DPlaneStrainLaw);

		/**
		 * Default constructor.
		 */
		JohnsonCookThermalPlastic2DPlaneStrainLaw();

		/**
		 * Copy constructor.
		 */
		JohnsonCookThermalPlastic2DPlaneStrainLaw(const JohnsonCookThermalPlastic2DPlaneStrainLaw& rOther);

		/**
		 * Assignment operator.
		 */
		JohnsonCookThermalPlastic2DPlaneStrainLaw& operator=(const JohnsonCookThermalPlastic2DPlaneStrainLaw& rOther);

		/**
		 * Clone function (has to be implemented by any derived class)
		 * @return a pointer to a new instance of this constitutive law
		 */
		ConstitutiveLaw::Pointer Clone() const override;

		/**
		 * Destructor.
		 */
		~JohnsonCookThermalPlastic2DPlaneStrainLaw() override;

		/// Dimension of the law:
		SizeType WorkingSpaceDimension() override
		{
			return 2;
		};

		/// Voigt tensor size:
		SizeType GetStrainSize() const override
		{
			return 3;
		};

	protected:

		void MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput) override;

		void MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput) override;

	private:

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, JohnsonCookThermalPlastic3DLaw);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, JohnsonCookThermalPlastic3DLaw);
		}
	}; // Class JohnsonCookThermalPlastic2DPlaneStrainLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED defined
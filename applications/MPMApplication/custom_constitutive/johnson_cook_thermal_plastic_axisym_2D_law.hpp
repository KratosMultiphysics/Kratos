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

#if !defined (KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_AXISYM_LAW_H_INCLUDED)
#define  KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_AXISYM_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_strain_2D_law.hpp"

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

	class KRATOS_API(MPM_APPLICATION) JohnsonCookThermalPlastic2DAxisymLaw : public JohnsonCookThermalPlastic2DPlaneStrainLaw
	{
	public:

		/// Type Definitions
		typedef ProcessInfo          ProcessInfoType;
		typedef JohnsonCookThermalPlastic2DPlaneStrainLaw         BaseType;
		typedef std::size_t             SizeType;
		typedef Properties::Pointer            PropertiesPointer;

		/// Counted pointer of JohnsonCookThermalPlastic2DAxisymLaw
		KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalPlastic2DAxisymLaw);

		/**
		 * Default constructor.
		 */
		JohnsonCookThermalPlastic2DAxisymLaw();

		/**
		 * Copy constructor.
		 */
		JohnsonCookThermalPlastic2DAxisymLaw(const JohnsonCookThermalPlastic2DAxisymLaw& rOther);

		/**
		 * Assignment operator.
		 */
		JohnsonCookThermalPlastic2DAxisymLaw& operator=(const JohnsonCookThermalPlastic2DAxisymLaw& rOther);

		/**
		 * Clone function (has to be implemented by any derived class)
		 * @return a pointer to a new instance of this constitutive law
		 */
		ConstitutiveLaw::Pointer Clone() const override;

		/**
		 * Destructor.
		 */
		~JohnsonCookThermalPlastic2DAxisymLaw() override;

		/// Voigt tensor size:
		SizeType GetStrainSize() const override
		{
			return 4;
		};

	protected:

		void MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput) override;

		void MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput) override;

	private:

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, JohnsonCookThermalPlastic2DPlaneStrainLaw);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, JohnsonCookThermalPlastic2DPlaneStrainLaw);
		}
	}; // Class JohnsonCookThermalPlastic2DAxisymLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_2D_AXISYM_LAW_H_INCLUDED defined
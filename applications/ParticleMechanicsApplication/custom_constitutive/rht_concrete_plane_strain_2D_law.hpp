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

#if !defined (KRATOS_RHT_CONCRETE_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_RHT_CONCRETE_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/rht_concrete_3D_law.hpp"

namespace Kratos
{
	/**
	 * The Riedel-Hiermaier-Thoma (RHT) strain-rate senstive plastic plane strain 2D material law.
	 * Requires a strain vector to be provided by the element, which
	 * should ideally be objective to enable large displacements.
	 * Only suitable for explicit time integration because calculate
	 * constitutive tensor is not implemented.
	 * Refer to RHTConcrete3DLaw for references.
	 */
	class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) RHTConcretePlaneStrain2DLaw : public RHTConcrete3DLaw
	{
	public:

		/// Type Definitions
		typedef ProcessInfo          ProcessInfoType;
		typedef RHTConcrete3DLaw         BaseType;
		typedef std::size_t             SizeType;
		typedef Properties::Pointer            PropertiesPointer;

		/// Counted pointer of RHTConcretePlaneStrain2DLaw
		KRATOS_CLASS_POINTER_DEFINITION(RHTConcretePlaneStrain2DLaw);

		/**
		 * Default constructor.
		 */
		RHTConcretePlaneStrain2DLaw();

		/**
		 * Copy constructor.
		 */
		RHTConcretePlaneStrain2DLaw(const RHTConcretePlaneStrain2DLaw& rOther);

		/**
		 * Assignment operator.
		 */
		RHTConcretePlaneStrain2DLaw& operator=(const RHTConcretePlaneStrain2DLaw& rOther);

		/**
		 * Clone function (has to be implemented by any derived class)
		 * @return a pointer to a new instance of this constitutive law
		 */
		ConstitutiveLaw::Pointer Clone() const override;

		/**
		 * Destructor.
		 */
		~RHTConcretePlaneStrain2DLaw() override;

		/// Dimension of the law:
		SizeType WorkingSpaceDimension() override
		{
			return 2;
		};

		/// Voigt tensor size:
		SizeType GetStrainSize() override
		{
			return 3;
		};

	protected:

	private:

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RHTConcrete3DLaw);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RHTConcrete3DLaw);
		}
	}; // Class RHTConcretePlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_RHT_CONCRETE_PLANE_STRAIN_2D_LAW_H_INCLUDED defined
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

#if !defined (KRATOS_RHT_CONCRETE_AXISYM_2D_LAW_H_INCLUDED)
#define  KRATOS_RHT_CONCRETE_AXISYM_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/rht_concrete_plane_strain_2D_law.hpp"

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
	class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) RHTConcreteAxisym2DLaw : public RHTConcretePlaneStrain2DLaw
	{
	public:

		/// Type Definitions
		typedef ProcessInfo          ProcessInfoType;
		typedef RHTConcretePlaneStrain2DLaw         BaseType;
		typedef std::size_t             SizeType;
		typedef Properties::Pointer            PropertiesPointer;

		/// Counted pointer of RHTConcreteAxisym2DLaw
		KRATOS_CLASS_POINTER_DEFINITION(RHTConcreteAxisym2DLaw);

		/**
		 * Default constructor.
		 */
		RHTConcreteAxisym2DLaw();

		/**
		 * Copy constructor.
		 */
		RHTConcreteAxisym2DLaw(const RHTConcreteAxisym2DLaw& rOther);

		/**
		 * Assignment operator.
		 */
		RHTConcreteAxisym2DLaw& operator=(const RHTConcreteAxisym2DLaw& rOther);

		/**
		 * Clone function (has to be implemented by any derived class)
		 * @return a pointer to a new instance of this constitutive law
		 */
		ConstitutiveLaw::Pointer Clone() const override;

		/**
		 * Destructor.
		 */
		~RHTConcreteAxisym2DLaw() override;

		/// Voigt tensor size:
		SizeType GetStrainSize() override
		{
			return 4;
		};

	protected:

	private:

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RHTConcretePlaneStrain2DLaw);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RHTConcretePlaneStrain2DLaw);
		}
	}; // Class RHTConcreteAxisym2DLaw
}  // namespace Kratos.
#endif // KRATOS_RHT_CONCRETE_AXISYM_2D_LAW_H_INCLUDED defined
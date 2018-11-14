//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

#if !defined(KRATOS_ZARATE_LAW_H_INCLUDED)
#define KRATOS_ZARATE_LAW_H_INCLUDED

#include "includes/constitutive_law.h"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"

namespace Kratos
{

class ZarateLaw : public LinearElasticPlaneStrain2DLaw
{
  public:
	KRATOS_CLASS_POINTER_DEFINITION(ZarateLaw);
	/**
		* Default constructor.
		*/
	ZarateLaw();

	/**
		* Clone function (has to be implemented by any derived class)
		* @return a pointer to a new instance of this constitutive law
		*/
	ConstitutiveLaw::Pointer Clone();

	/**
		* Copy constructor.
		*/
	ZarateLaw(const ZarateLaw &rOther);

	/**
		* Assignment operator.
		*/

	//LinearElasticPlaneStrain2DLaw& operator=(const LinearElasticPlaneStrain2DLaw& rOther);

	/**
		* Destructor.
		*/
	virtual ~ZarateLaw();

  private:
	friend class Serializer;

	virtual void save(Serializer &rSerializer)
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearElasticPlaneStrain2DLaw)
	}

	virtual void load(Serializer &rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearElasticPlaneStrain2DLaw)
	}
};
} // namespace Kratos
#endif // KRATOS_ZARATE_LAW_H_INCLUDED

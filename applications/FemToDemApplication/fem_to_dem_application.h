//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


#if !defined(KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

//#include "custom_elements/zaratipito_element.hpp"
#include "custom_constitutive/zarate_law.hpp"
#include "custom_constitutive/fem_dem_elastic_law.hpp"
#include "custom_elements/generic_small_strain_femdem_element.hpp"
#include "custom_elements/generic_large_displacement_femdem_element.hpp"

#include "fem_to_dem_application_variables.h"

namespace Kratos {

class KratosFemToDemApplication : public KratosApplication 

{

public:
	
	KRATOS_CLASS_POINTER_DEFINITION(KratosFemToDemApplication);

	/// Default constructor.
	KratosFemToDemApplication();

	/// Destructor.
	virtual ~KratosFemToDemApplication(){}

	virtual void Register() override; 



	///@}
	///@name Access
	///@{


	///@}
	///@name Inquiry
	///@{


	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const override {
		return "KratosFemToDemApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const override {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const override {
  		KRATOS_WATCH("in my application");
  		KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
    }


protected:

private:
	
	// Elements
	const GenericSmallStrainFemDemElement<2,0> mSmallStrainModifiedMohrCoulombFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,0> mSmallStrainModifiedMohrCoulombFemDemElement3D;
	const GenericSmallStrainFemDemElement<2,1> mSmallStrainRankineFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,1> mSmallStrainRankineFemDemElement3D;
	const GenericSmallStrainFemDemElement<2,2> mSmallStrainSimoJuFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,2> mSmallStrainSimoJuFemDemElement3D;
	const GenericSmallStrainFemDemElement<2,3> mSmallStrainDruckerPragerFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,3> mSmallStrainDruckerPragerFemDemElement3D;
	const GenericSmallStrainFemDemElement<2,4> mSmallStrainVonMisesFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,4> mSmallStrainVonMisesFemDemElement3D;
	const GenericSmallStrainFemDemElement<2,5> mSmallStrainTrescaFemDemElement2D;
	const GenericSmallStrainFemDemElement<3,5> mSmallStrainTrescaFemDemElement3D;

	const GenericLargeDisplacementFemDemElement<2,0> mLargeDisplacementModifiedMohrCoulombFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,0> mLargeDisplacementModifiedMohrCoulombFemDemElement3D;
	const GenericLargeDisplacementFemDemElement<2,1> mLargeDisplacementRankineFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,1> mLargeDisplacementRankineFemDemElement3D;
	const GenericLargeDisplacementFemDemElement<2,2> mLargeDisplacementSimoJuFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,2> mLargeDisplacementSimoJuFemDemElement3D;
	const GenericLargeDisplacementFemDemElement<2,3> mLargeDisplacementDruckerPragerFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,3> mLargeDisplacementDruckerPragerFemDemElement3D;
	const GenericLargeDisplacementFemDemElement<2,4> mLargeDisplacementVonMisesFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,4> mLargeDisplacementVonMisesFemDemElement3D;
	const GenericLargeDisplacementFemDemElement<2,5> mLargeDisplacementTrescaFemDemElement2D;
	const GenericLargeDisplacementFemDemElement<3,5> mLargeDisplacementTrescaFemDemElement3D;

	//elastic laws
   const ZarateLaw mZarateLaw;
   const FemDemElasticLaw mFemDemElasticLaw;
	

	/// Assignment operator.
	KratosFemToDemApplication& operator=(KratosFemToDemApplication const& rOther);

	/// Copy constructor.
	KratosFemToDemApplication(KratosFemToDemApplication const& rOther);

}; // Class KratosFemToDemApplication

}  // namespace Kratos.

#endif // KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED  defined

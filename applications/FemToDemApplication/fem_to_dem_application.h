//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
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
#include "custom_elements/femdem2d_element.hpp"
#include "custom_elements/femdem3d_element.hpp"
#include "custom_elements/romfemdem3d_element.hpp"
#include "custom_elements/femdem3d_large_displacement_element.hpp"
#include "custom_elements/femdem3d_hexahedron_element.hpp"

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

	virtual void Register(); 



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
	virtual std::string Info() const {
		return "KratosFemToDemApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
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
	const FemDem2DElement mFemDem2DElement;
	const FemDem3DElement   mFemDem3DElement;
	const RomFemDem3DElement mRomFemDem3DElement;
	const FemDem3DLargeDisplacementElement mFemDem3DLargeDisplacementElement;
	const FemDem3DHexahedronElement mFemDem3DHexahedronElement;


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

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
// 
//  Project Name:        Kratos       
//  Last Modified by:    $Author:  Miguel Masó Sotomayor
//  Date:                $Date:  april 26 2017
//  Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PURE_DIFFUSION_APPLICATION_H_INCLUDED )
#define  KRATOS_PURE_DIFFUSION_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "pure_diffusion_application_variables.h"
#include "custom_elements/poisson_2d.h"    //including the file for the element
#include "custom_elements/projected_swe.h" //including the file for the second element
#include "includes/condition.h"            //we'll also need conditions for the point heat loads
#include "custom_conditions/point_source.h"         
#include "includes/ublas_interface.h"


namespace Kratos
{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosPureDiffusionApplication : public KratosApplication
	{
	public:

		/// Pointer definition of KratosPureDiffusionApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPureDiffusionApplication);


		/// Default constructor.
		KratosPureDiffusionApplication();


		/// Destructor.
		virtual ~KratosPureDiffusionApplication(){}


		virtual void Register();


		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "KratosPureDiffusionApplication";
		}


		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}


		///// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
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

		const Poisson2D mPoisson2D;        // element
		const ProjectedSWE mProjectedSWE;  // element
		const PointSource mPointSource;    // condition

		/// Assignment operator.
		KratosPureDiffusionApplication& operator=(KratosPureDiffusionApplication const& rOther);

		/// Copy constructor.
		KratosPureDiffusionApplication(KratosPureDiffusionApplication const& rOther);



	}; // Class KratosPureDiffusionApplication 



}  // namespace Kratos.

#endif // KRATOS_PURE_DIFFUSION_APPLICATION_H_INCLUDED  defined 



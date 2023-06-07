//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//                   Suneth Warnakulasuriya

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/kratos_application.h"

/* ELEMENTS */
#include "custom_elements/helmholtz_surf_shape_element.h"
#include "custom_elements/helmholtz_surf_thickness_element.h"
#include "custom_elements/helmholtz_bulk_shape_element.h"
#include "custom_elements/helmholtz_bulk_element.h"
#include "custom_elements/helmholtz_element.h"

/* ELEMENT DATA CONTAINERS*/
#include "custom_elements/data_containers/helmholtz_surface_data_container.h"
#include "custom_elements/data_containers/helmholtz_solid_data_container.h"
#include "custom_elements/data_containers/helmholtz_solid_shape_data_container.h"

/* ADJOINT ELEMENTS */
#include "custom_elements/adjoint_small_displacement_element.h"

/* CONDITIONS */
#include "custom_conditions/helmholtz_surf_shape_condition.h"
#include "custom_conditions/helmholtz_surface_shape_condition.h"

/* CONSTITUTIVE LAWS */
#include "custom_constitutive/helmholtz_jacobian_stiffened_3d.h"

namespace Kratos
{
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KRATOS_API(OPTIMIZATION_APPLICATION) KratosOptimizationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of KratosOptimizationApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosOptimizationApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosOptimizationApplication();

		/// Destructor.
		~KratosOptimizationApplication() override {}

		KratosOptimizationApplication& operator=(KratosOptimizationApplication const& rOther) = delete;

		/// Copy constructor.
		KratosOptimizationApplication(KratosOptimizationApplication const& rOther) = delete;

		///@}
		///@name Operations
		///@{

	    void Register() override;

		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		std::string Info() const override
		{
			return "KratosOptimizationApplication";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream& rOStream) const override
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
		void PrintData(std::ostream& rOStream) const override
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

		///@}

	private:
		///@name Member Variables
		///@{

		/* ELEMENTS */

		const HelmholtzSurfShapeElement mHelmholtzSurfShape3D3N;
		const HelmholtzSurfThicknessElement mHelmholtzSurfThickness3D3N;
		const HelmholtzBulkShapeElement mHelmholtzBulkShape3D4N;
		const HelmholtzBulkElement mHelmholtzBulkTopology3D4N;

		/* ADJ ELEMENTS */
		const AdjointSmallDisplacementElement mAdjointSmallDisplacementElement3D4N;

		// Helmholtz elements
		const HelmholtzElement<HelmholtzSurfaceDataContainer<3, 3, 1>> mHelmholtzSurfaceElement3D3N;
		const HelmholtzElement<HelmholtzSurfaceDataContainer<3, 4, 1>> mHelmholtzSurfaceElement3D4N;
		const HelmholtzElement<HelmholtzSurfaceDataContainer<3, 3, 3>> mHelmholtzVectorSurfaceElement3D3N;
		const HelmholtzElement<HelmholtzSurfaceDataContainer<3, 4, 3>> mHelmholtzVectorSurfaceElement3D4N;

		const HelmholtzElement<HelmholtzSolidDataContainer<3, 4, 1>> mHelmholtzSolidElement3D4N;
		const HelmholtzElement<HelmholtzSolidDataContainer<3, 8, 1>> mHelmholtzSolidElement3D8N;
		const HelmholtzElement<HelmholtzSolidDataContainer<3, 4, 3>> mHelmholtzVectorSolidElement3D4N;
		const HelmholtzElement<HelmholtzSolidDataContainer<3, 8, 3>> mHelmholtzVectorSolidElement3D8N;

		const HelmholtzElement<HelmholtzSolidShapeDataContainer<3, 4>> mHelmholtzSolidShapeElement3D4N;
		const HelmholtzElement<HelmholtzSolidShapeDataContainer<3, 8>> mHelmholtzSolidShapeElement3D8N;

		/* CONDITIONS*/
		// Surface conditions
		const HelmholtzSurfShapeCondition mHelmholtzSurfShapeCondition3D3N;
		const HelmholtzSurfaceShapeCondition mHelmholtzSurfaceShapeCondition3D3N;
		const HelmholtzSurfaceShapeCondition mHelmholtzSurfaceShapeCondition3D4N;

		/* CONSTITUTIVE LAWS */
		const HelmholtzJacobianStiffened3D mHelmholtzJacobianStiffened3D;

		///@}

	}; // Class KratosOptimizationApplication

	///@}

}  // namespace Kratos.



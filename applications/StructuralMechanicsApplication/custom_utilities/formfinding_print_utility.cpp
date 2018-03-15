//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher 
//
	        

// System includes


// External includes 


// Project includes
#include "includes/checks.h"
#include "shell_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
    namespace ShellUtilities
	{
        typedef Properties PropertiesType; // TODO remove this?

		double dN_seren_dxi(const int actualNodeNumber,const double xi,
			const double eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1.0)*(-0.25*xi + 0.25) -
					0.25*(-eta + 1.0)*(-eta - xi - 1.0);
				break;
			case 2:
				returnValue = (-eta + 1.0)*(0.25*xi + 0.25) +
					0.25*(-eta + 1.0)*(-eta + xi - 1.0);
				break;
			case 3:
				returnValue = (eta + 1.0)*(0.25*xi + 0.25) +
					0.25*(eta + 1.0)*(eta + xi - 1.0);
				break;
			case 4:
				returnValue = -(eta + 1.0)*(-0.25*xi + 0.25) -
					0.25*(eta + 1.0)*(eta - xi - 1.0);
				break;
			case 5:
				returnValue = -1.0*xi*(-eta + 1.0);
				break;
			case 6:
				returnValue = -0.5*eta*eta + 0.5;
				break;
			case 7:
				returnValue = -1.0*xi*(eta + 1.0);
				break;
			case 8:
				returnValue = 0.5*eta*eta - 0.5;
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return returnValue;
		}

		double dN_seren_deta(const int actualNodeNumber,const double xi,
			const double eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1.0)*(-0.25*xi + 0.25) -
					(-0.25*xi + 0.25)*(-eta - xi - 1.0);
				break;
			case 2:
				returnValue = -(-eta + 1.0)*(0.25*xi + 0.25) -
					(0.25*xi + 0.25)*(-eta + xi - 1.0);
				break;
			case 3:
				returnValue = (eta + 1.0)*(0.25*xi + 0.25) +
					(0.25*xi + 0.25)*(eta + xi - 1.0);
				break;
			case 4:
				returnValue = (eta + 1.0)*(-0.25*xi + 0.25) +
					(-0.25*xi + 0.25)*(eta - xi - 1.0);
				break;
			case 5:
				returnValue = 0.5*xi*xi - 0.5;
				break;
			case 6:
				returnValue = -1.0*eta*(xi + 1.0);
				break;
			case 7:
				returnValue = -0.5*xi*xi + 0.5;
				break;
			case 8:
				returnValue = -1.0*eta*(-xi + 1.0);
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return returnValue;
		}
		
		void InterpToStandardGaussPoints(double& v1, double& v2,
			double& v3)
		{
			double vg1 = v1;
			double vg2 = v2;
			double vg3 = v3;
#ifdef OPT_AVERAGE_RESULTS
			v1 = (vg1 + vg2 + vg3) / 3.0;
			v2 = (vg1 + vg2 + vg3) / 3.0;
			v3 = (vg1 + vg2 + vg3) / 3.0;
#else
			v1 = (2.0*vg1) / 3.0 - vg2 / 3.0 + (2.0*vg3) / 3.0;
			v2 = (2.0*vg1) / 3.0 + (2.0*vg2) / 3.0 - vg3 / 3.0;
			v3 = (2.0*vg2) / 3.0 - vg1 / 3.0 + (2.0*vg3) / 3.0;
#endif // OPT_AVERAGE_RESULTS
		}

		void InterpToStandardGaussPoints(std::vector< double >& v)
		{
			if (v.size() != 3) return;
			InterpToStandardGaussPoints(v[0], v[1], v[2]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			3> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 3; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,
			6> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 6; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Vector >& v)
		{
			if (v.size() != 3) return;
			size_t ncomp = v[0].size();
			for (int i = 1; i < 3; i++)
				if (v[i].size() != ncomp)
					return;
			for (size_t i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Matrix >& v)
		{
			if (v.size() != 3) return;
			size_t nrows = v[0].size1();
			size_t ncols = v[0].size2();
			for (int i = 1; i < 3; i++)
				if (v[i].size1() != nrows || v[i].size2() != ncols)
					return;
			for (size_t i = 0; i < nrows; i++)
				for (size_t j = 0; j < ncols; j++)
					InterpToStandardGaussPoints
					(v[0](i, j), v[1](i, j), v[2](i, j));
		}



        void CheckVariables()
        {
            KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
            KRATOS_CHECK_VARIABLE_KEY(ROTATION);
            KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
            KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
            KRATOS_CHECK_VARIABLE_KEY(DENSITY);
            KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
            KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
            KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);        
        }

        void CheckDofs(GeometryType& rGeom)
        {
            // verify that the dofs exist
            for (unsigned int i = 0; i < rGeom.size(); i++)
            {
                auto& r_node = rGeom[i];
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);

                KRATOS_CHECK_DOF_IN_NODE(ROTATION_X, r_node);
                KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y, r_node);
                KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z, r_node);

                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);

                if (r_node.GetBufferSize() < 2)
                    KRATOS_ERROR << "This Element needs at least a buffer size = 2" << std::endl;
            }
        }

        void CheckProperties(const Element* pTheElement, const ProcessInfo& rCurrentProcessInfo, const bool IsThickShell)
        {
            // check properties
            if(pTheElement->pGetProperties() == nullptr)
                KRATOS_ERROR << "Properties not provided for element " << pTheElement->Id() << std::endl;

            const PropertiesType & props = pTheElement->GetProperties();

            const GeometryType& geom = pTheElement->GetGeometry(); // TODO check if this can be const

            if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
            {
                const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
                if(section == nullptr)
                    KRATOS_ERROR << "SHELL_CROSS_SECTION not provided for element " << pTheElement->Id() << std::endl;
        
                section->Check(props, geom, rCurrentProcessInfo);
            }
            else if (props.Has(SHELL_ORTHOTROPIC_LAYERS))
            {
                CheckSpecificProperties(pTheElement, props, IsThickShell);
        
                // perform detailed orthotropic check later in shell_cross_section
            }
            else // ... allow the automatic creation of a homogeneous section from a material and a thickness
            {
                CheckSpecificProperties(pTheElement, props, IsThickShell);
        
                ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
                dummySection->BeginStack();
                dummySection->AddPly(props[THICKNESS], 0.0, 5, pTheElement->pGetProperties());
                dummySection->EndStack();
                dummySection->SetSectionBehavior(ShellCrossSection::Thick);
                dummySection->Check(props, geom, rCurrentProcessInfo);
            }

        }

        void CheckSpecificProperties(const Element* pTheElement, const PropertiesType & rProps, const bool IsThickShell)
        {
            if (!rProps.Has(CONSTITUTIVE_LAW))
                KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << pTheElement->Id() << std::endl;
            const ConstitutiveLaw::Pointer& claw = rProps[CONSTITUTIVE_LAW];
            if (claw == nullptr)
                KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << pTheElement->Id() << std::endl;

            if(!rProps.Has(THICKNESS))
                KRATOS_ERROR << "THICKNESS not provided for element " << pTheElement->Id() << std::endl;
            if(rProps[THICKNESS] <= 0.0)
                KRATOS_ERROR << "wrong THICKNESS value provided for element " << pTheElement->Id() << std::endl;
                
            if(!rProps.Has(DENSITY))
                KRATOS_ERROR << "DENSITY not provided for element " << pTheElement->Id() << std::endl;
            if(rProps[DENSITY] < 0.0)
                KRATOS_ERROR << "wrong DENSITY value provided for element " << pTheElement->Id() << std::endl;
            
            if(IsThickShell)
            {
                // Check constitutive law has been verified with Stenberg stabilization
                // applicable for 5-parameter shells only.
				bool stenberg_stabilization_suitable = false;
				claw->GetValue(STENBERG_SHEAR_STABILIZATION_SUITABLE, stenberg_stabilization_suitable);
                if (!stenberg_stabilization_suitable)
                {
                    std::cout << "\nWARNING:\nThe current constitutive law has not been checked with Stenberg shear stabilization."
                        << "\nPlease check results carefully."
                        << std::endl;
                }
            }
        }  
    } // namespace ShellUtilities
}  // namespace Kratos.



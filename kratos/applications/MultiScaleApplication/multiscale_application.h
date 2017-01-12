//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_MULTISCALE_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTISCALE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

#include "includes/kratos_application.h"

#include "includes/variables.h"
#include "includes/ublas_interface.h"

#include "custom_elements/small_displacement_interface_element.hpp"
#include "custom_elements/shell_thick_interface_element.hpp"
#include "custom_elements/opt_triangle_element.hpp"
#include "custom_elements/eas_quad_element_v2.hpp"
#include "custom_elements/q4ristab_element.hpp"
#include "custom_elements/convdiff_interface_element.hpp"
#include "custom_elements/convdiff_element.hpp"

#include "custom_elements/ebst_element_2d3n.h"
#include "custom_elements/agq4_element.hpp"
#include "custom_elements/small_displacement_elastic_link_element.hpp"
#include "custom_conditions/periodic_condition_lm_2D2N.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/// MultiScale Application for KRATOS.
/**
 * This application is a Demo
 */
class KratosMultiScaleApplication : public KratosApplication
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMultiScaleApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultiScaleApplication);
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMultiScaleApplication();

    /// Destructor.
    virtual ~KratosMultiScaleApplication() {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * Registers the structural application in the KRATOS kernel
     */
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
    virtual std::string Info() const
    {
        return "KratosMultiScaleApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in KratosMultiScaleApplication application");
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

    ///@name Friends
    ///@{
    ///@}

public:

	///@name Member Variables
	///@{

	struct KeyComparor
	{
		bool operator()(const vector<double>& lhs, const vector<double>& rhs) const
		{
			if (lhs.size() != rhs.size())
				return false;

			for (unsigned int i = 0; i<lhs.size(); i++)
			{
				if (lhs[i] != rhs[i]) return false;
			}

			return true;
		}
	};

	struct KeyHasher
	{
		std::size_t operator()(const vector<int>& k) const
		{
			return boost::hash_range(k.begin(), k.end());
		}
	};

	/*Unordered maps are associative containers that store elements
	* formed by the combination of a key value and a mapped value,
	* and which allows for fast retrieval of individual elements based on their keys.
	*
	* In an unordered_map, the key value is generally used to uniquely identify the element,
	* while the mapped value is an object with the content associated to this key.
	* Types of key and mapped value may differ.
	*
	* Internally, the elements in the unordered_map are not sorted in any particular order
	* with respect to either their key or mapped values, but organized into buckets depending
	* on their hash values to allow for fast access to individual elements directly by their
	* key values (with a constant average time complexity on average).
	*
	* unordered_map containers are faster than map containers to access individual elements
	* by their key, although they are generally less efficient for range iteration through
	* a subset of their elements.
	*
	* Unordered maps implement the direct access operator (operator[]) which allows for
	* direct access of the mapped value using its key value as argument.
	* Ref.: http://www.cplusplus.com/reference/unordered_map/unordered_map/
	*/
	typedef boost::unordered_map<vector<int>, Vector, KeyHasher, KeyComparor > UnorderedVectorVectorMap;

	typedef boost::unordered_map<vector<int>, Matrix, KeyHasher, KeyComparor> UnorderedVectorMatrixMap;

	///@}

private:

    ///@name Member Variables
    ///@{

	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement2D4N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D6N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D8N;

	const ShellThickInterfaceElement mShellThickInterfaceElement3D4N;
	
	const OptTriangleElement mOptTriangleElement2D3N;
	const EASQuadElementV2 mEASQuadElementV22D4N;
	const Q4RIStabElement mQ4RIStabElement2D4N;

	const ConvDiffInterfaceElement mConvDiffInterfaceElement2D4N;
	const ConvDiffInterfaceElement mConvDiffInterfaceElement3D6N;
	const ConvDiffInterfaceElement mConvDiffInterfaceElement3D8N;

	const ConvDiffElement mConvDiffElement2D3N;
	const ConvDiffElement mConvDiffElement2D4N;
	const ConvDiffElement mConvDiffElement3D4N;
	const ConvDiffElement mConvDiffElement3D8N;
	
	const EBSTElement2D3N mEBSTElement2D3N;
	const AGQ4Element mAGQ4Element2D4N;
	const SmallDisplacementElasticLinkElement mSmallDisplacementElasticLinkElement3D2;

	const PeriodicConditionLM2D2N mPeriodicConditionLM2D2N;

    ///@}

    ///@name Private Operators
    ///@{
    ///@}

    ///@name Private Operations
    ///@{
    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosMultiScaleApplication& operator=(KratosMultiScaleApplication const& rOther);

    /// Copy constructor.
    KratosMultiScaleApplication(KratosMultiScaleApplication const& rOther);

    ///@}

};

///@}

}
#endif // KRATOS_MULTISCALE_APPLICATION_H_INCLUDED  defined 

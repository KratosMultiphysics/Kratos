// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//

#ifndef THREAD_SPECIFICATION_UTILITY_H
#define THREAD_SPECIFICATION_UTILITY_H


// ==============================================================================

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

/// Short class definition.
/** Detail class definition.

 */

class ThreadSpecificationUtility
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of ThreadSpecificationUtility
	KRATOS_CLASS_POINTER_DEFINITION(ThreadSpecificationUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	ThreadSpecificationUtility()
	{
	}

	/// Destructor.
	virtual ~ThreadSpecificationUtility()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	// ==============================================================================
	void SetNumberOfThreads(int num_threads)
	{
        omp_set_num_threads(num_threads);
	}

	// --------------------------------------------------------------------------
	int GetMaxNumberOfThreads()
	{
		return omp_get_max_threads();
	}

	// ==============================================================================

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
	std::string Info() const
	{
		return "ThreadSpecificationUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "ThreadSpecificationUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}
	///@name Friends
	///@{

	///@}

protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	// ==============================================================================

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

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
	//      ThreadSpecificationUtility& operator=(ThreadSpecificationUtility const& rOther);

	/// Copy constructor.
	//      ThreadSpecificationUtility(ThreadSpecificationUtility const& rOther);

	///@}

}; // Class ThreadSpecificationUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // THREAD_SPECIFICATION_UTILITY_H

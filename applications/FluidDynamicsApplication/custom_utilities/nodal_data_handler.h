//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_NODAL_DATA_HANDLER_H_INCLUDED)
#define KRATOS_NODAL_DATA_HANDLER_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/data_handler.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template <class TDataType, unsigned int TNumNodes, class TStorageType>
class NodalDataHandler : public DataHandler<TDataType, TStorageType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalDataHandler
    KRATOS_CLASS_POINTER_DEFINITION(NodalDataHandler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NodalDataHandler(const Variable<TDataType>& rVariable);

    /// Destructor.
    ~NodalDataHandler() override;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override;

	TDataType Interpolate(const boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement) override;

    int Check(const Element& rElement) override;

    ///@}
    ///@name Access
    ///@{

    TStorageType& Get() override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NodalDataHandler for " << this->mrVariable.Name() << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NodalDataHandler for " << this->mrVariable.Name() << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    TStorageType mValues;

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodalDataHandler& operator=(NodalDataHandler const& rOther);

    /// Copy constructor.
    NodalDataHandler(NodalDataHandler const& rOther);

    ///@}

}; // Class NodalDataHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TDataType, unsigned int TNumNodes, class TStorageType>
inline std::istream& operator>>(std::istream& rIStream, NodalDataHandler<TDataType,TNumNodes,TStorageType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TDataType, unsigned int TNumNodes, class TStorageType>
inline std::ostream& operator<<(std::ostream& rOStream, const NodalDataHandler<TDataType,TNumNodes,TStorageType>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

///@} addtogroup block

namespace Internals {
template <size_t TNumNodes>
void SpecializedInitialization(array_1d<double, TNumNodes>& rValues,
                               const Variable<double>& rVariable,
                               const Element& rElement,
                               const ProcessInfo& rProcessInfo)
{
    rValues = array_1d<double, TNumNodes>(TNumNodes, 0.0);
    const Geometry<Node<3>>& r_geometry = rElement.GetGeometry();
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rValues[i] = r_geometry[i].FastGetSolutionStepValue(rVariable);
    }
}

template <size_t TNumNodes, size_t TDim>
void SpecializedInitialization(boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rValues,
                               const Variable<array_1d<double,3>>& rVariable,
                               const Element& rElement,
                               const ProcessInfo& rProcessInfo)
{
    const Geometry<Node<3>>& r_geometry = rElement.GetGeometry();
    for (size_t i = 0; i < TNumNodes; i++)
    {
        const array_1d<double,3>& r_nodal_values = r_geometry[i].FastGetSolutionStepValue(rVariable);
        for (size_t j = 0; j < rValues.size2(); j++)
        {
            rValues(i,j) = r_nodal_values[j];
        }
    }
}

template < size_t TNumNodes >
double SpecializedInterpolation(array_1d<double,TNumNodes>& rValues, const boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement)
{
    double result = 0.0;
	for (size_t i = 0; i < TNumNodes; i++)
	{
		result += rN[i] * rValues[i];
	}

	return result;
}

template < size_t TNumNodes, size_t TDim >
array_1d<double,3> SpecializedInterpolation(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim>& rValues, const boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement)
{
    array_1d<double, 3> result(3, 0.0);

    for (size_t i = 0; i < TNumNodes; i++) {
        for (size_t j = 0; j < TDim; j++) {
            result[j] += rN[i] * rValues(i,j);
        }
    }

	return result;
}

}

} // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined

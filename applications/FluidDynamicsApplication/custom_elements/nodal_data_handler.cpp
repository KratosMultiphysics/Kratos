#include "nodal_data_handler.h"

#include "containers/array_1d.h"
#include "boost/numeric/ublas/matrix.hpp"

namespace Kratos
{
template class NodalDataHandler<double, 3, array_1d<double, 3>>;
template class NodalDataHandler<array_1d<double, 3>, 3, boost::numeric::ublas::bounded_matrix<double, 2, 3>>;

template <class TDataType, unsigned int TNumNodes, class TStorageType>
NodalDataHandler<TDataType, TNumNodes, TStorageType>::NodalDataHandler(const Variable<TDataType>& rVariable):
    DataHandler<TDataType,TStorageType>(rVariable)
{
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
NodalDataHandler<TDataType, TNumNodes, TStorageType>::~NodalDataHandler()
{
}

// Variable<double> version ///////////////////////////////////////////////////

template <>
void NodalDataHandler<double, 3, array_1d<double, 3>>::Initialize(const Element& rElement,
                                                                  const ProcessInfo& rProcessInfo)
{
    mValues = array_1d<double, 3>(3, 0.0); // NOTE: The 3 here is TNumNodes!!!
    const Geometry<Node<3>> r_geometry = rElement.GetGeometry();
    for (unsigned int i = 0; i < 3; i++) // NOTE: TNumNodes
    {
        mValues[i] = r_geometry[i].FastGetSolutionStepValue(this->mrVariable);
    }
}

template<>
double NodalDataHandler<double, 3, array_1d<double, 3>>::Interpolate(boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement)
{
	double result = 0.0;
	for (unsigned int i = 0; i < 3; i++) // NOTE: The 3 is TNumNodes
	{
		result += rN[i] * this->mValues[i];
	}

	return result;
}


template <>
array_1d<double, 3>& NodalDataHandler<double, 3, array_1d<double, 3>>::Get()
{
    return mValues;
}

// Variable< array_1d<double,3> > version /////////////////////////////////////

template <>
void NodalDataHandler< array_1d<double, 3>, 3, boost::numeric::ublas::bounded_matrix<double, 2, 3>>::Initialize(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry<Node<3>> r_geometry = rElement.GetGeometry();
    for (unsigned int i = 0; i < 3; i++) // NOTE: TNumNodes
    {
        const array_1d<double,3>& r_nodal_values = r_geometry[i].FastGetSolutionStepValue(this->mrVariable);
        for (unsigned int j = 0; j < mValues.size1(); j++)
        {
            mValues(i,j) = r_nodal_values[j];
        }
    }
}

template <>
array_1d<double, 3> NodalDataHandler< array_1d<double, 3>, 3, boost::numeric::ublas::bounded_matrix<double, 2, 3>>::Interpolate(boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement)
{
	array_1d<double, 3> result(3, 0.0);
	noalias(result) = boost::numeric::ublas::prod(rN, mValues);

	return result;
}

template <>
boost::numeric::ublas::bounded_matrix<double, 2, 3>& NodalDataHandler<array_1d<double, 3>, 3, boost::numeric::ublas::bounded_matrix<double, 2, 3>>::Get()
{
    return mValues;
}

}
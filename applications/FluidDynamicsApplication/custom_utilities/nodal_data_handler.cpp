#include "nodal_data_handler.h"

#include "includes/checks.h"
#include "containers/array_1d.h"
#include "boost/numeric/ublas/matrix.hpp"

namespace Kratos
{

// Variable<double> version
template class NodalDataHandler<double, 3, array_1d<double, 3>>;
template class NodalDataHandler<double, 4, array_1d<double, 4>>;

// Variable< array_1d<double,3> > version (2D case)
template class NodalDataHandler<array_1d<double, 3>, 3, boost::numeric::ublas::bounded_matrix<double, 3, 2>>;

// Variable< array_1d<double,3> > version (3D case)
template class NodalDataHandler<array_1d<double, 3>, 4, boost::numeric::ublas::bounded_matrix<double, 4, 3>>;

template <class TDataType, unsigned int TNumNodes, class TStorageType>
NodalDataHandler<TDataType, TNumNodes, TStorageType>::NodalDataHandler(const Variable<TDataType>& rVariable):
    DataHandler<TDataType,TStorageType>(rVariable)
{
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
NodalDataHandler<TDataType, TNumNodes, TStorageType>::~NodalDataHandler()
{
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
void NodalDataHandler<TDataType, TNumNodes, TStorageType>::Initialize(const Element& rElement,
                                                                      const ProcessInfo& rProcessInfo)
{
    Internals::SpecializedInitialization(mValues,this->mrVariable,rElement,rProcessInfo);
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
TDataType NodalDataHandler<TDataType, TNumNodes, TStorageType>::Interpolate(
    const boost::numeric::ublas::matrix_row<Matrix>& rN, Element* pElement)
{
    return Internals::SpecializedInterpolation(mValues, rN, pElement);
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
int NodalDataHandler<TDataType,TNumNodes,TStorageType >::Check(const Element& rElement) {

    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    
    KRATOS_CHECK_VARIABLE_KEY(this->mrVariable);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(this->mrVariable,r_geometry[i]);
    }

    return 0;
}

template <class TDataType, unsigned int TNumNodes, class TStorageType>
TStorageType& NodalDataHandler<TDataType,TNumNodes,TStorageType >::Data() {
    return mValues;
}

}
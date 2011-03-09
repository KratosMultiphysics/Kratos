#include "bingham_vms.h"

namespace Kratos
{
    ///@name Specialized implementation of VMS for functions that depend on TDim
    ///@{

    template <>
    void VMS<2,3>::save(Serializer& rSerializer) const
    {
        typedef VMS<2,3> BaseClassType;
        rSerializer.save("Name","BinghamVMS2D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
    }

    template <>
    void VMS<3,4>::save(Serializer& rSerializer) const
    {
        typedef VMS<3,4> BaseClassType;
        rSerializer.save("Name","BinghamVMS3D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType );
    }

    ///@} // Specialized implementations
}

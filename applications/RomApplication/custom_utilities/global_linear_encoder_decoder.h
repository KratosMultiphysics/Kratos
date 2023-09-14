//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

#if !defined(GLOBAL_LINEAR_ENCODER_DECODER)
#define GLOBAL_LINEAR_ENCODER_DECODER

#include "custom_utilities/base_encoder_decoder.h"

namespace Kratos{

class GlobalLinearEncoderDecoder: public BaseEncoderDecoder{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(GlobalLinearEncoderDecoder);

        GlobalLinearEncoderDecoder(Parameters ThisParameters): BaseEncoderDecoder(ThisParameters)
        {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "some_default_parameter" : "to_define_later"
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        }

        ~GlobalLinearEncoderDecoder()= default;

        virtual void Encode() override{
        }

        virtual void Decode() override{
        }

        virtual void GetEncoderDerivative() override{
        }

        virtual void GetDecoderDerivative() override{
        }

        };

}

#endif // GLOBAL_LINEAR_ENCODER_DECODER
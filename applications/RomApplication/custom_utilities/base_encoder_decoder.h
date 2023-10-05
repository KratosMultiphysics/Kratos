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

#if !defined(BASE_ENCODER_DECODER)
#define BASE_ENCODER_DECODER


namespace Kratos{


class BaseEncoderDecoder{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(BaseEncoderDecoder);

        BaseEncoderDecoder(){}

        // BaseEncoderDecoder(Parameters ThisParameters)
        // {
        // // Validate default parameters
        // Parameters default_parameters = Parameters(R"(
        // {
        //     "some_default_parameter" : "to_define_later"
        // })" );

        // ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // }

        ~BaseEncoderDecoder()= default;

        virtual void Encode(){};

        virtual void Decode(){};

        virtual void GetEncoderDerivative(){};

        virtual Matrix GetDecoderDerivative(int node_id){};

    };


class LinearEncoderDecoder: public BaseEncoderDecoder{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(LinearEncoderDecoder);

        LinearEncoderDecoder(){}

        // LinearEncoderDecoder(Parameters ThisParameters): BaseEncoderDecoder(ThisParameters)
        // {
        // // Validate default parameters
        // Parameters default_parameters = Parameters(R"(
        // {
        //     "some_default_parameter" : "to_define_later"
        // })" );

        // ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // }

        ~LinearEncoderDecoder()= default;

        virtual void Encode() override {
        }

        virtual void Decode() override {
        }

        virtual void GetEncoderDerivative() override {
        }

        virtual Matrix GetDecoderDerivative(int node_id) override {
            return mMapPhi[node_id];
        }

        void SetNodalBasis(int node_id, Matrix &nodal_basis){
            mMapPhi[node_id] = nodal_basis;
        }

        protected:
        std::unordered_map<int, Matrix> mMapPhi;

    };


}


#endif // BASE_ENCODER_DECODER


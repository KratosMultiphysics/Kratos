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


    // class RomBasis
    // {
    //     public:

    //     KRATOS_CLASS_POINTER_DEFINITION(RomBasis);

    //     ~RomBasis()= default;

    //     void SetNodalBasis(int node_id, Matrix &nodal_basis){
    //         //mMapPhi[node_id] = the_basis;
    //         std::shared_ptr<Matrix> pnodal_basis{new Matrix{nodal_basis}};
    //         mMapPhi[node_id] = pnodal_basis;
    //         }

    //     std::shared_ptr<Matrix> GetNodalBasis(int node_id){
    //         return mMapPhi[node_id];
    //     }

    //     protected:
    //     //std::unordered_map<int, Matrix> mMapPhi;
    //     std::unordered_map<int, std::shared_ptr<Matrix>> mMapPhi;

    // };

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

        virtual void Encode(){
        }

        virtual void Decode(){
        }

        virtual void GetEncoderDerivative(){
        }

        virtual std::shared_ptr<Matrix> GetDecoderDerivative(int node_id){
        }
    };


class GlobalLinearEncoderDecoder: public BaseEncoderDecoder{
    public:

        KRATOS_CLASS_POINTER_DEFINITION(GlobalLinearEncoderDecoder);

        GlobalLinearEncoderDecoder(){}

        // GlobalLinearEncoderDecoder(Parameters ThisParameters): BaseEncoderDecoder(ThisParameters)
        // {
        // // Validate default parameters
        // Parameters default_parameters = Parameters(R"(
        // {
        //     "some_default_parameter" : "to_define_later"
        // })" );

        // ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // }

        ~GlobalLinearEncoderDecoder()= default;

        virtual void Encode() {
        }

        virtual void Decode() {
        }

        virtual void GetEncoderDerivative() {
        }

        virtual std::shared_ptr<Matrix> GetDecoderDerivative(int node_id) {
            return mMapPhi[node_id];
        }

        void SetNodalBasis(int node_id, Matrix &nodal_basis){
            //mMapPhi[node_id] = the_basis;
            std::shared_ptr<Matrix> pnodal_basis{new Matrix{nodal_basis}};
            mMapPhi[node_id] = pnodal_basis;
            }

        // std::shared_ptr<Matrix> GetNodalBasis(int node_id){
        //     return mMapPhi[node_id];
        // }

        protected:
        //std::unordered_map<int, Matrix> mMapPhi;
        std::unordered_map<int, std::shared_ptr<Matrix>> mMapPhi;


    };




}


#endif // BASE_ENCODER_DECODER
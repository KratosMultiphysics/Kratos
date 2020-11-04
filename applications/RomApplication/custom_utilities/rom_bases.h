//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    RAUL BRAVO
//

#if !defined( ROM_BASES_H_INCLUDED )
#define  ROM_BASES_H_INCLUDED

/* Project includes */


/* Application includes */

namespace Kratos
{

    class RomBasis
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomBasis);

        // RomBasis(){
        //     KRATOS_WATCH('Basis CREATED');
        // }

        // RomBasis(const RomBasis &source){
        //     std::cout<<'Copy constructor called'<<std::endl;
        // }

        // RomBasis& operator=(const RomBasis &source)
        // {
        //     KRATOS_WATCH('Assignment operator called');
        //     return *this;
        // }

        ~RomBasis()= default;

        void SetNodalBasis(int &node_id, Matrix &the_basis){
            mMapPhi[node_id] = the_basis;
            }

        Matrix GetNodalBasis(int node_id){
            return mMapPhi[node_id];
        }

        protected:
        std::unordered_map<int, Matrix> mMapPhi;


    };

    class RomBases
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomBases);

        // RomBases(const RomBases &source){
        //     std::cout<<'Copy constructor called'<<std::endl;
        // }

        // RomBases& operator=(const RomBases &source)
        // {
        //     KRATOS_WATCH('Assignment operator called');
        //     return *this;
        // }

        // RomBases(){
        //     KRATOS_WATCH('Bases CREATED')
        // }

        ~RomBases()= default;


    void AddBasis(int &basis_id, RomBasis &new_basis){
        // bases should be added in the correct order
        //mBases.push_back(new_basis);
        mMapBases[basis_id] = new_basis;
        }

    RomBasis GetBasis(int k){
        return mMapBases[k];
        //return mBases.at(k);
    }

        protected:
        std::unordered_map<int, RomBasis> mMapBases;
        //std::vector<RomBasis> mBases;

    };






} // namespace Kratos



#endif // ROM_BASES_H_INCLUDED  defined

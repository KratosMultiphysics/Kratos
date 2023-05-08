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

        ~RomBasis()= default;

        void SetNodalBasis(const int &node_id, const Matrix &nodal_basis){
            //mMapPhi[node_id] = the_basis;
            std::shared_ptr<Matrix> pnodal_basis{new Matrix{nodal_basis}};
            mMapPhi[node_id] = pnodal_basis;
            }

        std::shared_ptr<Matrix> GetNodalBasis(const int &node_id){
            return mMapPhi[node_id];
        }

        protected:
        //std::unordered_map<int, Matrix> mMapPhi;
        std::unordered_map<int, std::shared_ptr<Matrix>> mMapPhi;

    };

    class RomBases
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomBases);

        ~RomBases()= default;


    void AddBasis(const int &basis_id, const RomBasis &basis){
        // bases should be added in the correct order
        //mBases.push_back(new_basis);
        std::shared_ptr<RomBasis> pbasis{new RomBasis{basis}};
        mMapBases[basis_id] = pbasis;
        }

    std::shared_ptr<RomBasis> GetBasis(const int &k){
        return mMapBases[k];
        //return mBases.at(k);
    }

        protected:
        std::unordered_map<int, std::shared_ptr<RomBasis>> mMapBases;
        //std::vector<RomBasis> mBases;

    };



} // namespace Kratos



#endif // ROM_BASES_H_INCLUDED  defined

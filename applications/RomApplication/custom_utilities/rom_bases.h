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




    class MultipleBasesManager
    {
        // this naive implementation functions as an interface to impose the correct cluster from Python without having to call the RomB&S
        // (as there is not possible to access it as easily as before)
        public:
            KRATOS_CLASS_POINTER_DEFINITION(MultipleBasesManager);

            MultipleBasesManager(): current_basis_index(0){
                std::cout<<"no args constructor called"<<std::endl<<std::endl;
            }

            MultipleBasesManager(const int &this_many_clusters): current_basis_index(0) {
                std::cout<<"args constructor called"<<std::endl<<std::endl;
                number_of_bases = this_many_clusters;

            }

            void HardSetCurrentCluster(int this_index){
                current_basis_index = this_index;
            }

            int GetCurrentCluster(){
                return current_basis_index;
            }

            int current_basis_index;

        protected:


            // int current_basis_index{0};
            int number_of_bases;

    };









    class DistanceToClusters // deprecated implementation. Moreover, it seems to have a little bug in the update of the distance :/
    {
        // this implemenation of Local POD is based on (Nonlinear model order reduction based on local reduced‐order bases, David Amsallem et. al, 2012)
        public:
            KRATOS_CLASS_POINTER_DEFINITION(DistanceToClusters);

            DistanceToClusters(){
                std::cout<<"no args constructor called"<<std::endl<<std::endl;
            }

            DistanceToClusters(const int &this_many_clusters){
                std::cout<<"args constructor called"<<std::endl<<std::endl;
                number_of_bases = this_many_clusters;
                SetUpWVector();
                SetUpZMatrix();
            }


            void SetWEntry(const Vector &this_vector,const int &cluster, const int &m, const int &p){
                w[cluster][m][p] = this_vector;
            }

            Vector GetWEntry(const int &cluster, const int &m, const int &p){
                return w[cluster][m][p];
            }

            void SetZEntry(const double &this_double,const int &m, const int &p){
                z(m,p) = this_double;
            }

            Matrix GetZMatrix(){
                return z;
            }

            void SetUpWVector(){
                //creating 3-D std vector
                for(int i = 0; i < number_of_bases; i++){
                    std::vector<std::vector<Vector>> temporary;
                    for(int j = 0; j < number_of_bases; j++){
                        std::vector<Vector> temporary2;
                        for(int k = 0; k< number_of_bases; k++){
                            temporary2.push_back(Vector());
                        };
                        temporary.push_back(temporary2);
                    };
                    w.push_back(temporary);
                };
            }

            void SetUpZMatrix(){
                z.resize(number_of_bases, number_of_bases,false);
                //filling matrix manually, there might be a better way (although it's small)
                for(int i=0; i<number_of_bases; i++){
                    for(int j=0; j<number_of_bases; j++){
                        z(i,j) = -1;
                    };
                };
            }

            void UpdateZMatrix(const Vector &current_dq)
            {
                for(int m=0; m<number_of_bases; m++){
                    int k=m+1;
                    for(int p=k; p<number_of_bases; p++){
                        z(m,p) = z(m,p) + inner_prod(w[current_basis_index][m][p], current_dq);
                        z(p,m) = -z(m,p);
                    };
                };

            }

            void UpdateCurrentCluster(){
                // std::cout<<"current basis index before update is: "<<current_basis_index<<std::endl;
                //KRATOS_WATCH(z);
                bool index_found = false;
                int the_index = -1;
                while(index_found == false){
                    the_index++;
                    index_found = true;
                    for(int i=0;i<number_of_bases; i++){
                        if(z(the_index,i)>0){
                            index_found = false;
                        };
                    };
                };
                current_basis_index = the_index;
                std::cout<<"current basis index after update is: "<<current_basis_index<<std::endl;
            }


            void HardSetCurrentCluster(int this_index){
                current_basis_index = this_index;
            }


            int GetCurrentCluster(){
                return current_basis_index;
            }

        protected:


            int current_basis_index{5};
            int number_of_bases;
            std::vector<std::vector<std::vector<Vector>>> w;
            Matrix z;


    };








} // namespace Kratos



#endif // ROM_BASES_H_INCLUDED  defined



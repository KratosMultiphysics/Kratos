//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "processes/tetrahedra_mesh_edge_swapping_process.h"
#include "geometries/tetrahedra_3d_4.h"
#include "processes/element_erase_process.h"



namespace Kratos
{
	namespace Internals {
		class EdgeSwappingCase {
		public:
			EdgeSwappingCase() {};
			EdgeSwappingCase(std::vector<std::size_t> const&& TheTrianglesIndices) : mTriangleIndices(TheTrianglesIndices), mMinQuality(2.00) {}
			std::size_t GetTringleIndex(std::size_t Index) const { return mTriangleIndices[Index]; }
			double GetMinQuality() const { return mMinQuality; }
			void SetMinQuality(double NewMinQuality) { mMinQuality = NewMinQuality; }
		private:
			std::vector<std::size_t> mTriangleIndices;
			double mMinQuality;
		};
		template<std::size_t TNumberOfCases, std::size_t TNumberOfTriangles, std::size_t TNumberOfTrianglesPerCase>
		class EdgeSwappingCases {
		public:
			using TriangleConnectivityType = std::array<std::size_t, 3>;
			static std::size_t GetNumberOfCases() { return TNumberOfCases; }
			static std::size_t NumberOfTriangles() { return TNumberOfTriangles; }
			static std::size_t NumberOfTrianglesPerCase() { return TNumberOfTrianglesPerCase; }
			std::array<EdgeSwappingCase, TNumberOfCases> const& GetCases() { return mCases; }
			TriangleConnectivityType const& GetTriangleConectivity(std::size_t TheIndex) { return mTriangles[TheIndex]; }
			void SetTetrahedraForCase(EdgeSwappingCase const& TheCase, std::size_t TriangleIndex, TetrahedraEdgeShell& EdgeShell, Tetrahedra3D4<Node<3>>& rTetrahedra1, Tetrahedra3D4<Node<3>>& rTetrahedra2) {
				auto const& triangle = GetTriangleConectivity(TheCase.GetTringleIndex(TriangleIndex));
				rTetrahedra1(0) = EdgeShell.Point1();
				rTetrahedra1(1) = EdgeShell.ShellPoint(triangle[0]);
				rTetrahedra1(2) = EdgeShell.ShellPoint(triangle[1]);
				rTetrahedra1(3) = EdgeShell.ShellPoint(triangle[2]);

				rTetrahedra2(0) = EdgeShell.Point2();
				rTetrahedra2(1) = EdgeShell.ShellPoint(triangle[0]);
				rTetrahedra2(2) = EdgeShell.ShellPoint(triangle[2]);
				rTetrahedra2(3) = EdgeShell.ShellPoint(triangle[1]);
			}
		protected:
			EdgeSwappingCases() {}
			std::array<TriangleConnectivityType, TNumberOfTriangles> mTriangles;
			std::array<EdgeSwappingCase, TNumberOfCases> mCases;
		};

		class EdgeSwappingCases3 : public EdgeSwappingCases<1, 1, 1> {
		public:
			EdgeSwappingCases3() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0 });
				mTriangles[0] = { 0,1,2 };
			}
		};

		class EdgeSwappingCases4 : public EdgeSwappingCases<2, 4, 2> {
		public:
			EdgeSwappingCases4() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0,1 });
				mCases[1] = EdgeSwappingCase({ 2,3 });
				mTriangles[0] = { 0,1,2 };
				mTriangles[1] = { 0,2,3 };
				mTriangles[2] = { 1,2,3 };
				mTriangles[3] = { 0,1,3 };
			}
		};

		class EdgeSwappingCases5 : public EdgeSwappingCases<5, 10, 3> {
		public:
			EdgeSwappingCases5() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({  0, 1, 2  });
				mCases[1] = EdgeSwappingCase({  0, 3, 4  });
				mCases[2] = EdgeSwappingCase({  5, 6, 2  });
				mCases[3] = EdgeSwappingCase({  5, 7, 8  });
				mCases[4] = EdgeSwappingCase({  3, 9, 8  });
				mTriangles[0] = { 0, 1, 2 };
				mTriangles[1] = { 0, 2, 3 };
				mTriangles[2] = {  0, 3, 4  };
				mTriangles[3] = {  2, 3, 4  };
				mTriangles[4] = {  0, 2, 4  };
				mTriangles[5] = {  1, 2, 3  };
				mTriangles[6] = {  0, 1, 3  };
				mTriangles[7] = {  1, 3, 4  };
				mTriangles[8] = {  0, 1, 4  };
				mTriangles[9] = {  1, 2, 4  };
			}
		};

		class  EdgeSwappingCases6 : public EdgeSwappingCases< 14, 20, 4 > {
		public:
			EdgeSwappingCases6() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0, 1, 2, 3 });
				mCases[1] = EdgeSwappingCase({ 0, 1, 4, 5 });
				mCases[2] = EdgeSwappingCase({ 0, 6, 7, 3 });
				mCases[3] = EdgeSwappingCase({ 0, 6, 8, 9 });
				mCases[4] = EdgeSwappingCase({ 0, 4, 10, 9 });
				mCases[5] = EdgeSwappingCase({ 11, 12, 2, 3 });
				mCases[6] = EdgeSwappingCase({ 11, 12, 4, 5 });
				mCases[7] = EdgeSwappingCase({ 11, 13, 14, 3 });
				mCases[8] = EdgeSwappingCase({ 11, 13, 15, 16 });
				mCases[9] = EdgeSwappingCase({ 11, 4, 17, 16 });
				mCases[10] = EdgeSwappingCase({ 6, 18, 14, 3 });
				mCases[11] = EdgeSwappingCase({ 6, 18, 15, 16 });
				mCases[12] = EdgeSwappingCase({ 6, 8, 19, 16 });
				mCases[13] = EdgeSwappingCase({ 4, 10, 19, 16 });
				mTriangles[0] = { 0 , 1 , 2 };
				mTriangles[1] = { 0 , 2 , 3 };
				mTriangles[2] = { 0 , 3 , 4 };
				mTriangles[3] = { 0 , 4 , 5 };
				mTriangles[4] = { 3 , 4 , 5 };
				mTriangles[5] = { 0 , 3 , 5 };
				mTriangles[6] = { 2 , 3 , 4 };
				mTriangles[7] = { 0 , 2 , 4 };
				mTriangles[8] = { 2 , 4 , 5 };
				mTriangles[9] = { 0 , 2 , 5 };
				mTriangles[10] = { 2 , 3 , 5 };
				mTriangles[11] = { 1 , 2 , 3 };
				mTriangles[12] = { 0 , 1 , 3 };
				mTriangles[13] = { 1 , 3 , 4 };
				mTriangles[14] = { 0 , 1 , 4 };
				mTriangles[15] = { 1 , 4 , 5 };
				mTriangles[16] = { 0 , 1 , 5 };
				mTriangles[17] = { 1 , 3 , 5 };
				mTriangles[18] = { 1 , 2 , 4 };
				mTriangles[19] = { 1 , 2 , 5 };
			}
		};

		class  EdgeSwappingCases7 : public EdgeSwappingCases< 42, 35, 5 > {
		public:
			EdgeSwappingCases7() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0, 1, 2, 3, 4 });
				mCases[1] = EdgeSwappingCase({ 0, 1, 2, 5, 6 });
				mCases[2] = EdgeSwappingCase({ 0, 1, 7, 8, 4 });
				mCases[3] = EdgeSwappingCase({ 0, 1, 7, 9, 10 });
				mCases[4] = EdgeSwappingCase({ 0, 1, 5, 11, 10 });
				mCases[5] = EdgeSwappingCase({ 0, 12, 13, 3, 4 });
				mCases[6] = EdgeSwappingCase({ 0, 12, 13, 5, 6 });
				mCases[7] = EdgeSwappingCase({ 0, 12, 14, 15, 4 });
				mCases[8] = EdgeSwappingCase({ 0, 12, 14, 16, 17 });
				mCases[9] = EdgeSwappingCase({ 0, 12, 5, 18, 17 });
				mCases[10] = EdgeSwappingCase({ 0, 7, 19, 15, 4 });
				mCases[11] = EdgeSwappingCase({ 0, 7, 19, 16, 17 });
				mCases[12] = EdgeSwappingCase({ 0, 7, 9, 20, 17 });
				mCases[13] = EdgeSwappingCase({ 0, 5, 11, 20, 17 });
				mCases[14] = EdgeSwappingCase({ 21, 22, 2, 3, 4 });
				mCases[15] = EdgeSwappingCase({ 21, 22, 2, 5, 6 });
				mCases[16] = EdgeSwappingCase({ 21, 22, 7, 8, 4 });
				mCases[17] = EdgeSwappingCase({ 21, 22, 7, 9, 10 });
				mCases[18] = EdgeSwappingCase({ 21, 22, 5, 11, 10 });
				mCases[19] = EdgeSwappingCase({ 21, 23, 24, 3, 4 });
				mCases[20] = EdgeSwappingCase({ 21, 23, 24, 5, 6 });
				mCases[21] = EdgeSwappingCase({ 21, 23, 25, 26, 4 });
				mCases[22] = EdgeSwappingCase({ 21, 23, 25, 27, 28 });
				mCases[23] = EdgeSwappingCase({ 21, 23, 5, 29, 28 });
				mCases[24] = EdgeSwappingCase({ 21, 7, 30, 26, 4 });
				mCases[25] = EdgeSwappingCase({ 21, 7, 30, 27, 28 });
				mCases[26] = EdgeSwappingCase({ 21, 7, 9, 31, 28 });
				mCases[27] = EdgeSwappingCase({ 21, 5, 11, 31, 28 });
				mCases[28] = EdgeSwappingCase({ 12, 32, 24, 3, 4 });
				mCases[29] = EdgeSwappingCase({ 12, 32, 24, 5, 6 });
				mCases[30] = EdgeSwappingCase({ 12, 32, 25, 26, 4 });
				mCases[31] = EdgeSwappingCase({ 12, 32, 25, 27, 28 });
				mCases[32] = EdgeSwappingCase({ 12, 32, 5, 29, 28 });
				mCases[33] = EdgeSwappingCase({ 12, 14, 33, 26, 4 });
				mCases[34] = EdgeSwappingCase({ 12, 14, 33, 27, 28 });
				mCases[35] = EdgeSwappingCase({ 12, 14, 16, 34, 28 });
				mCases[36] = EdgeSwappingCase({ 12, 5, 18, 34, 28 });
				mCases[37] = EdgeSwappingCase({ 7, 19, 33, 26, 4 });
				mCases[38] = EdgeSwappingCase({ 7, 19, 33, 27, 28 });
				mCases[39] = EdgeSwappingCase({ 7, 19, 16, 34, 28 });
				mCases[40] = EdgeSwappingCase({ 7, 9, 20, 34, 28 });
				mCases[41] = EdgeSwappingCase({ 5, 11, 20, 34, 28 });
				mTriangles[0] = { 0 , 1 , 2 };
				mTriangles[1] = { 0 , 2 , 3 };
				mTriangles[2] = { 0 , 3 , 4 };
				mTriangles[3] = { 0 , 4 , 5 };
				mTriangles[4] = { 0 , 5 , 6 };
				mTriangles[5] = { 4 , 5 , 6 };
				mTriangles[6] = { 0 , 4 , 6 };
				mTriangles[7] = { 3 , 4 , 5 };
				mTriangles[8] = { 0 , 3 , 5 };
				mTriangles[9] = { 3 , 5 , 6 };
				mTriangles[10] = { 0 , 3 , 6 };
				mTriangles[11] = { 3 , 4 , 6 };
				mTriangles[12] = { 2 , 3 , 4 };
				mTriangles[13] = { 0 , 2 , 4 };
				mTriangles[14] = { 2 , 4 , 5 };
				mTriangles[15] = { 0 , 2 , 5 };
				mTriangles[16] = { 2 , 5 , 6 };
				mTriangles[17] = { 0 , 2 , 6 };
				mTriangles[18] = { 2 , 4 , 6 };
				mTriangles[19] = { 2 , 3 , 5 };
				mTriangles[20] = { 2 , 3 , 6 };
				mTriangles[21] = { 1 , 2 , 3 };
				mTriangles[22] = { 0 , 1 , 3 };
				mTriangles[23] = { 1 , 3 , 4 };
				mTriangles[24] = { 0 , 1 , 4 };
				mTriangles[25] = { 1 , 4 , 5 };
				mTriangles[26] = { 0 , 1 , 5 };
				mTriangles[27] = { 1 , 5 , 6 };
				mTriangles[28] = { 0 , 1 , 6 };
				mTriangles[29] = { 1 , 4 , 6 };
				mTriangles[30] = { 1 , 3 , 5 };
				mTriangles[31] = { 1 , 3 , 6 };
				mTriangles[32] = { 1 , 2 , 4 };
				mTriangles[33] = { 1 , 2 , 5 };
				mTriangles[34] = { 1 , 2 , 6 };
			}
		};

		class  EdgeSwappingCases8 : public EdgeSwappingCases< 132, 56, 6 > {
		public:
			EdgeSwappingCases8() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0, 1, 2, 3, 4, 5 });
				mCases[1] = EdgeSwappingCase({ 0, 1, 2, 3, 6, 7 });
				mCases[2] = EdgeSwappingCase({ 0, 1, 2, 8, 9, 5 });
				mCases[3] = EdgeSwappingCase({ 0, 1, 2, 8, 10, 11 });
				mCases[4] = EdgeSwappingCase({ 0, 1, 2, 6, 12, 11 });
				mCases[5] = EdgeSwappingCase({ 0, 1, 13, 14, 4, 5 });
				mCases[6] = EdgeSwappingCase({ 0, 1, 13, 14, 6, 7 });
				mCases[7] = EdgeSwappingCase({ 0, 1, 13, 15, 16, 5 });
				mCases[8] = EdgeSwappingCase({ 0, 1, 13, 15, 17, 18 });
				mCases[9] = EdgeSwappingCase({ 0, 1, 13, 6, 19, 18 });
				mCases[10] = EdgeSwappingCase({ 0, 1, 8, 20, 16, 5 });
				mCases[11] = EdgeSwappingCase({ 0, 1, 8, 20, 17, 18 });
				mCases[12] = EdgeSwappingCase({ 0, 1, 8, 10, 21, 18 });
				mCases[13] = EdgeSwappingCase({ 0, 1, 6, 12, 21, 18 });
				mCases[14] = EdgeSwappingCase({ 0, 22, 23, 3, 4, 5 });
				mCases[15] = EdgeSwappingCase({ 0, 22, 23, 3, 6, 7 });
				mCases[16] = EdgeSwappingCase({ 0, 22, 23, 8, 9, 5 });
				mCases[17] = EdgeSwappingCase({ 0, 22, 23, 8, 10, 11 });
				mCases[18] = EdgeSwappingCase({ 0, 22, 23, 6, 12, 11 });
				mCases[19] = EdgeSwappingCase({ 0, 22, 24, 25, 4, 5 });
				mCases[20] = EdgeSwappingCase({ 0, 22, 24, 25, 6, 7 });
				mCases[21] = EdgeSwappingCase({ 0, 22, 24, 26, 27, 5 });
				mCases[22] = EdgeSwappingCase({ 0, 22, 24, 26, 28, 29 });
				mCases[23] = EdgeSwappingCase({ 0, 22, 24, 6, 30, 29 });
				mCases[24] = EdgeSwappingCase({ 0, 22, 8, 31, 27, 5 });
				mCases[25] = EdgeSwappingCase({ 0, 22, 8, 31, 28, 29 });
				mCases[26] = EdgeSwappingCase({ 0, 22, 8, 10, 32, 29 });
				mCases[27] = EdgeSwappingCase({ 0, 22, 6, 12, 32, 29 });
				mCases[28] = EdgeSwappingCase({ 0, 13, 33, 25, 4, 5 });
				mCases[29] = EdgeSwappingCase({ 0, 13, 33, 25, 6, 7 });
				mCases[30] = EdgeSwappingCase({ 0, 13, 33, 26, 27, 5 });
				mCases[31] = EdgeSwappingCase({ 0, 13, 33, 26, 28, 29 });
				mCases[32] = EdgeSwappingCase({ 0, 13, 33, 6, 30, 29 });
				mCases[33] = EdgeSwappingCase({ 0, 13, 15, 34, 27, 5 });
				mCases[34] = EdgeSwappingCase({ 0, 13, 15, 34, 28, 29 });
				mCases[35] = EdgeSwappingCase({ 0, 13, 15, 17, 35, 29 });
				mCases[36] = EdgeSwappingCase({ 0, 13, 6, 19, 35, 29 });
				mCases[37] = EdgeSwappingCase({ 0, 8, 20, 34, 27, 5 });
				mCases[38] = EdgeSwappingCase({ 0, 8, 20, 34, 28, 29 });
				mCases[39] = EdgeSwappingCase({ 0, 8, 20, 17, 35, 29 });
				mCases[40] = EdgeSwappingCase({ 0, 8, 10, 21, 35, 29 });
				mCases[41] = EdgeSwappingCase({ 0, 6, 12, 21, 35, 29 });
				mCases[42] = EdgeSwappingCase({ 36, 37, 2, 3, 4, 5 });
				mCases[43] = EdgeSwappingCase({ 36, 37, 2, 3, 6, 7 });
				mCases[44] = EdgeSwappingCase({ 36, 37, 2, 8, 9, 5 });
				mCases[45] = EdgeSwappingCase({ 36, 37, 2, 8, 10, 11 });
				mCases[46] = EdgeSwappingCase({ 36, 37, 2, 6, 12, 11 });
				mCases[47] = EdgeSwappingCase({ 36, 37, 13, 14, 4, 5 });
				mCases[48] = EdgeSwappingCase({ 36, 37, 13, 14, 6, 7 });
				mCases[49] = EdgeSwappingCase({ 36, 37, 13, 15, 16, 5 });
				mCases[50] = EdgeSwappingCase({ 36, 37, 13, 15, 17, 18 });
				mCases[51] = EdgeSwappingCase({ 36, 37, 13, 6, 19, 18 });
				mCases[52] = EdgeSwappingCase({ 36, 37, 8, 20, 16, 5 });
				mCases[53] = EdgeSwappingCase({ 36, 37, 8, 20, 17, 18 });
				mCases[54] = EdgeSwappingCase({ 36, 37, 8, 10, 21, 18 });
				mCases[55] = EdgeSwappingCase({ 36, 37, 6, 12, 21, 18 });
				mCases[56] = EdgeSwappingCase({ 36, 38, 39, 3, 4, 5 });
				mCases[57] = EdgeSwappingCase({ 36, 38, 39, 3, 6, 7 });
				mCases[58] = EdgeSwappingCase({ 36, 38, 39, 8, 9, 5 });
				mCases[59] = EdgeSwappingCase({ 36, 38, 39, 8, 10, 11 });
				mCases[60] = EdgeSwappingCase({ 36, 38, 39, 6, 12, 11 });
				mCases[61] = EdgeSwappingCase({ 36, 38, 40, 41, 4, 5 });
				mCases[62] = EdgeSwappingCase({ 36, 38, 40, 41, 6, 7 });
				mCases[63] = EdgeSwappingCase({ 36, 38, 40, 42, 43, 5 });
				mCases[64] = EdgeSwappingCase({ 36, 38, 40, 42, 44, 45 });
				mCases[65] = EdgeSwappingCase({ 36, 38, 40, 6, 46, 45 });
				mCases[66] = EdgeSwappingCase({ 36, 38, 8, 47, 43, 5 });
				mCases[67] = EdgeSwappingCase({ 36, 38, 8, 47, 44, 45 });
				mCases[68] = EdgeSwappingCase({ 36, 38, 8, 10, 48, 45 });
				mCases[69] = EdgeSwappingCase({ 36, 38, 6, 12, 48, 45 });
				mCases[70] = EdgeSwappingCase({ 36, 13, 49, 41, 4, 5 });
				mCases[71] = EdgeSwappingCase({ 36, 13, 49, 41, 6, 7 });
				mCases[72] = EdgeSwappingCase({ 36, 13, 49, 42, 43, 5 });
				mCases[73] = EdgeSwappingCase({ 36, 13, 49, 42, 44, 45 });
				mCases[74] = EdgeSwappingCase({ 36, 13, 49, 6, 46, 45 });
				mCases[75] = EdgeSwappingCase({ 36, 13, 15, 50, 43, 5 });
				mCases[76] = EdgeSwappingCase({ 36, 13, 15, 50, 44, 45 });
				mCases[77] = EdgeSwappingCase({ 36, 13, 15, 17, 51, 45 });
				mCases[78] = EdgeSwappingCase({ 36, 13, 6, 19, 51, 45 });
				mCases[79] = EdgeSwappingCase({ 36, 8, 20, 50, 43, 5 });
				mCases[80] = EdgeSwappingCase({ 36, 8, 20, 50, 44, 45 });
				mCases[81] = EdgeSwappingCase({ 36, 8, 20, 17, 51, 45 });
				mCases[82] = EdgeSwappingCase({ 36, 8, 10, 21, 51, 45 });
				mCases[83] = EdgeSwappingCase({ 36, 6, 12, 21, 51, 45 });
				mCases[84] = EdgeSwappingCase({ 22, 52, 39, 3, 4, 5 });
				mCases[85] = EdgeSwappingCase({ 22, 52, 39, 3, 6, 7 });
				mCases[86] = EdgeSwappingCase({ 22, 52, 39, 8, 9, 5 });
				mCases[87] = EdgeSwappingCase({ 22, 52, 39, 8, 10, 11 });
				mCases[88] = EdgeSwappingCase({ 22, 52, 39, 6, 12, 11 });
				mCases[89] = EdgeSwappingCase({ 22, 52, 40, 41, 4, 5 });
				mCases[90] = EdgeSwappingCase({ 22, 52, 40, 41, 6, 7 });
				mCases[91] = EdgeSwappingCase({ 22, 52, 40, 42, 43, 5 });
				mCases[92] = EdgeSwappingCase({ 22, 52, 40, 42, 44, 45 });
				mCases[93] = EdgeSwappingCase({ 22, 52, 40, 6, 46, 45 });
				mCases[94] = EdgeSwappingCase({ 22, 52, 8, 47, 43, 5 });
				mCases[95] = EdgeSwappingCase({ 22, 52, 8, 47, 44, 45 });
				mCases[96] = EdgeSwappingCase({ 22, 52, 8, 10, 48, 45 });
				mCases[97] = EdgeSwappingCase({ 22, 52, 6, 12, 48, 45 });
				mCases[98] = EdgeSwappingCase({ 22, 24, 53, 41, 4, 5 });
				mCases[99] = EdgeSwappingCase({ 22, 24, 53, 41, 6, 7 });
				mCases[100] = EdgeSwappingCase({ 22, 24, 53, 42, 43, 5 });
				mCases[101] = EdgeSwappingCase({ 22, 24, 53, 42, 44, 45 });
				mCases[102] = EdgeSwappingCase({ 22, 24, 53, 6, 46, 45 });
				mCases[103] = EdgeSwappingCase({ 22, 24, 26, 54, 43, 5 });
				mCases[104] = EdgeSwappingCase({ 22, 24, 26, 54, 44, 45 });
				mCases[105] = EdgeSwappingCase({ 22, 24, 26, 28, 55, 45 });
				mCases[106] = EdgeSwappingCase({ 22, 24, 6, 30, 55, 45 });
				mCases[107] = EdgeSwappingCase({ 22, 8, 31, 54, 43, 5 });
				mCases[108] = EdgeSwappingCase({ 22, 8, 31, 54, 44, 45 });
				mCases[109] = EdgeSwappingCase({ 22, 8, 31, 28, 55, 45 });
				mCases[110] = EdgeSwappingCase({ 22, 8, 10, 32, 55, 45 });
				mCases[111] = EdgeSwappingCase({ 22, 6, 12, 32, 55, 45 });
				mCases[112] = EdgeSwappingCase({ 13, 33, 53, 41, 4, 5 });
				mCases[113] = EdgeSwappingCase({ 13, 33, 53, 41, 6, 7 });
				mCases[114] = EdgeSwappingCase({ 13, 33, 53, 42, 43, 5 });
				mCases[115] = EdgeSwappingCase({ 13, 33, 53, 42, 44, 45 });
				mCases[116] = EdgeSwappingCase({ 13, 33, 53, 6, 46, 45 });
				mCases[117] = EdgeSwappingCase({ 13, 33, 26, 54, 43, 5 });
				mCases[118] = EdgeSwappingCase({ 13, 33, 26, 54, 44, 45 });
				mCases[119] = EdgeSwappingCase({ 13, 33, 26, 28, 55, 45 });
				mCases[120] = EdgeSwappingCase({ 13, 33, 6, 30, 55, 45 });
				mCases[121] = EdgeSwappingCase({ 13, 15, 34, 54, 43, 5 });
				mCases[122] = EdgeSwappingCase({ 13, 15, 34, 54, 44, 45 });
				mCases[123] = EdgeSwappingCase({ 13, 15, 34, 28, 55, 45 });
				mCases[124] = EdgeSwappingCase({ 13, 15, 17, 35, 55, 45 });
				mCases[125] = EdgeSwappingCase({ 13, 6, 19, 35, 55, 45 });
				mCases[126] = EdgeSwappingCase({ 8, 20, 34, 54, 43, 5 });
				mCases[127] = EdgeSwappingCase({ 8, 20, 34, 54, 44, 45 });
				mCases[128] = EdgeSwappingCase({ 8, 20, 34, 28, 55, 45 });
				mCases[129] = EdgeSwappingCase({ 8, 20, 17, 35, 55, 45 });
				mCases[130] = EdgeSwappingCase({ 8, 10, 21, 35, 55, 45 });
				mCases[131] = EdgeSwappingCase({ 6, 12, 21, 35, 55, 45 });
				mTriangles[0] = { 0 , 1 , 2 };
				mTriangles[1] = { 0 , 2 , 3 };
				mTriangles[2] = { 0 , 3 , 4 };
				mTriangles[3] = { 0 , 4 , 5 };
				mTriangles[4] = { 0 , 5 , 6 };
				mTriangles[5] = { 0 , 6 , 7 };
				mTriangles[6] = { 5 , 6 , 7 };
				mTriangles[7] = { 0 , 5 , 7 };
				mTriangles[8] = { 4 , 5 , 6 };
				mTriangles[9] = { 0 , 4 , 6 };
				mTriangles[10] = { 4 , 6 , 7 };
				mTriangles[11] = { 0 , 4 , 7 };
				mTriangles[12] = { 4 , 5 , 7 };
				mTriangles[13] = { 3 , 4 , 5 };
				mTriangles[14] = { 0 , 3 , 5 };
				mTriangles[15] = { 3 , 5 , 6 };
				mTriangles[16] = { 0 , 3 , 6 };
				mTriangles[17] = { 3 , 6 , 7 };
				mTriangles[18] = { 0 , 3 , 7 };
				mTriangles[19] = { 3 , 5 , 7 };
				mTriangles[20] = { 3 , 4 , 6 };
				mTriangles[21] = { 3 , 4 , 7 };
				mTriangles[22] = { 2 , 3 , 4 };
				mTriangles[23] = { 0 , 2 , 4 };
				mTriangles[24] = { 2 , 4 , 5 };
				mTriangles[25] = { 0 , 2 , 5 };
				mTriangles[26] = { 2 , 5 , 6 };
				mTriangles[27] = { 0 , 2 , 6 };
				mTriangles[28] = { 2 , 6 , 7 };
				mTriangles[29] = { 0 , 2 , 7 };
				mTriangles[30] = { 2 , 5 , 7 };
				mTriangles[31] = { 2 , 4 , 6 };
				mTriangles[32] = { 2 , 4 , 7 };
				mTriangles[33] = { 2 , 3 , 5 };
				mTriangles[34] = { 2 , 3 , 6 };
				mTriangles[35] = { 2 , 3 , 7 };
				mTriangles[36] = { 1 , 2 , 3 };
				mTriangles[37] = { 0 , 1 , 3 };
				mTriangles[38] = { 1 , 3 , 4 };
				mTriangles[39] = { 0 , 1 , 4 };
				mTriangles[40] = { 1 , 4 , 5 };
				mTriangles[41] = { 0 , 1 , 5 };
				mTriangles[42] = { 1 , 5 , 6 };
				mTriangles[43] = { 0 , 1 , 6 };
				mTriangles[44] = { 1 , 6 , 7 };
				mTriangles[45] = { 0 , 1 , 7 };
				mTriangles[46] = { 1 , 5 , 7 };
				mTriangles[47] = { 1 , 4 , 6 };
				mTriangles[48] = { 1 , 4 , 7 };
				mTriangles[49] = { 1 , 3 , 5 };
				mTriangles[50] = { 1 , 3 , 6 };
				mTriangles[51] = { 1 , 3 , 7 };
				mTriangles[52] = { 1 , 2 , 4 };
				mTriangles[53] = { 1 , 2 , 5 };
				mTriangles[54] = { 1 , 2 , 6 };
				mTriangles[55] = { 1 , 2 , 7 };
			}
		};

	}


TetrahedraMeshEdgeSwappingProcess::TetrahedraMeshEdgeSwappingProcess(ModelPart & rModelPart): mrModelPart(rModelPart), mEdges(){

}

TetrahedraMeshEdgeSwappingProcess::~TetrahedraMeshEdgeSwappingProcess(){

}

void TetrahedraMeshEdgeSwappingProcess::Execute(){
	std::cout << std::endl;

	constexpr int tetrahedra_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
	for(auto i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++){
		auto& element_geometry = i_element->GetGeometry();
	 	//std::cout << "Processing element #" << i_element->Id() << "[" << element_geometry[0].Id() << "," << element_geometry[1].Id() << ","
	  //<< element_geometry[2].Id() << "," << element_geometry[3].Id()	<< "]" << std::endl;
		for(int i = 0 ; i < 6 ; i++){
				auto i_edge = mEdges.find(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])));
				if(i_edge == mEdges.end())
					i_edge = mEdges.emplace(std::make_pair(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])), TetrahedraEdgeShell(element_geometry(tetrahedra_edges[i][0]),element_geometry(tetrahedra_edges[i][1])))).first;

				//std::cout << "Before: edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has ";
				//i_edge->second.PrintData(std::cout);
				//std::cout << std::endl;
				i_edge->second.AddElement((*i_element.base()).get(),i);
				//std::cout << "After : edge " << i_edge->first.GetPoint1()->Id() << " -> " << i_edge->first.GetPoint2()->Id() << " has ";
				//i_edge->second.PrintData(std::cout);
				//std::cout << std::endl;
				auto i_edge1 = mEdges.find(Edge(&(element_geometry[tetrahedra_edges[i][0]]), &(element_geometry[tetrahedra_edges[i][1]])));
				//std::cout << "After :edge1 " << i_edge1->first.GetPoint1()->Id() << " -> " << i_edge1->first.GetPoint2()->Id() << " has ";
				//i_edge1->second.PrintData(std::cout);
				//std::cout << std::endl;
			}
	}
	for (auto& edge : mEdges) {
		//std::cout << "Before: edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has ";
		//edge.second.PrintData(std::cout);
		//std::cout << std::endl;
		edge.second.AddShellPoints();
		//std::cout << "After : edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has ";
		//edge.second.PrintData(std::cout);
		//std::cout << std::endl;
	}

	std::array<int, 100> edge_counter;
	for (auto& i : edge_counter)
		i = 0;

	for(auto& edge : mEdges){
		auto size = edge.second.GetNumberOfShellPoints();
		auto tet_numbers = edge.second.GetNumberOfTetrahedra();
		if(size < 100)
			edge_counter[size]++;
		//if ((tet_numbers == 1 && size !=2) || (tet_numbers == 2 && size != 3) || (tet_numbers == 3 && size != 3) || (tet_numbers == 4 && size != 4) || (tet_numbers == 5 && size != 5) || (tet_numbers == 6 && size != 6))
		//{
		//	std::cout << "edge " << edge.first.GetPoint1()->Id() << " -> " << edge.first.GetPoint2()->Id() << " has " << edge.second.GetNumberOfTetrahedra() << " tetrahedras and " << size << " points ";
		//	edge.second.PrintData(std::cout);
		//	std::cout << std::endl;
		//}
	}
	for(std::size_t i = 0 ; i < edge_counter.size() ; i++)
		if(edge_counter[i] > 0)
			std::cout << edge_counter[i] << " edges with " << i << " points" << std::endl;
	KRATOS_WATCH(mEdges.size());

	for (auto& edge : mEdges) {
		if (edge.second.IsClosed()) {
			if (!(edge.second.IsModified())) {
				if (edge.second.GetNumberOfShellPoints() == 3)
					EdgeSwapping3(edge.second);
				if (edge.second.GetNumberOfShellPoints() == 4)
					EdgeSwapping<Internals::EdgeSwappingCases4>(edge.second);
				if (edge.second.GetNumberOfShellPoints() == 5)
					EdgeSwapping<Internals::EdgeSwappingCases5>(edge.second);
				if (edge.second.GetNumberOfShellPoints() == 6)
					EdgeSwapping<Internals::EdgeSwappingCases6>(edge.second);
				if (edge.second.GetNumberOfShellPoints() == 7)
					EdgeSwapping<Internals::EdgeSwappingCases7>(edge.second);
			}
		}

	}
	ElementEraseProcess(mrModelPart).Execute();
}

std::string TetrahedraMeshEdgeSwappingProcess::Info() const{
	 return "TetrahedraMeshEdgeSwappingProcess";
 }

/// Print information about this object.
void TetrahedraMeshEdgeSwappingProcess::PrintInfo(std::ostream& rOStream) const {
	rOStream << Info();
}

/// Print object's data.
void TetrahedraMeshEdgeSwappingProcess::PrintData(std::ostream& rOStream) const {

}

void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping3(TetrahedraEdgeShell & EdgeShell) {
	Internals::EdgeSwappingCases3 SwappingCases;
	auto swapping_case = SwappingCases.GetCases()[0];
	auto const& triangle = SwappingCases.GetTriangleConectivity(swapping_case.GetTringleIndex(0));
	Tetrahedra3D4<Node<3>> tetrahedra_1(EdgeShell.Point1(), EdgeShell.ShellPoint(triangle[0]), EdgeShell.ShellPoint(triangle[1]), EdgeShell.ShellPoint(triangle[2]));
	Tetrahedra3D4<Node<3>> tetrahedra_2(EdgeShell.Point2(), EdgeShell.ShellPoint(triangle[0]), EdgeShell.ShellPoint(triangle[2]), EdgeShell.ShellPoint(triangle[1]));
	auto quality_criteria = Geometry<Node<3> >::QualityCriteria::VOLUME_TO_EDGE_LENGTH;

	double original_min_quality = EdgeShell.CalculateMinQuality(quality_criteria);
	double min_quality = std::min(tetrahedra_1.Quality(quality_criteria), tetrahedra_2.Quality(quality_criteria));
	if (min_quality > original_min_quality) {
		EdgeShell.pGetElement(0)->GetGeometry() = tetrahedra_1;
		EdgeShell.pGetElement(1)->GetGeometry() = tetrahedra_2;
		EdgeShell.pGetElement(2)->Set(TO_ERASE);
	}
	//else
	//	std::cout << min_quality << " is worst respect to " << original_min_quality << std::endl;
}

void TetrahedraMeshEdgeSwappingProcess::EdgeSwapping4(TetrahedraEdgeShell & EdgeShell) {
	EdgeSwapping<Internals::EdgeSwappingCases4>(EdgeShell);
}


}  // namespace Kratos.

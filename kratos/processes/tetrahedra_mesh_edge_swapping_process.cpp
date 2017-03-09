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
			using EdgeSwappingCaseType = EdgeSwappingCase;
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

		class  EdgeSwappingCases9 : public EdgeSwappingCases< 429, 84, 7 > {
		public:
			EdgeSwappingCases9() : EdgeSwappingCases() {
				mCases[0] = EdgeSwappingCase({ 0, 1, 2, 3, 4, 5, 6 });
				mCases[1] = EdgeSwappingCase({ 0, 1, 2, 3, 4, 7, 8 });
				mCases[2] = EdgeSwappingCase({ 0, 1, 2, 3, 9, 10, 6 });
				mCases[3] = EdgeSwappingCase({ 0, 1, 2, 3, 9, 11, 12 });
				mCases[4] = EdgeSwappingCase({ 0, 1, 2, 3, 7, 13, 12 });
				mCases[5] = EdgeSwappingCase({ 0, 1, 2, 14, 15, 5, 6 });
				mCases[6] = EdgeSwappingCase({ 0, 1, 2, 14, 15, 7, 8 });
				mCases[7] = EdgeSwappingCase({ 0, 1, 2, 14, 16, 17, 6 });
				mCases[8] = EdgeSwappingCase({ 0, 1, 2, 14, 16, 18, 19 });
				mCases[9] = EdgeSwappingCase({ 0, 1, 2, 14, 7, 20, 19 });
				mCases[10] = EdgeSwappingCase({ 0, 1, 2, 9, 21, 17, 6 });
				mCases[11] = EdgeSwappingCase({ 0, 1, 2, 9, 21, 18, 19 });
				mCases[12] = EdgeSwappingCase({ 0, 1, 2, 9, 11, 22, 19 });
				mCases[13] = EdgeSwappingCase({ 0, 1, 2, 7, 13, 22, 19 });
				mCases[14] = EdgeSwappingCase({ 0, 1, 23, 24, 4, 5, 6 });
				mCases[15] = EdgeSwappingCase({ 0, 1, 23, 24, 4, 7, 8 });
				mCases[16] = EdgeSwappingCase({ 0, 1, 23, 24, 9, 10, 6 });
				mCases[17] = EdgeSwappingCase({ 0, 1, 23, 24, 9, 11, 12 });
				mCases[18] = EdgeSwappingCase({ 0, 1, 23, 24, 7, 13, 12 });
				mCases[19] = EdgeSwappingCase({ 0, 1, 23, 25, 26, 5, 6 });
				mCases[20] = EdgeSwappingCase({ 0, 1, 23, 25, 26, 7, 8 });
				mCases[21] = EdgeSwappingCase({ 0, 1, 23, 25, 27, 28, 6 });
				mCases[22] = EdgeSwappingCase({ 0, 1, 23, 25, 27, 29, 30 });
				mCases[23] = EdgeSwappingCase({ 0, 1, 23, 25, 7, 31, 30 });
				mCases[24] = EdgeSwappingCase({ 0, 1, 23, 9, 32, 28, 6 });
				mCases[25] = EdgeSwappingCase({ 0, 1, 23, 9, 32, 29, 30 });
				mCases[26] = EdgeSwappingCase({ 0, 1, 23, 9, 11, 33, 30 });
				mCases[27] = EdgeSwappingCase({ 0, 1, 23, 7, 13, 33, 30 });
				mCases[28] = EdgeSwappingCase({ 0, 1, 14, 34, 26, 5, 6 });
				mCases[29] = EdgeSwappingCase({ 0, 1, 14, 34, 26, 7, 8 });
				mCases[30] = EdgeSwappingCase({ 0, 1, 14, 34, 27, 28, 6 });
				mCases[31] = EdgeSwappingCase({ 0, 1, 14, 34, 27, 29, 30 });
				mCases[32] = EdgeSwappingCase({ 0, 1, 14, 34, 7, 31, 30 });
				mCases[33] = EdgeSwappingCase({ 0, 1, 14, 16, 35, 28, 6 });
				mCases[34] = EdgeSwappingCase({ 0, 1, 14, 16, 35, 29, 30 });
				mCases[35] = EdgeSwappingCase({ 0, 1, 14, 16, 18, 36, 30 });
				mCases[36] = EdgeSwappingCase({ 0, 1, 14, 7, 20, 36, 30 });
				mCases[37] = EdgeSwappingCase({ 0, 1, 9, 21, 35, 28, 6 });
				mCases[38] = EdgeSwappingCase({ 0, 1, 9, 21, 35, 29, 30 });
				mCases[39] = EdgeSwappingCase({ 0, 1, 9, 21, 18, 36, 30 });
				mCases[40] = EdgeSwappingCase({ 0, 1, 9, 11, 22, 36, 30 });
				mCases[41] = EdgeSwappingCase({ 0, 1, 7, 13, 22, 36, 30 });
				mCases[42] = EdgeSwappingCase({ 0, 37, 38, 3, 4, 5, 6 });
				mCases[43] = EdgeSwappingCase({ 0, 37, 38, 3, 4, 7, 8 });
				mCases[44] = EdgeSwappingCase({ 0, 37, 38, 3, 9, 10, 6 });
				mCases[45] = EdgeSwappingCase({ 0, 37, 38, 3, 9, 11, 12 });
				mCases[46] = EdgeSwappingCase({ 0, 37, 38, 3, 7, 13, 12 });
				mCases[47] = EdgeSwappingCase({ 0, 37, 38, 14, 15, 5, 6 });
				mCases[48] = EdgeSwappingCase({ 0, 37, 38, 14, 15, 7, 8 });
				mCases[49] = EdgeSwappingCase({ 0, 37, 38, 14, 16, 17, 6 });
				mCases[50] = EdgeSwappingCase({ 0, 37, 38, 14, 16, 18, 19 });
				mCases[51] = EdgeSwappingCase({ 0, 37, 38, 14, 7, 20, 19 });
				mCases[52] = EdgeSwappingCase({ 0, 37, 38, 9, 21, 17, 6 });
				mCases[53] = EdgeSwappingCase({ 0, 37, 38, 9, 21, 18, 19 });
				mCases[54] = EdgeSwappingCase({ 0, 37, 38, 9, 11, 22, 19 });
				mCases[55] = EdgeSwappingCase({ 0, 37, 38, 7, 13, 22, 19 });
				mCases[56] = EdgeSwappingCase({ 0, 37, 39, 40, 4, 5, 6 });
				mCases[57] = EdgeSwappingCase({ 0, 37, 39, 40, 4, 7, 8 });
				mCases[58] = EdgeSwappingCase({ 0, 37, 39, 40, 9, 10, 6 });
				mCases[59] = EdgeSwappingCase({ 0, 37, 39, 40, 9, 11, 12 });
				mCases[60] = EdgeSwappingCase({ 0, 37, 39, 40, 7, 13, 12 });
				mCases[61] = EdgeSwappingCase({ 0, 37, 39, 41, 42, 5, 6 });
				mCases[62] = EdgeSwappingCase({ 0, 37, 39, 41, 42, 7, 8 });
				mCases[63] = EdgeSwappingCase({ 0, 37, 39, 41, 43, 44, 6 });
				mCases[64] = EdgeSwappingCase({ 0, 37, 39, 41, 43, 45, 46 });
				mCases[65] = EdgeSwappingCase({ 0, 37, 39, 41, 7, 47, 46 });
				mCases[66] = EdgeSwappingCase({ 0, 37, 39, 9, 48, 44, 6 });
				mCases[67] = EdgeSwappingCase({ 0, 37, 39, 9, 48, 45, 46 });
				mCases[68] = EdgeSwappingCase({ 0, 37, 39, 9, 11, 49, 46 });
				mCases[69] = EdgeSwappingCase({ 0, 37, 39, 7, 13, 49, 46 });
				mCases[70] = EdgeSwappingCase({ 0, 37, 14, 50, 42, 5, 6 });
				mCases[71] = EdgeSwappingCase({ 0, 37, 14, 50, 42, 7, 8 });
				mCases[72] = EdgeSwappingCase({ 0, 37, 14, 50, 43, 44, 6 });
				mCases[73] = EdgeSwappingCase({ 0, 37, 14, 50, 43, 45, 46 });
				mCases[74] = EdgeSwappingCase({ 0, 37, 14, 50, 7, 47, 46 });
				mCases[75] = EdgeSwappingCase({ 0, 37, 14, 16, 51, 44, 6 });
				mCases[76] = EdgeSwappingCase({ 0, 37, 14, 16, 51, 45, 46 });
				mCases[77] = EdgeSwappingCase({ 0, 37, 14, 16, 18, 52, 46 });
				mCases[78] = EdgeSwappingCase({ 0, 37, 14, 7, 20, 52, 46 });
				mCases[79] = EdgeSwappingCase({ 0, 37, 9, 21, 51, 44, 6 });
				mCases[80] = EdgeSwappingCase({ 0, 37, 9, 21, 51, 45, 46 });
				mCases[81] = EdgeSwappingCase({ 0, 37, 9, 21, 18, 52, 46 });
				mCases[82] = EdgeSwappingCase({ 0, 37, 9, 11, 22, 52, 46 });
				mCases[83] = EdgeSwappingCase({ 0, 37, 7, 13, 22, 52, 46 });
				mCases[84] = EdgeSwappingCase({ 0, 23, 53, 40, 4, 5, 6 });
				mCases[85] = EdgeSwappingCase({ 0, 23, 53, 40, 4, 7, 8 });
				mCases[86] = EdgeSwappingCase({ 0, 23, 53, 40, 9, 10, 6 });
				mCases[87] = EdgeSwappingCase({ 0, 23, 53, 40, 9, 11, 12 });
				mCases[88] = EdgeSwappingCase({ 0, 23, 53, 40, 7, 13, 12 });
				mCases[89] = EdgeSwappingCase({ 0, 23, 53, 41, 42, 5, 6 });
				mCases[90] = EdgeSwappingCase({ 0, 23, 53, 41, 42, 7, 8 });
				mCases[91] = EdgeSwappingCase({ 0, 23, 53, 41, 43, 44, 6 });
				mCases[92] = EdgeSwappingCase({ 0, 23, 53, 41, 43, 45, 46 });
				mCases[93] = EdgeSwappingCase({ 0, 23, 53, 41, 7, 47, 46 });
				mCases[94] = EdgeSwappingCase({ 0, 23, 53, 9, 48, 44, 6 });
				mCases[95] = EdgeSwappingCase({ 0, 23, 53, 9, 48, 45, 46 });
				mCases[96] = EdgeSwappingCase({ 0, 23, 53, 9, 11, 49, 46 });
				mCases[97] = EdgeSwappingCase({ 0, 23, 53, 7, 13, 49, 46 });
				mCases[98] = EdgeSwappingCase({ 0, 23, 25, 54, 42, 5, 6 });
				mCases[99] = EdgeSwappingCase({ 0, 23, 25, 54, 42, 7, 8 });
				mCases[100] = EdgeSwappingCase({ 0, 23, 25, 54, 43, 44, 6 });
				mCases[101] = EdgeSwappingCase({ 0, 23, 25, 54, 43, 45, 46 });
				mCases[102] = EdgeSwappingCase({ 0, 23, 25, 54, 7, 47, 46 });
				mCases[103] = EdgeSwappingCase({ 0, 23, 25, 27, 55, 44, 6 });
				mCases[104] = EdgeSwappingCase({ 0, 23, 25, 27, 55, 45, 46 });
				mCases[105] = EdgeSwappingCase({ 0, 23, 25, 27, 29, 56, 46 });
				mCases[106] = EdgeSwappingCase({ 0, 23, 25, 7, 31, 56, 46 });
				mCases[107] = EdgeSwappingCase({ 0, 23, 9, 32, 55, 44, 6 });
				mCases[108] = EdgeSwappingCase({ 0, 23, 9, 32, 55, 45, 46 });
				mCases[109] = EdgeSwappingCase({ 0, 23, 9, 32, 29, 56, 46 });
				mCases[110] = EdgeSwappingCase({ 0, 23, 9, 11, 33, 56, 46 });
				mCases[111] = EdgeSwappingCase({ 0, 23, 7, 13, 33, 56, 46 });
				mCases[112] = EdgeSwappingCase({ 0, 14, 34, 54, 42, 5, 6 });
				mCases[113] = EdgeSwappingCase({ 0, 14, 34, 54, 42, 7, 8 });
				mCases[114] = EdgeSwappingCase({ 0, 14, 34, 54, 43, 44, 6 });
				mCases[115] = EdgeSwappingCase({ 0, 14, 34, 54, 43, 45, 46 });
				mCases[116] = EdgeSwappingCase({ 0, 14, 34, 54, 7, 47, 46 });
				mCases[117] = EdgeSwappingCase({ 0, 14, 34, 27, 55, 44, 6 });
				mCases[118] = EdgeSwappingCase({ 0, 14, 34, 27, 55, 45, 46 });
				mCases[119] = EdgeSwappingCase({ 0, 14, 34, 27, 29, 56, 46 });
				mCases[120] = EdgeSwappingCase({ 0, 14, 34, 7, 31, 56, 46 });
				mCases[121] = EdgeSwappingCase({ 0, 14, 16, 35, 55, 44, 6 });
				mCases[122] = EdgeSwappingCase({ 0, 14, 16, 35, 55, 45, 46 });
				mCases[123] = EdgeSwappingCase({ 0, 14, 16, 35, 29, 56, 46 });
				mCases[124] = EdgeSwappingCase({ 0, 14, 16, 18, 36, 56, 46 });
				mCases[125] = EdgeSwappingCase({ 0, 14, 7, 20, 36, 56, 46 });
				mCases[126] = EdgeSwappingCase({ 0, 9, 21, 35, 55, 44, 6 });
				mCases[127] = EdgeSwappingCase({ 0, 9, 21, 35, 55, 45, 46 });
				mCases[128] = EdgeSwappingCase({ 0, 9, 21, 35, 29, 56, 46 });
				mCases[129] = EdgeSwappingCase({ 0, 9, 21, 18, 36, 56, 46 });
				mCases[130] = EdgeSwappingCase({ 0, 9, 11, 22, 36, 56, 46 });
				mCases[131] = EdgeSwappingCase({ 0, 7, 13, 22, 36, 56, 46 });
				mCases[132] = EdgeSwappingCase({ 57, 58, 2, 3, 4, 5, 6 });
				mCases[133] = EdgeSwappingCase({ 57, 58, 2, 3, 4, 7, 8 });
				mCases[134] = EdgeSwappingCase({ 57, 58, 2, 3, 9, 10, 6 });
				mCases[135] = EdgeSwappingCase({ 57, 58, 2, 3, 9, 11, 12 });
				mCases[136] = EdgeSwappingCase({ 57, 58, 2, 3, 7, 13, 12 });
				mCases[137] = EdgeSwappingCase({ 57, 58, 2, 14, 15, 5, 6 });
				mCases[138] = EdgeSwappingCase({ 57, 58, 2, 14, 15, 7, 8 });
				mCases[139] = EdgeSwappingCase({ 57, 58, 2, 14, 16, 17, 6 });
				mCases[140] = EdgeSwappingCase({ 57, 58, 2, 14, 16, 18, 19 });
				mCases[141] = EdgeSwappingCase({ 57, 58, 2, 14, 7, 20, 19 });
				mCases[142] = EdgeSwappingCase({ 57, 58, 2, 9, 21, 17, 6 });
				mCases[143] = EdgeSwappingCase({ 57, 58, 2, 9, 21, 18, 19 });
				mCases[144] = EdgeSwappingCase({ 57, 58, 2, 9, 11, 22, 19 });
				mCases[145] = EdgeSwappingCase({ 57, 58, 2, 7, 13, 22, 19 });
				mCases[146] = EdgeSwappingCase({ 57, 58, 23, 24, 4, 5, 6 });
				mCases[147] = EdgeSwappingCase({ 57, 58, 23, 24, 4, 7, 8 });
				mCases[148] = EdgeSwappingCase({ 57, 58, 23, 24, 9, 10, 6 });
				mCases[149] = EdgeSwappingCase({ 57, 58, 23, 24, 9, 11, 12 });
				mCases[150] = EdgeSwappingCase({ 57, 58, 23, 24, 7, 13, 12 });
				mCases[151] = EdgeSwappingCase({ 57, 58, 23, 25, 26, 5, 6 });
				mCases[152] = EdgeSwappingCase({ 57, 58, 23, 25, 26, 7, 8 });
				mCases[153] = EdgeSwappingCase({ 57, 58, 23, 25, 27, 28, 6 });
				mCases[154] = EdgeSwappingCase({ 57, 58, 23, 25, 27, 29, 30 });
				mCases[155] = EdgeSwappingCase({ 57, 58, 23, 25, 7, 31, 30 });
				mCases[156] = EdgeSwappingCase({ 57, 58, 23, 9, 32, 28, 6 });
				mCases[157] = EdgeSwappingCase({ 57, 58, 23, 9, 32, 29, 30 });
				mCases[158] = EdgeSwappingCase({ 57, 58, 23, 9, 11, 33, 30 });
				mCases[159] = EdgeSwappingCase({ 57, 58, 23, 7, 13, 33, 30 });
				mCases[160] = EdgeSwappingCase({ 57, 58, 14, 34, 26, 5, 6 });
				mCases[161] = EdgeSwappingCase({ 57, 58, 14, 34, 26, 7, 8 });
				mCases[162] = EdgeSwappingCase({ 57, 58, 14, 34, 27, 28, 6 });
				mCases[163] = EdgeSwappingCase({ 57, 58, 14, 34, 27, 29, 30 });
				mCases[164] = EdgeSwappingCase({ 57, 58, 14, 34, 7, 31, 30 });
				mCases[165] = EdgeSwappingCase({ 57, 58, 14, 16, 35, 28, 6 });
				mCases[166] = EdgeSwappingCase({ 57, 58, 14, 16, 35, 29, 30 });
				mCases[167] = EdgeSwappingCase({ 57, 58, 14, 16, 18, 36, 30 });
				mCases[168] = EdgeSwappingCase({ 57, 58, 14, 7, 20, 36, 30 });
				mCases[169] = EdgeSwappingCase({ 57, 58, 9, 21, 35, 28, 6 });
				mCases[170] = EdgeSwappingCase({ 57, 58, 9, 21, 35, 29, 30 });
				mCases[171] = EdgeSwappingCase({ 57, 58, 9, 21, 18, 36, 30 });
				mCases[172] = EdgeSwappingCase({ 57, 58, 9, 11, 22, 36, 30 });
				mCases[173] = EdgeSwappingCase({ 57, 58, 7, 13, 22, 36, 30 });
				mCases[174] = EdgeSwappingCase({ 57, 59, 60, 3, 4, 5, 6 });
				mCases[175] = EdgeSwappingCase({ 57, 59, 60, 3, 4, 7, 8 });
				mCases[176] = EdgeSwappingCase({ 57, 59, 60, 3, 9, 10, 6 });
				mCases[177] = EdgeSwappingCase({ 57, 59, 60, 3, 9, 11, 12 });
				mCases[178] = EdgeSwappingCase({ 57, 59, 60, 3, 7, 13, 12 });
				mCases[179] = EdgeSwappingCase({ 57, 59, 60, 14, 15, 5, 6 });
				mCases[180] = EdgeSwappingCase({ 57, 59, 60, 14, 15, 7, 8 });
				mCases[181] = EdgeSwappingCase({ 57, 59, 60, 14, 16, 17, 6 });
				mCases[182] = EdgeSwappingCase({ 57, 59, 60, 14, 16, 18, 19 });
				mCases[183] = EdgeSwappingCase({ 57, 59, 60, 14, 7, 20, 19 });
				mCases[184] = EdgeSwappingCase({ 57, 59, 60, 9, 21, 17, 6 });
				mCases[185] = EdgeSwappingCase({ 57, 59, 60, 9, 21, 18, 19 });
				mCases[186] = EdgeSwappingCase({ 57, 59, 60, 9, 11, 22, 19 });
				mCases[187] = EdgeSwappingCase({ 57, 59, 60, 7, 13, 22, 19 });
				mCases[188] = EdgeSwappingCase({ 57, 59, 61, 62, 4, 5, 6 });
				mCases[189] = EdgeSwappingCase({ 57, 59, 61, 62, 4, 7, 8 });
				mCases[190] = EdgeSwappingCase({ 57, 59, 61, 62, 9, 10, 6 });
				mCases[191] = EdgeSwappingCase({ 57, 59, 61, 62, 9, 11, 12 });
				mCases[192] = EdgeSwappingCase({ 57, 59, 61, 62, 7, 13, 12 });
				mCases[193] = EdgeSwappingCase({ 57, 59, 61, 63, 64, 5, 6 });
				mCases[194] = EdgeSwappingCase({ 57, 59, 61, 63, 64, 7, 8 });
				mCases[195] = EdgeSwappingCase({ 57, 59, 61, 63, 65, 66, 6 });
				mCases[196] = EdgeSwappingCase({ 57, 59, 61, 63, 65, 67, 68 });
				mCases[197] = EdgeSwappingCase({ 57, 59, 61, 63, 7, 69, 68 });
				mCases[198] = EdgeSwappingCase({ 57, 59, 61, 9, 70, 66, 6 });
				mCases[199] = EdgeSwappingCase({ 57, 59, 61, 9, 70, 67, 68 });
				mCases[200] = EdgeSwappingCase({ 57, 59, 61, 9, 11, 71, 68 });
				mCases[201] = EdgeSwappingCase({ 57, 59, 61, 7, 13, 71, 68 });
				mCases[202] = EdgeSwappingCase({ 57, 59, 14, 72, 64, 5, 6 });
				mCases[203] = EdgeSwappingCase({ 57, 59, 14, 72, 64, 7, 8 });
				mCases[204] = EdgeSwappingCase({ 57, 59, 14, 72, 65, 66, 6 });
				mCases[205] = EdgeSwappingCase({ 57, 59, 14, 72, 65, 67, 68 });
				mCases[206] = EdgeSwappingCase({ 57, 59, 14, 72, 7, 69, 68 });
				mCases[207] = EdgeSwappingCase({ 57, 59, 14, 16, 73, 66, 6 });
				mCases[208] = EdgeSwappingCase({ 57, 59, 14, 16, 73, 67, 68 });
				mCases[209] = EdgeSwappingCase({ 57, 59, 14, 16, 18, 74, 68 });
				mCases[210] = EdgeSwappingCase({ 57, 59, 14, 7, 20, 74, 68 });
				mCases[211] = EdgeSwappingCase({ 57, 59, 9, 21, 73, 66, 6 });
				mCases[212] = EdgeSwappingCase({ 57, 59, 9, 21, 73, 67, 68 });
				mCases[213] = EdgeSwappingCase({ 57, 59, 9, 21, 18, 74, 68 });
				mCases[214] = EdgeSwappingCase({ 57, 59, 9, 11, 22, 74, 68 });
				mCases[215] = EdgeSwappingCase({ 57, 59, 7, 13, 22, 74, 68 });
				mCases[216] = EdgeSwappingCase({ 57, 23, 75, 62, 4, 5, 6 });
				mCases[217] = EdgeSwappingCase({ 57, 23, 75, 62, 4, 7, 8 });
				mCases[218] = EdgeSwappingCase({ 57, 23, 75, 62, 9, 10, 6 });
				mCases[219] = EdgeSwappingCase({ 57, 23, 75, 62, 9, 11, 12 });
				mCases[220] = EdgeSwappingCase({ 57, 23, 75, 62, 7, 13, 12 });
				mCases[221] = EdgeSwappingCase({ 57, 23, 75, 63, 64, 5, 6 });
				mCases[222] = EdgeSwappingCase({ 57, 23, 75, 63, 64, 7, 8 });
				mCases[223] = EdgeSwappingCase({ 57, 23, 75, 63, 65, 66, 6 });
				mCases[224] = EdgeSwappingCase({ 57, 23, 75, 63, 65, 67, 68 });
				mCases[225] = EdgeSwappingCase({ 57, 23, 75, 63, 7, 69, 68 });
				mCases[226] = EdgeSwappingCase({ 57, 23, 75, 9, 70, 66, 6 });
				mCases[227] = EdgeSwappingCase({ 57, 23, 75, 9, 70, 67, 68 });
				mCases[228] = EdgeSwappingCase({ 57, 23, 75, 9, 11, 71, 68 });
				mCases[229] = EdgeSwappingCase({ 57, 23, 75, 7, 13, 71, 68 });
				mCases[230] = EdgeSwappingCase({ 57, 23, 25, 76, 64, 5, 6 });
				mCases[231] = EdgeSwappingCase({ 57, 23, 25, 76, 64, 7, 8 });
				mCases[232] = EdgeSwappingCase({ 57, 23, 25, 76, 65, 66, 6 });
				mCases[233] = EdgeSwappingCase({ 57, 23, 25, 76, 65, 67, 68 });
				mCases[234] = EdgeSwappingCase({ 57, 23, 25, 76, 7, 69, 68 });
				mCases[235] = EdgeSwappingCase({ 57, 23, 25, 27, 77, 66, 6 });
				mCases[236] = EdgeSwappingCase({ 57, 23, 25, 27, 77, 67, 68 });
				mCases[237] = EdgeSwappingCase({ 57, 23, 25, 27, 29, 78, 68 });
				mCases[238] = EdgeSwappingCase({ 57, 23, 25, 7, 31, 78, 68 });
				mCases[239] = EdgeSwappingCase({ 57, 23, 9, 32, 77, 66, 6 });
				mCases[240] = EdgeSwappingCase({ 57, 23, 9, 32, 77, 67, 68 });
				mCases[241] = EdgeSwappingCase({ 57, 23, 9, 32, 29, 78, 68 });
				mCases[242] = EdgeSwappingCase({ 57, 23, 9, 11, 33, 78, 68 });
				mCases[243] = EdgeSwappingCase({ 57, 23, 7, 13, 33, 78, 68 });
				mCases[244] = EdgeSwappingCase({ 57, 14, 34, 76, 64, 5, 6 });
				mCases[245] = EdgeSwappingCase({ 57, 14, 34, 76, 64, 7, 8 });
				mCases[246] = EdgeSwappingCase({ 57, 14, 34, 76, 65, 66, 6 });
				mCases[247] = EdgeSwappingCase({ 57, 14, 34, 76, 65, 67, 68 });
				mCases[248] = EdgeSwappingCase({ 57, 14, 34, 76, 7, 69, 68 });
				mCases[249] = EdgeSwappingCase({ 57, 14, 34, 27, 77, 66, 6 });
				mCases[250] = EdgeSwappingCase({ 57, 14, 34, 27, 77, 67, 68 });
				mCases[251] = EdgeSwappingCase({ 57, 14, 34, 27, 29, 78, 68 });
				mCases[252] = EdgeSwappingCase({ 57, 14, 34, 7, 31, 78, 68 });
				mCases[253] = EdgeSwappingCase({ 57, 14, 16, 35, 77, 66, 6 });
				mCases[254] = EdgeSwappingCase({ 57, 14, 16, 35, 77, 67, 68 });
				mCases[255] = EdgeSwappingCase({ 57, 14, 16, 35, 29, 78, 68 });
				mCases[256] = EdgeSwappingCase({ 57, 14, 16, 18, 36, 78, 68 });
				mCases[257] = EdgeSwappingCase({ 57, 14, 7, 20, 36, 78, 68 });
				mCases[258] = EdgeSwappingCase({ 57, 9, 21, 35, 77, 66, 6 });
				mCases[259] = EdgeSwappingCase({ 57, 9, 21, 35, 77, 67, 68 });
				mCases[260] = EdgeSwappingCase({ 57, 9, 21, 35, 29, 78, 68 });
				mCases[261] = EdgeSwappingCase({ 57, 9, 21, 18, 36, 78, 68 });
				mCases[262] = EdgeSwappingCase({ 57, 9, 11, 22, 36, 78, 68 });
				mCases[263] = EdgeSwappingCase({ 57, 7, 13, 22, 36, 78, 68 });
				mCases[264] = EdgeSwappingCase({ 37, 79, 60, 3, 4, 5, 6 });
				mCases[265] = EdgeSwappingCase({ 37, 79, 60, 3, 4, 7, 8 });
				mCases[266] = EdgeSwappingCase({ 37, 79, 60, 3, 9, 10, 6 });
				mCases[267] = EdgeSwappingCase({ 37, 79, 60, 3, 9, 11, 12 });
				mCases[268] = EdgeSwappingCase({ 37, 79, 60, 3, 7, 13, 12 });
				mCases[269] = EdgeSwappingCase({ 37, 79, 60, 14, 15, 5, 6 });
				mCases[270] = EdgeSwappingCase({ 37, 79, 60, 14, 15, 7, 8 });
				mCases[271] = EdgeSwappingCase({ 37, 79, 60, 14, 16, 17, 6 });
				mCases[272] = EdgeSwappingCase({ 37, 79, 60, 14, 16, 18, 19 });
				mCases[273] = EdgeSwappingCase({ 37, 79, 60, 14, 7, 20, 19 });
				mCases[274] = EdgeSwappingCase({ 37, 79, 60, 9, 21, 17, 6 });
				mCases[275] = EdgeSwappingCase({ 37, 79, 60, 9, 21, 18, 19 });
				mCases[276] = EdgeSwappingCase({ 37, 79, 60, 9, 11, 22, 19 });
				mCases[277] = EdgeSwappingCase({ 37, 79, 60, 7, 13, 22, 19 });
				mCases[278] = EdgeSwappingCase({ 37, 79, 61, 62, 4, 5, 6 });
				mCases[279] = EdgeSwappingCase({ 37, 79, 61, 62, 4, 7, 8 });
				mCases[280] = EdgeSwappingCase({ 37, 79, 61, 62, 9, 10, 6 });
				mCases[281] = EdgeSwappingCase({ 37, 79, 61, 62, 9, 11, 12 });
				mCases[282] = EdgeSwappingCase({ 37, 79, 61, 62, 7, 13, 12 });
				mCases[283] = EdgeSwappingCase({ 37, 79, 61, 63, 64, 5, 6 });
				mCases[284] = EdgeSwappingCase({ 37, 79, 61, 63, 64, 7, 8 });
				mCases[285] = EdgeSwappingCase({ 37, 79, 61, 63, 65, 66, 6 });
				mCases[286] = EdgeSwappingCase({ 37, 79, 61, 63, 65, 67, 68 });
				mCases[287] = EdgeSwappingCase({ 37, 79, 61, 63, 7, 69, 68 });
				mCases[288] = EdgeSwappingCase({ 37, 79, 61, 9, 70, 66, 6 });
				mCases[289] = EdgeSwappingCase({ 37, 79, 61, 9, 70, 67, 68 });
				mCases[290] = EdgeSwappingCase({ 37, 79, 61, 9, 11, 71, 68 });
				mCases[291] = EdgeSwappingCase({ 37, 79, 61, 7, 13, 71, 68 });
				mCases[292] = EdgeSwappingCase({ 37, 79, 14, 72, 64, 5, 6 });
				mCases[293] = EdgeSwappingCase({ 37, 79, 14, 72, 64, 7, 8 });
				mCases[294] = EdgeSwappingCase({ 37, 79, 14, 72, 65, 66, 6 });
				mCases[295] = EdgeSwappingCase({ 37, 79, 14, 72, 65, 67, 68 });
				mCases[296] = EdgeSwappingCase({ 37, 79, 14, 72, 7, 69, 68 });
				mCases[297] = EdgeSwappingCase({ 37, 79, 14, 16, 73, 66, 6 });
				mCases[298] = EdgeSwappingCase({ 37, 79, 14, 16, 73, 67, 68 });
				mCases[299] = EdgeSwappingCase({ 37, 79, 14, 16, 18, 74, 68 });
				mCases[300] = EdgeSwappingCase({ 37, 79, 14, 7, 20, 74, 68 });
				mCases[301] = EdgeSwappingCase({ 37, 79, 9, 21, 73, 66, 6 });
				mCases[302] = EdgeSwappingCase({ 37, 79, 9, 21, 73, 67, 68 });
				mCases[303] = EdgeSwappingCase({ 37, 79, 9, 21, 18, 74, 68 });
				mCases[304] = EdgeSwappingCase({ 37, 79, 9, 11, 22, 74, 68 });
				mCases[305] = EdgeSwappingCase({ 37, 79, 7, 13, 22, 74, 68 });
				mCases[306] = EdgeSwappingCase({ 37, 39, 80, 62, 4, 5, 6 });
				mCases[307] = EdgeSwappingCase({ 37, 39, 80, 62, 4, 7, 8 });
				mCases[308] = EdgeSwappingCase({ 37, 39, 80, 62, 9, 10, 6 });
				mCases[309] = EdgeSwappingCase({ 37, 39, 80, 62, 9, 11, 12 });
				mCases[310] = EdgeSwappingCase({ 37, 39, 80, 62, 7, 13, 12 });
				mCases[311] = EdgeSwappingCase({ 37, 39, 80, 63, 64, 5, 6 });
				mCases[312] = EdgeSwappingCase({ 37, 39, 80, 63, 64, 7, 8 });
				mCases[313] = EdgeSwappingCase({ 37, 39, 80, 63, 65, 66, 6 });
				mCases[314] = EdgeSwappingCase({ 37, 39, 80, 63, 65, 67, 68 });
				mCases[315] = EdgeSwappingCase({ 37, 39, 80, 63, 7, 69, 68 });
				mCases[316] = EdgeSwappingCase({ 37, 39, 80, 9, 70, 66, 6 });
				mCases[317] = EdgeSwappingCase({ 37, 39, 80, 9, 70, 67, 68 });
				mCases[318] = EdgeSwappingCase({ 37, 39, 80, 9, 11, 71, 68 });
				mCases[319] = EdgeSwappingCase({ 37, 39, 80, 7, 13, 71, 68 });
				mCases[320] = EdgeSwappingCase({ 37, 39, 41, 81, 64, 5, 6 });
				mCases[321] = EdgeSwappingCase({ 37, 39, 41, 81, 64, 7, 8 });
				mCases[322] = EdgeSwappingCase({ 37, 39, 41, 81, 65, 66, 6 });
				mCases[323] = EdgeSwappingCase({ 37, 39, 41, 81, 65, 67, 68 });
				mCases[324] = EdgeSwappingCase({ 37, 39, 41, 81, 7, 69, 68 });
				mCases[325] = EdgeSwappingCase({ 37, 39, 41, 43, 82, 66, 6 });
				mCases[326] = EdgeSwappingCase({ 37, 39, 41, 43, 82, 67, 68 });
				mCases[327] = EdgeSwappingCase({ 37, 39, 41, 43, 45, 83, 68 });
				mCases[328] = EdgeSwappingCase({ 37, 39, 41, 7, 47, 83, 68 });
				mCases[329] = EdgeSwappingCase({ 37, 39, 9, 48, 82, 66, 6 });
				mCases[330] = EdgeSwappingCase({ 37, 39, 9, 48, 82, 67, 68 });
				mCases[331] = EdgeSwappingCase({ 37, 39, 9, 48, 45, 83, 68 });
				mCases[332] = EdgeSwappingCase({ 37, 39, 9, 11, 49, 83, 68 });
				mCases[333] = EdgeSwappingCase({ 37, 39, 7, 13, 49, 83, 68 });
				mCases[334] = EdgeSwappingCase({ 37, 14, 50, 81, 64, 5, 6 });
				mCases[335] = EdgeSwappingCase({ 37, 14, 50, 81, 64, 7, 8 });
				mCases[336] = EdgeSwappingCase({ 37, 14, 50, 81, 65, 66, 6 });
				mCases[337] = EdgeSwappingCase({ 37, 14, 50, 81, 65, 67, 68 });
				mCases[338] = EdgeSwappingCase({ 37, 14, 50, 81, 7, 69, 68 });
				mCases[339] = EdgeSwappingCase({ 37, 14, 50, 43, 82, 66, 6 });
				mCases[340] = EdgeSwappingCase({ 37, 14, 50, 43, 82, 67, 68 });
				mCases[341] = EdgeSwappingCase({ 37, 14, 50, 43, 45, 83, 68 });
				mCases[342] = EdgeSwappingCase({ 37, 14, 50, 7, 47, 83, 68 });
				mCases[343] = EdgeSwappingCase({ 37, 14, 16, 51, 82, 66, 6 });
				mCases[344] = EdgeSwappingCase({ 37, 14, 16, 51, 82, 67, 68 });
				mCases[345] = EdgeSwappingCase({ 37, 14, 16, 51, 45, 83, 68 });
				mCases[346] = EdgeSwappingCase({ 37, 14, 16, 18, 52, 83, 68 });
				mCases[347] = EdgeSwappingCase({ 37, 14, 7, 20, 52, 83, 68 });
				mCases[348] = EdgeSwappingCase({ 37, 9, 21, 51, 82, 66, 6 });
				mCases[349] = EdgeSwappingCase({ 37, 9, 21, 51, 82, 67, 68 });
				mCases[350] = EdgeSwappingCase({ 37, 9, 21, 51, 45, 83, 68 });
				mCases[351] = EdgeSwappingCase({ 37, 9, 21, 18, 52, 83, 68 });
				mCases[352] = EdgeSwappingCase({ 37, 9, 11, 22, 52, 83, 68 });
				mCases[353] = EdgeSwappingCase({ 37, 7, 13, 22, 52, 83, 68 });
				mCases[354] = EdgeSwappingCase({ 23, 53, 80, 62, 4, 5, 6 });
				mCases[355] = EdgeSwappingCase({ 23, 53, 80, 62, 4, 7, 8 });
				mCases[356] = EdgeSwappingCase({ 23, 53, 80, 62, 9, 10, 6 });
				mCases[357] = EdgeSwappingCase({ 23, 53, 80, 62, 9, 11, 12 });
				mCases[358] = EdgeSwappingCase({ 23, 53, 80, 62, 7, 13, 12 });
				mCases[359] = EdgeSwappingCase({ 23, 53, 80, 63, 64, 5, 6 });
				mCases[360] = EdgeSwappingCase({ 23, 53, 80, 63, 64, 7, 8 });
				mCases[361] = EdgeSwappingCase({ 23, 53, 80, 63, 65, 66, 6 });
				mCases[362] = EdgeSwappingCase({ 23, 53, 80, 63, 65, 67, 68 });
				mCases[363] = EdgeSwappingCase({ 23, 53, 80, 63, 7, 69, 68 });
				mCases[364] = EdgeSwappingCase({ 23, 53, 80, 9, 70, 66, 6 });
				mCases[365] = EdgeSwappingCase({ 23, 53, 80, 9, 70, 67, 68 });
				mCases[366] = EdgeSwappingCase({ 23, 53, 80, 9, 11, 71, 68 });
				mCases[367] = EdgeSwappingCase({ 23, 53, 80, 7, 13, 71, 68 });
				mCases[368] = EdgeSwappingCase({ 23, 53, 41, 81, 64, 5, 6 });
				mCases[369] = EdgeSwappingCase({ 23, 53, 41, 81, 64, 7, 8 });
				mCases[370] = EdgeSwappingCase({ 23, 53, 41, 81, 65, 66, 6 });
				mCases[371] = EdgeSwappingCase({ 23, 53, 41, 81, 65, 67, 68 });
				mCases[372] = EdgeSwappingCase({ 23, 53, 41, 81, 7, 69, 68 });
				mCases[373] = EdgeSwappingCase({ 23, 53, 41, 43, 82, 66, 6 });
				mCases[374] = EdgeSwappingCase({ 23, 53, 41, 43, 82, 67, 68 });
				mCases[375] = EdgeSwappingCase({ 23, 53, 41, 43, 45, 83, 68 });
				mCases[376] = EdgeSwappingCase({ 23, 53, 41, 7, 47, 83, 68 });
				mCases[377] = EdgeSwappingCase({ 23, 53, 9, 48, 82, 66, 6 });
				mCases[378] = EdgeSwappingCase({ 23, 53, 9, 48, 82, 67, 68 });
				mCases[379] = EdgeSwappingCase({ 23, 53, 9, 48, 45, 83, 68 });
				mCases[380] = EdgeSwappingCase({ 23, 53, 9, 11, 49, 83, 68 });
				mCases[381] = EdgeSwappingCase({ 23, 53, 7, 13, 49, 83, 68 });
				mCases[382] = EdgeSwappingCase({ 23, 25, 54, 81, 64, 5, 6 });
				mCases[383] = EdgeSwappingCase({ 23, 25, 54, 81, 64, 7, 8 });
				mCases[384] = EdgeSwappingCase({ 23, 25, 54, 81, 65, 66, 6 });
				mCases[385] = EdgeSwappingCase({ 23, 25, 54, 81, 65, 67, 68 });
				mCases[386] = EdgeSwappingCase({ 23, 25, 54, 81, 7, 69, 68 });
				mCases[387] = EdgeSwappingCase({ 23, 25, 54, 43, 82, 66, 6 });
				mCases[388] = EdgeSwappingCase({ 23, 25, 54, 43, 82, 67, 68 });
				mCases[389] = EdgeSwappingCase({ 23, 25, 54, 43, 45, 83, 68 });
				mCases[390] = EdgeSwappingCase({ 23, 25, 54, 7, 47, 83, 68 });
				mCases[391] = EdgeSwappingCase({ 23, 25, 27, 55, 82, 66, 6 });
				mCases[392] = EdgeSwappingCase({ 23, 25, 27, 55, 82, 67, 68 });
				mCases[393] = EdgeSwappingCase({ 23, 25, 27, 55, 45, 83, 68 });
				mCases[394] = EdgeSwappingCase({ 23, 25, 27, 29, 56, 83, 68 });
				mCases[395] = EdgeSwappingCase({ 23, 25, 7, 31, 56, 83, 68 });
				mCases[396] = EdgeSwappingCase({ 23, 9, 32, 55, 82, 66, 6 });
				mCases[397] = EdgeSwappingCase({ 23, 9, 32, 55, 82, 67, 68 });
				mCases[398] = EdgeSwappingCase({ 23, 9, 32, 55, 45, 83, 68 });
				mCases[399] = EdgeSwappingCase({ 23, 9, 32, 29, 56, 83, 68 });
				mCases[400] = EdgeSwappingCase({ 23, 9, 11, 33, 56, 83, 68 });
				mCases[401] = EdgeSwappingCase({ 23, 7, 13, 33, 56, 83, 68 });
				mCases[402] = EdgeSwappingCase({ 14, 34, 54, 81, 64, 5, 6 });
				mCases[403] = EdgeSwappingCase({ 14, 34, 54, 81, 64, 7, 8 });
				mCases[404] = EdgeSwappingCase({ 14, 34, 54, 81, 65, 66, 6 });
				mCases[405] = EdgeSwappingCase({ 14, 34, 54, 81, 65, 67, 68 });
				mCases[406] = EdgeSwappingCase({ 14, 34, 54, 81, 7, 69, 68 });
				mCases[407] = EdgeSwappingCase({ 14, 34, 54, 43, 82, 66, 6 });
				mCases[408] = EdgeSwappingCase({ 14, 34, 54, 43, 82, 67, 68 });
				mCases[409] = EdgeSwappingCase({ 14, 34, 54, 43, 45, 83, 68 });
				mCases[410] = EdgeSwappingCase({ 14, 34, 54, 7, 47, 83, 68 });
				mCases[411] = EdgeSwappingCase({ 14, 34, 27, 55, 82, 66, 6 });
				mCases[412] = EdgeSwappingCase({ 14, 34, 27, 55, 82, 67, 68 });
				mCases[413] = EdgeSwappingCase({ 14, 34, 27, 55, 45, 83, 68 });
				mCases[414] = EdgeSwappingCase({ 14, 34, 27, 29, 56, 83, 68 });
				mCases[415] = EdgeSwappingCase({ 14, 34, 7, 31, 56, 83, 68 });
				mCases[416] = EdgeSwappingCase({ 14, 16, 35, 55, 82, 66, 6 });
				mCases[417] = EdgeSwappingCase({ 14, 16, 35, 55, 82, 67, 68 });
				mCases[418] = EdgeSwappingCase({ 14, 16, 35, 55, 45, 83, 68 });
				mCases[419] = EdgeSwappingCase({ 14, 16, 35, 29, 56, 83, 68 });
				mCases[420] = EdgeSwappingCase({ 14, 16, 18, 36, 56, 83, 68 });
				mCases[421] = EdgeSwappingCase({ 14, 7, 20, 36, 56, 83, 68 });
				mCases[422] = EdgeSwappingCase({ 9, 21, 35, 55, 82, 66, 6 });
				mCases[423] = EdgeSwappingCase({ 9, 21, 35, 55, 82, 67, 68 });
				mCases[424] = EdgeSwappingCase({ 9, 21, 35, 55, 45, 83, 68 });
				mCases[425] = EdgeSwappingCase({ 9, 21, 35, 29, 56, 83, 68 });
				mCases[426] = EdgeSwappingCase({ 9, 21, 18, 36, 56, 83, 68 });
				mCases[427] = EdgeSwappingCase({ 9, 11, 22, 36, 56, 83, 68 });
				mCases[428] = EdgeSwappingCase({ 7, 13, 22, 36, 56, 83, 68 });
				mTriangles[0] = { 0 , 1 , 2 };
				mTriangles[1] = { 0 , 2 , 3 };
				mTriangles[2] = { 0 , 3 , 4 };
				mTriangles[3] = { 0 , 4 , 5 };
				mTriangles[4] = { 0 , 5 , 6 };
				mTriangles[5] = { 0 , 6 , 7 };
				mTriangles[6] = { 0 , 7 , 8 };
				mTriangles[7] = { 6 , 7 , 8 };
				mTriangles[8] = { 0 , 6 , 8 };
				mTriangles[9] = { 5 , 6 , 7 };
				mTriangles[10] = { 0 , 5 , 7 };
				mTriangles[11] = { 5 , 7 , 8 };
				mTriangles[12] = { 0 , 5 , 8 };
				mTriangles[13] = { 5 , 6 , 8 };
				mTriangles[14] = { 4 , 5 , 6 };
				mTriangles[15] = { 0 , 4 , 6 };
				mTriangles[16] = { 4 , 6 , 7 };
				mTriangles[17] = { 0 , 4 , 7 };
				mTriangles[18] = { 4 , 7 , 8 };
				mTriangles[19] = { 0 , 4 , 8 };
				mTriangles[20] = { 4 , 6 , 8 };
				mTriangles[21] = { 4 , 5 , 7 };
				mTriangles[22] = { 4 , 5 , 8 };
				mTriangles[23] = { 3 , 4 , 5 };
				mTriangles[24] = { 0 , 3 , 5 };
				mTriangles[25] = { 3 , 5 , 6 };
				mTriangles[26] = { 0 , 3 , 6 };
				mTriangles[27] = { 3 , 6 , 7 };
				mTriangles[28] = { 0 , 3 , 7 };
				mTriangles[29] = { 3 , 7 , 8 };
				mTriangles[30] = { 0 , 3 , 8 };
				mTriangles[31] = { 3 , 6 , 8 };
				mTriangles[32] = { 3 , 5 , 7 };
				mTriangles[33] = { 3 , 5 , 8 };
				mTriangles[34] = { 3 , 4 , 6 };
				mTriangles[35] = { 3 , 4 , 7 };
				mTriangles[36] = { 3 , 4 , 8 };
				mTriangles[37] = { 2 , 3 , 4 };
				mTriangles[38] = { 0 , 2 , 4 };
				mTriangles[39] = { 2 , 4 , 5 };
				mTriangles[40] = { 0 , 2 , 5 };
				mTriangles[41] = { 2 , 5 , 6 };
				mTriangles[42] = { 0 , 2 , 6 };
				mTriangles[43] = { 2 , 6 , 7 };
				mTriangles[44] = { 0 , 2 , 7 };
				mTriangles[45] = { 2 , 7 , 8 };
				mTriangles[46] = { 0 , 2 , 8 };
				mTriangles[47] = { 2 , 6 , 8 };
				mTriangles[48] = { 2 , 5 , 7 };
				mTriangles[49] = { 2 , 5 , 8 };
				mTriangles[50] = { 2 , 4 , 6 };
				mTriangles[51] = { 2 , 4 , 7 };
				mTriangles[52] = { 2 , 4 , 8 };
				mTriangles[53] = { 2 , 3 , 5 };
				mTriangles[54] = { 2 , 3 , 6 };
				mTriangles[55] = { 2 , 3 , 7 };
				mTriangles[56] = { 2 , 3 , 8 };
				mTriangles[57] = { 1 , 2 , 3 };
				mTriangles[58] = { 0 , 1 , 3 };
				mTriangles[59] = { 1 , 3 , 4 };
				mTriangles[60] = { 0 , 1 , 4 };
				mTriangles[61] = { 1 , 4 , 5 };
				mTriangles[62] = { 0 , 1 , 5 };
				mTriangles[63] = { 1 , 5 , 6 };
				mTriangles[64] = { 0 , 1 , 6 };
				mTriangles[65] = { 1 , 6 , 7 };
				mTriangles[66] = { 0 , 1 , 7 };
				mTriangles[67] = { 1 , 7 , 8 };
				mTriangles[68] = { 0 , 1 , 8 };
				mTriangles[69] = { 1 , 6 , 8 };
				mTriangles[70] = { 1 , 5 , 7 };
				mTriangles[71] = { 1 , 5 , 8 };
				mTriangles[72] = { 1 , 4 , 6 };
				mTriangles[73] = { 1 , 4 , 7 };
				mTriangles[74] = { 1 , 4 , 8 };
				mTriangles[75] = { 1 , 3 , 5 };
				mTriangles[76] = { 1 , 3 , 6 };
				mTriangles[77] = { 1 , 3 , 7 };
				mTriangles[78] = { 1 , 3 , 8 };
				mTriangles[79] = { 1 , 2 , 4 };
				mTriangles[80] = { 1 , 2 , 5 };
				mTriangles[81] = { 1 , 2 , 6 };
				mTriangles[82] = { 1 , 2 , 7 };
				mTriangles[83] = { 1 , 2 , 8 };
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
				//if (edge.second.GetNumberOfShellPoints() == 3)
				//	EdgeSwapping3(edge.second);
				//if (edge.second.GetNumberOfShellPoints() == 4)
				//	EdgeSwapping<Internals::EdgeSwappingCases4>(edge.second);
				if (edge.second.GetNumberOfShellPoints() == 5)
					EdgeSwapping<Internals::EdgeSwappingCases5>(edge.second);
				//if (edge.second.GetNumberOfShellPoints() == 6)
				//	EdgeSwapping<Internals::EdgeSwappingCases6>(edge.second);
				//if (edge.second.GetNumberOfShellPoints() == 7)
				//	EdgeSwapping<Internals::EdgeSwappingCases7>(edge.second);
				//if (edge.second.GetNumberOfShellPoints() == 8)
				//	EdgeSwapping<Internals::EdgeSwappingCases8>(edge.second);
				//if (edge.second.GetNumberOfShellPoints() == 9)
				//	EdgeSwapping<Internals::EdgeSwappingCases9>(edge.second);
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
	mStatus.PrintData(rOStream);
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

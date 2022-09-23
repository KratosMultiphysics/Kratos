//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "includes/table.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include <chrono>  // for high_resolution_clock

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(BaseTable, KratosCoreFastSuite)
        {
            Table<double> table;
            for (std::size_t i = 0; i < 6; ++i)
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));

            double nearest = (table.GetNearestRow(2.1))[0];
            KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetValue(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(1.0), 2.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 15.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), -3.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(2.1), 2.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(1.0), 2.0);
            auto& r_data = table.Data();
            KRATOS_CHECK_EQUAL(r_data.size(), 6);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            for (std::size_t i = 0; i < 6; ++i){
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
            }
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(5.0), 10.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 10.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), 0.0);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            for (std::size_t i = 0; i < 6; ++i){
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i)+1.0);
            }
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 5.1);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 5.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(5.0), 11.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 11.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), 0.0);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            // Inverse filling with insert
            for (std::size_t i = 6; i > 0; --i)
                table.insert(static_cast<double>(i), 2.0 * static_cast<double>(i));
            KRATOS_CHECK_EQUAL(r_data.size(), 6);
            nearest = (table.GetNearestRow(2.1))[0];
            KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetValue(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(2.1), 2.0);
        }

        KRATOS_TEST_CASE_IN_SUITE(PerformanceTable, KratosCoreFastSuite)
        {
            Table<double> table;
            table.PushBack(24.99, 7817.74);
            table.PushBack(25.0, 7817.74);
            table.PushBack(50.0, 7810.84);
            table.PushBack(100.0, 7796.76);
            table.PushBack(200.0, 7767.52);
            table.PushBack(400.0, 7704.83);
            table.PushBack(600.0, 7636.78);
            table.PushBack(800.0, 7564.02);
            table.PushBack(1000.0, 7555.95);
            table.PushBack(1200.0, 7429.19);
            table.PushBack(1400.0, 7311.05);
            table.PushBack(1600.0, 6997.9);
            table.PushBack(1600.01, 6997.9);
            double d = 0.0;
            constexpr std::size_t nloop = 5;
            constexpr std::size_t nel = 1000000; //10M access
            constexpr std::size_t neval = 5;
            std::cout << std::endl;
            double sumd = 0;
            for(std::size_t l=0; l<nloop; ++l){
                // New One
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.GetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd = 0;
                }
                // GetValueOldInterpolation
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.GetValueOldInterpolation(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables GetValueOldInterpolation: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd = 0;
                }
                // Old One
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.OldGetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables OldGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd=0;
                }
                // New proposal
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.NewProposalGetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables NewProposalGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd=0;
                }
                std::cout << " End of loop ------------------- "<< l << std::endl;
            }


                std::cout << " ------------------- Now once more outside the loop, it is copy-pasted " << std::endl;
                // New One
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.GetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd=0;
                }
                // GetValueOldInterpolation
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.GetValueOldInterpolation(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables GetValueOldInterpolation: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd=0;
                }
                // Old One
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.OldGetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    sumd=0;
                    std::cout << i << " C++ Time of tables OldGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                }
                // New proposal
                for (std::size_t i = 0; i < neval; ++i){
                    auto tic = std::chrono::high_resolution_clock::now();
                    for (std::size_t j = 0; j < nel; ++j){
                        for (double t = 0.0; t<1625.0; t+=1.0){
                            d = table.NewProposalGetValue(t);
                            sumd += d;
                        }
                    }
                    auto toc = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = toc - tic;
                    std::cout << i << " C++ Time of tables NewProposalGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                    sumd=0;
                }
            std::cout << " ------------------- Now with EnthalpyTables " << std::endl;
            // Now with an enthalpy table
            Table<double> table2;
            table2.PushBack(293.14	, 0.0);
            table2.PushBack(293.15	, 32890);
            table2.PushBack(373.14	, 2.66078e+08);
            table2.PushBack(373.15	, 2.66112e+08);
            table2.PushBack(473.15	, 6.06791e+08);
            table2.PushBack(673.15	, 1.31472e+09);
            table2.PushBack(873.15	, 2.0582e+09);
            table2.PushBack(873.16	, 2.05824e+09);
            table2.PushBack(1073.15, 2.8187e+09);
            table2.PushBack(1273.15, 3.57781e+09);
            table2.PushBack(1373.15, 3.95678e+09);
            table2.PushBack(1387.15, 4.00973e+09);
            table2.PushBack(1390.65, 4.02292e+09);
            table2.PushBack(1393.65, 4.03416e+09);
            table2.PushBack(1395.65, 4.04161e+09);
            table2.PushBack(1397.15, 4.04716e+09);
            table2.PushBack(1399.15, 4.05446e+09);
            table2.PushBack(1403.15, 4.06888e+09);
            table2.PushBack(1413.15, 4.10484e+09);
            table2.PushBack(1423.15, 4.1407e+09);
            table2.PushBack(1433.15, 4.17648e+09);
            table2.PushBack(1443.15, 4.21219e+09);
            table2.PushBack(1453.15, 4.24783e+09);
            table2.PushBack(1467.15, 4.29761e+09);
            table2.PushBack(1471.15, 4.3118e+09);
            table2.PushBack(1473.15, 4.31888e+09);
            table2.PushBack(1478.15, 4.33652e+09);
            table2.PushBack(1873.15, 5.69989e+09);
            table2.PushBack(1873.16, 5.69992e+09);
            // New One
            for (std::size_t i = 0; i < neval; ++i){
                auto tic = std::chrono::high_resolution_clock::now();
                for (std::size_t j = 0; j < nel; ++j){
                    for (double t = 0.0; t<1625.0; t+=1.0){
                        d = table2.GetValue(t);
                        sumd += d;
                    }
                }
                auto toc = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = toc - tic;
                std::cout << i << " C++ Time of tables: " << elapsed.count() << " - "<< sumd <<std::endl;
                sumd=0;
            }
   
            // GetValueOldInterpolation
            for (std::size_t i = 0; i < neval; ++i){
                auto tic = std::chrono::high_resolution_clock::now();
                for (std::size_t j = 0; j < nel; ++j){
                    for (double t = 0.0; t<1625.0; t+=1.0){
                        d = table2.GetValueOldInterpolation(t);
                        sumd += d;
                    }
                }
                auto toc = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = toc - tic;
                std::cout << i << " C++ Time of tables GetValueOldInterpolation: " << elapsed.count() << " - "<< sumd <<std::endl;
                sumd=0;
            }
            // Old One
            for (std::size_t i = 0; i < neval; ++i){
                auto tic = std::chrono::high_resolution_clock::now();
                for (std::size_t j = 0; j < nel; ++j){
                    for (double t = 0.0; t<1625.0; t+=1.0){
                        d = table2.OldGetValue(t);
                        sumd += d;
                    }
                }
                auto toc = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = toc - tic;
                std::cout << i << " C++ Time of tables OldGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                sumd=0;
            }
            // New proposal
            for (std::size_t i = 0; i < neval; ++i){
                auto tic = std::chrono::high_resolution_clock::now();
                for (std::size_t j = 0; j < nel; ++j){
                    for (double t = 0.0; t<1625.0; t+=1.0){
                        d = table2.NewProposalGetValue(t);
                        sumd += d;
                    }
                }
                auto toc = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = toc - tic;
                std::cout << i << " C++ Time of tables NewProposalGetValue: " << elapsed.count() << " - "<< sumd <<std::endl;
                sumd=0;
            }


        }
    }
}  // namespace Kratos.

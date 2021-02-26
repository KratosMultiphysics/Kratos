


// int main() {



//     std::cout << " n threads = " << omp_get_max_threads() << std::endl;




//     std::size_t N =1e4;

//     std::vector<double> aux(N);

//     for(unsigned int i=0; i<aux.size(); ++i)

//         aux[i] = i*0.1;




//     double reference_value = 0.0;

//     for( int i=0; i<aux.size(); ++i)

//         reference_value += aux[i];

//     std::cout << "reference_value = " << reference_value << std::endl;




//     double value = 0.0;

//     #pragma omp parallel for

//     for( int i=0; i<aux.size(); ++i)

//         AtomicAdd(value,aux[i]);




//     std::cout << "value = " << value << std::endl;

//     assert( std::abs(value - reference_value)/reference_value < 1e-14);




//     return 0;

// }



// int main() {



//     std::cout << " n threads = " << omp_get_max_threads() << std::endl;



//     A a;

//     unsigned int tot=0;

//     #pragma omp parallel

//     {

//         a.GetTLS().vec.resize(omp_get_thread_num());

//     }







//     std::cout << a.GetTlsSize() << std::endl;



//     #pragma omp parallel

//     {

//         assert(a.GetTlsSize()==omp_get_thread_num());



//         #pragma omp atomic

//         tot += a.GetTLS().vec.size();

//     }



//     int reference_tot = 0;

//     for(unsigned int i=0;i<omp_get_max_threads(); ++i)

//         reference_tot+=i;

//     assert(tot==reference_tot);

//     std::cout << tot << std::endl;



//     return 0;

// }
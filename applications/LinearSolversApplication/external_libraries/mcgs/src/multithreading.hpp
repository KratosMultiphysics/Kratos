#pragma once

// --- External Includes ---
#ifdef MCGS_OPENMP
    #include <omp.h> // omp_lock_t
#endif

// --- STL Includes ---
#ifdef MCGS_OPENMP
    #include <vector> // std::vector
#else
    #include <cstddef> // std::size_t, std::ptrdiff_t
    #include <iterator> // std::random_access_iterator_tag

    namespace mcgs::detail {
        struct Dummy {};

        template <class T>
        struct DummyVector
        {
            using value_type = T;
            using pointer = value_type*;
            using reference = value_type&;
            using const_reference = const value_type&;
            using iterator = value_type*;
            using const_iterator = const value_type*;
            using size_type = std::size_t;
            using difference_type = std::ptrdiff_t;
            DummyVector() noexcept = default;
            DummyVector(std::size_t) noexcept {}
            DummyVector(std::size_t, const T&) noexcept {}
            void resize(std::size_t) noexcept {}
            void operator[](std::size_t) noexcept {}
            iterator begin() noexcept {return nullptr;}
            iterator end() noexcept {return nullptr;}
            size_type size() const noexcept {return 0;}
        }; // struct DummyVector
    } // namespace mcgs::detail
#endif

#ifdef MCGS_OPENMP
    #define MCGS_MUTEX omp_lock_t
    #define MCGS_MUTEX_ARRAY std::vector<MCGS_MUTEX>
    #define MCGS_INITIALIZE_MUTEX(rMutex) omp_init_lock(&rMutex)
    #define MCGS_ACQUIRE_MUTEX(rMutex) omp_set_lock(&rMutex)
    #define MCGS_RELEASE_MUTEX(rMutex) omp_unset_lock(&rMutex)
    #define MCGS_DEINITIALIZE_MUTEX(rMutex) omp_destroy_lock(&rMutex)
#else
    #define MCGS_MUTEX mcgs::detail::Dummy
    #define MCGS_MUTEX_ARRAY mcgs::detail::DummyVector<MCGS_MUTEX>
    #define MCGS_INITIALIZE_MUTEX(rMutex)
    #define MCGS_ACQUIRE_MUTEX(rMutex)
    #define MCGS_RELEASE_MUTEX(rMutex)
    #define MCGS_DEINITIALIZE_MUTEX(rMutex)
#endif

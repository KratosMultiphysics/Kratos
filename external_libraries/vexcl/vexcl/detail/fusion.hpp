#ifndef VEXCL_DETAIL_FUSION_HPP
#define VEXCL_DETAIL_FUSION_HPP

#include <boost/mpl/vector.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/sizeof.hpp>

#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/join.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/fusion/include/swap.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/zip_view.hpp>
#include <boost/fusion/include/vector_tie.hpp>
#ifndef BOOST_NO_VARIADIC_TEMPLATES
#  include <boost/fusion/adapted/std_tuple.hpp>
#endif


namespace vex {
namespace detail {

// Transform tuple of vex::vectors into mpl::vector of value_types
template <class Tuple>
struct extract_value_types {
    typedef typename std::decay<Tuple>::type T;

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4348)
#endif

    template <size_t I, size_t N, class Enable = void>
    struct loop;

    template <size_t I, size_t N>
    struct loop<I, N, typename std::enable_if<I + 1 < N>::type> {
        typedef typename boost::mpl::push_front<
            typename loop<I + 1, N>::type,
                typename std::decay<
                    typename boost::fusion::result_of::at_c<T, I>::type
                >::type::value_type
            >::type type;
    };

    template <size_t I, size_t N>
    struct loop<I, N, typename std::enable_if<I + 1 == N>::type> {
        typedef boost::mpl::vector<
            typename std::decay<
                typename boost::fusion::result_of::at_c<T, I>::type
            >::type::value_type
        > type;
    };

    typedef typename loop<0, boost::fusion::result_of::size<T>::value>::type type;

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

};

struct type_iterator {
    int pos;
    std::function<void(int, std::string)> f;

    template <class Function>
    type_iterator(Function f) : pos(0), f(f) {}

    template <class T>
    void operator()(T) {
        f(pos++, type_name<T>());
    }
};

template <class T>
void print_types(std::ostringstream &s) {
    boost::mpl::for_each<T>(
            type_iterator([&](int, std::string t) { s << "_" << t; })
            );
}

template <template<class> class Address, bool Const = false>
struct pointer_param {
    backend::source_generator &src;
    const char *name;
    int pos;

    pointer_param(backend::source_generator &src, const char *name)
        : src(src), name(name), pos(0) {}

    template <typename T>
        void operator()(T) {
            if (Const)
                src.template parameter< Address<const T> >(name + std::to_string(pos++));
            else
                src.template parameter< Address<T> >(name + std::to_string(pos++));
        }
};

template <class T>
typename std::enable_if<
    boost::fusion::traits::is_sequence<
        typename std::decay<T>::type
    >::value,
    T
>::type
forward_as_sequence(T &&t) {
    return std::forward<T>(t);
}

template <class T>
typename std::enable_if<
    !boost::fusion::traits::is_sequence<T>::value,
    boost::fusion::vector<T&>
>::type
forward_as_sequence(T &t) {
    return boost::fusion::vector<T&>(t);
}


template <class T>
typename std::enable_if<
    !boost::fusion::traits::is_sequence<T>::value,
    boost::fusion::vector<const T&>
>::type
forward_as_sequence(const T &t) {
    return boost::fusion::vector<const T&>(t);
}

//---------------------------------------------------------------------------
template <class K, size_t I = 0, class Enable = void>
struct temp_storage {
    device_vector< typename boost::mpl::at_c<K, I>::type > head;
    temp_storage<K, I + 1> tail;

    temp_storage(const backend::command_queue &queue, size_t n)
        : head(queue, n), tail(queue, n) {}

    template <size_t J>
    typename std::enable_if<
        (J > I),
        device_vector< typename boost::mpl::at_c<K, J>::type >
    >::type&
    get() {
        return tail.template get<J>();
    }

    template <size_t J>
    typename std::enable_if<
        (J == I),
        device_vector< typename boost::mpl::at_c<K, J>::type >
    >::type&
    get() {
        return head;
    }

    template <class Tuple>
    typename std::enable_if<
        boost::fusion::result_of::size<Tuple>::value == boost::mpl::size<K>::value &&
        I + 1 < static_cast<size_t>(boost::mpl::size<K>::value),
        void
    >::type
    swap(Tuple &t) {
        std::swap(head, boost::fusion::at_c<I>(t));
        tail.swap(t);
    }

    template <class Tuple>
    typename std::enable_if<
        boost::fusion::result_of::size<Tuple>::value == boost::mpl::size<K>::value &&
        I + 1 == boost::mpl::size<K>::value,
        void
    >::type
    swap(Tuple &t) {
        std::swap(head, boost::fusion::at_c<I>(t));
    }
};


template <class K, size_t I>
struct temp_storage<K, I,
    typename std::enable_if<I == boost::mpl::size<K>::value>::type >
{
    temp_storage(const backend::command_queue&, size_t) {}
};

template <size_t I, size_t N, class Enable = void>
struct arg_pusher {
    template <class T>
    static void push(backend::kernel &krn, T &&t) {
        krn.push_arg(boost::fusion::at_c<I>(t));
        arg_pusher<I + 1, N>::push(krn, std::forward<T>(t));
    }

    template <class T>
    static void push(backend::kernel &krn, temp_storage<T> &t) {
        krn.push_arg(t.template get<I>());
        arg_pusher<I + 1, N>::push(krn, t);
    }
};

template <size_t I, size_t N>
struct arg_pusher<I, N, typename std::enable_if<I == N>::type> {
    template <class T> static void push(backend::kernel&, T&&) { }
};

template <size_t N, class T>
void push_args(backend::kernel &krn, T &&t) {
    arg_pusher<0, N>::push(krn, t);
}

template <size_t I, class K>
device_vector< typename boost::mpl::at_c<K, I>::type >&
at_c(temp_storage<K, 0> &s) {
    return s.template get<I>();
}

struct extract_device_vector {
    uint d;

    extract_device_vector(uint d) : d(d) {}

    template <class T> struct result;

    template <class This, class T>
    struct result< This(T) > {
        typedef typename std::conditional<
            std::is_const<typename std::remove_reference<T>::type>::value,
            const device_vector<typename std::decay<T>::type::value_type>&,
            device_vector<typename std::decay<T>::type::value_type>&
            >::type type;
    };

    template <class T>
    typename result<extract_device_vector(T)>::type operator()(T &&t) const {
        return t(d);
    }
};

} // namespace detail
} // namespace vex

#endif

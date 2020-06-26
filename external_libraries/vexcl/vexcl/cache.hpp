#ifndef VEXCL_CACHE_HPP
#define VEXCL_CACHE_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   vexcl/cache.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Generic online cache implementation.
 */

#include <set>
#include <map>

#include <boost/thread.hpp>
#include <boost/utility.hpp>

#include <vexcl/backend.hpp>

namespace vex {
namespace detail {

// Abstract base class for object cache.
struct object_cache_base;

// List of all active caches.
template <bool dummy = true>
struct cache_register {
    static_assert(dummy, "Dummy parameter should be true");

    static std::set<object_cache_base*> caches;
    static boost::mutex caches_mx;

    static void add(object_cache_base *cache) {
        boost::lock_guard<boost::mutex> lock(caches_mx);
        caches.insert(cache);
    }

    static void remove(object_cache_base *cache) {
        boost::lock_guard<boost::mutex> lock(caches_mx);
        caches.erase(cache);
    }

    static void clear();
    static void erase(const backend::command_queue &q);
};

template <bool dummy>
std::set<object_cache_base*> cache_register<dummy>::caches;

template <bool dummy>
boost::mutex cache_register<dummy>::caches_mx;

// Abstract base class for object cache.
struct object_cache_base {
    virtual void clear() = 0;
    virtual void erase(const backend::command_queue &q) = 0;
    virtual ~object_cache_base() {}
};

template <bool dummy>
void cache_register<dummy>::clear() {
    for(auto c = caches.begin(); c != caches.end(); ++c)
        (*c)->clear();
}

template <bool dummy>
void cache_register<dummy>::erase(const backend::command_queue &q) {
    for(auto c = caches.begin(); c != caches.end(); ++c)
        (*c)->erase(q);
}

// Indexes cache objects by context
struct index_by_context {
    typedef backend::context          type;
    typedef backend::compare_contexts compare;

    static type get(const backend::command_queue &q) {
        return backend::get_context(q);
    }
};

// Indexes cache objects by command queue
struct index_by_queue {
    typedef backend::command_queue  type;
    typedef backend::compare_queues compare;

    static type get(const backend::command_queue &q) {
        return q;
    }
};

// Online cache. Stores Objects indexed by Key::type.
// Note that from the user standpoint everything is indexed by
// `const backend::command_queue&`.
template <class Key, class Object>
struct object_cache : public object_cache_base, boost::noncopyable {
    typedef std::map<typename Key::type, Object, typename Key::compare> store_type;

    store_type store;
    mutable boost::mutex store_mx;

    object_cache() {
        cache_register<true>::add(this);
    }

    ~object_cache() {
        cache_register<true>::remove(this);
    }

    template <class I>
    typename store_type::iterator
    insert(const backend::command_queue &q, I &&item) {
        boost::lock_guard<boost::mutex> lock(store_mx);

        return store.insert( std::make_pair(
                    Key::get(q), std::forward<I>(item)
                    ) ).first;
    }

    typename store_type::const_iterator end() const {
        boost::lock_guard<boost::mutex> lock(store_mx);
        return store.end();
    }

    typename store_type::iterator find(const backend::command_queue &q) {
        boost::lock_guard<boost::mutex> lock(store_mx);
        return store.find( Key::get(q) );
    }

    void clear() override {
        boost::lock_guard<boost::mutex> lock(store_mx);
        store.clear();
    }

    void erase(const backend::command_queue &q) override {
        boost::lock_guard<boost::mutex> lock(store_mx);
        store.erase( Key::get(q) );
    }
};

// The most common type of object cache is kernel cache:
typedef object_cache<index_by_context, backend::kernel> kernel_cache;

}

/// Clears cached objects, allowing to cleanly release contexts.
inline void purge_caches() {
    detail::cache_register<true>::clear();
}

/// Clears cached objects, allowing to cleanly release contexts.
inline void purge_caches(const backend::command_queue &q) {
    detail::cache_register<true>::erase(q);
}

/// Clears cached objects, allowing to cleanly release contexts.
inline void purge_caches(const std::vector<backend::command_queue> &queue) {
    for(auto q = queue.begin(); q != queue.end(); ++q)
        detail::cache_register<true>::erase( *q );
}

}


#endif

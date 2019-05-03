/******************************************************************************

The MIT License (MIT)

Copyright (c) 2017 LH_Mouse

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

******************************************************************************/

#ifndef __STD_INTRUSIVE_PTR_HPP
#define __STD_INTRUSIVE_PTR_HPP

#include <cassert> // assert
#include <cstddef> // nullptr_t
#include <atomic> // atomic<long>, atomic<T *>
#include <memory> // default_delete, unique_ptr
#include <exception> // terminate
#include <utility> // declval
#include <typeinfo> // bad_cast
#include <type_traits>
#include <mutex>

namespace std {

template<typename _T, class _D = default_delete<_T>>
class intrusive_base;

template<typename _T>
class intrusive_ptr;
template<typename _T>
class intrusive_weak_ptr;

namespace _Impl_intrusive_ptr {
	// Implementation of an atomic reference counter.
	class _Ref_count_base {
	private:
		mutable atomic<long> __x_ref;

	protected:
		constexpr _Ref_count_base() noexcept
			: __x_ref(1)
		{
		}
		constexpr _Ref_count_base(const _Ref_count_base &) noexcept
			: _Ref_count_base()
		{
		}
		_Ref_count_base &operator=(const _Ref_count_base &) noexcept {
			return *this;
		}
		~_Ref_count_base(){
			if(__x_ref.load(memory_order_relaxed) > 1){
				terminate();
			}
		}

	public:
		bool __unique() const volatile noexcept {
			return __x_ref.load(memory_order_relaxed) == 1;
		}
		long __get_ref() const volatile noexcept {
			return __x_ref.load(memory_order_relaxed);
		}
		bool __try_add_ref() const volatile noexcept {
			assert(__x_ref.load(memory_order_relaxed) >= 0);

			auto __old = __x_ref.load(memory_order_relaxed);
			do {
				if(__old == 0){
					return false;
				}
			} while(!__x_ref.compare_exchange_strong(__old, __old + 1, memory_order_relaxed));
			return true;
		}
		void __add_ref() const volatile noexcept {
			assert(__x_ref.load(memory_order_relaxed) > 0);

			__x_ref.fetch_add(1, memory_order_relaxed);
		}
		bool __drop_ref() const volatile noexcept {
			assert(__x_ref.load(memory_order_relaxed) > 0);

			return __x_ref.fetch_sub(1, memory_order_acq_rel) == 1;
		}
	};

	// Helper function template for casting from base to derived classes.
	template<typename _R, typename _S, typename = void>
	struct _Static_cast_or_dynamic_cast_helper {
		constexpr _R operator()(_S & __s) const {
			return dynamic_cast<_R>(forward<_S>(__s));
		}
	};
	template<typename _R, typename _S>
	struct _Static_cast_or_dynamic_cast_helper<_R, _S, decltype(static_cast<void>(static_cast<_R>(declval<_S>())))> {
		constexpr _R operator()(_S & __s) const {
			return static_cast<_R>(forward<_S>(__s));
		}
	};
	template<typename _R, typename _S>
	constexpr _R __static_cast_or_dynamic_cast(_S && __s){
		return _Static_cast_or_dynamic_cast_helper<_R, _S>()(__s);
	}

	// Observer for intrusive_weak_ptr.
	template<typename _T>
	class _Weak_view_template : public _Ref_count_base {
	private:
		mutable mutex __x_mutex;
		_T * __x_parent;

	public:
		explicit constexpr _Weak_view_template(_T * __parent) noexcept
			: _Ref_count_base()
			, __x_mutex(), __x_parent(__parent)
		{
		}

		_Weak_view_template(const _Weak_view_template &) = delete;
		_Weak_view_template &operator=(const _Weak_view_template &) = delete;

	public:
		bool __expired() const noexcept {
			const lock_guard<mutex> __mutex_guard(__x_mutex);
			const auto __p = __x_parent;
			if(!__p){
				return true;
			}
			return __p->_Ref_count_base::__get_ref() == 0;
		}
		template<typename _U>
		intrusive_ptr<_U> __lock() const noexcept {
			const lock_guard<mutex> __mutex_guard(__x_mutex);
			const auto __u = __static_cast_or_dynamic_cast<_U *>(__x_parent);
			if(!__u){
				return nullptr;
			}
			if(!__u->_Ref_count_base::__try_add_ref()){
				return nullptr;
			}
			return intrusive_ptr<_U>(__u);
		}
		void __unlink() noexcept {
			const lock_guard<mutex> __mutex_guard(__x_mutex);
			__x_parent = nullptr;
		}
	};
}

template<typename _T, class _D>
class intrusive_base : public _Impl_intrusive_ptr::_Ref_count_base, private _D {
	template<typename>
	friend class intrusive_ptr;
	template<typename>
	friend class intrusive_weak_ptr;

public:
	using pointer      = _T *;
	using element_type = _T;
	using deleter_type = _D;

private:
	using _Weak_view = _Impl_intrusive_ptr::_Weak_view_template<typename std::remove_cv<_T>::type>;

protected:
	template<typename _Cv_U, typename _Cv_T>
	static intrusive_ptr<_Cv_U> __fork_strong(_Cv_T * __t) noexcept {
		const auto __u = _Impl_intrusive_ptr::__static_cast_or_dynamic_cast<_Cv_U *>(__t);
		if(!__u){
			return nullptr;
		}
		__u->_Ref_count_base::__add_ref();
		return intrusive_ptr<_Cv_U>(__u);
	}
	template<typename _Cv_U, typename _Cv_T>
	static intrusive_weak_ptr<_Cv_U> __fork_weak(_Cv_T * __t){
		const auto __u = _Impl_intrusive_ptr::__static_cast_or_dynamic_cast<_Cv_U *>(__t);
		if(!__u){
			return nullptr;
		}
		return intrusive_weak_ptr<_Cv_U>(__u);
	}

private:
	mutable atomic<_Weak_view *> __x_view;

public:
	constexpr intrusive_base() noexcept
		: _Impl_intrusive_ptr::_Ref_count_base(), _D()
		, __x_view(nullptr)
	{
	}
	constexpr intrusive_base(const intrusive_base &) noexcept
		: intrusive_base()
	{
	}
	intrusive_base & operator=(const intrusive_base &) noexcept {
		return *this;
	}
	~intrusive_base(){
		const auto __v = __x_view.load(memory_order_consume);
		if(__v){
			if(__v->_Impl_intrusive_ptr::_Ref_count_base::__drop_ref()){
				delete __v;
			} else {
				__v->__unlink();
			}
		}
	}

private:
	_Weak_view * __get_weak_view() const volatile noexcept {
		auto __v = __x_view.load(memory_order_consume);
		return __v;
	}
	_Weak_view * __require_weak_view() const volatile {
		auto __v = __x_view.load(memory_order_consume);
		if(!__v){
			const auto __t = _Impl_intrusive_ptr::__static_cast_or_dynamic_cast<const volatile _T *>(this);
			if(!__t){
				throw bad_cast();
			}
			const auto __new_v = new _Weak_view(const_cast<_T *>(__t));
			if(__x_view.compare_exchange_strong(__v, __new_v, memory_order_release, memory_order_consume)){
				__v = __new_v;
			} else {
				delete __new_v;
			}
		}
		return __v;
	}

public:
	const deleter_type & get_deleter() const noexcept {
		return *this;
	}
	deleter_type & get_deleter() noexcept {
		return *this;
	}

	bool unique() const volatile noexcept {
		return _Impl_intrusive_ptr::_Ref_count_base::__unique();
	}
	long use_count() const volatile noexcept {
		return _Impl_intrusive_ptr::_Ref_count_base::__get_ref();
	}
	long weak_count() const volatile noexcept {
		const auto __v = __get_weak_view();
		if(!__v){
			return 0;
		}
		return __v->_Impl_intrusive_ptr::_Ref_count_base::__get_ref() - 1;
	}
	// Reserve the weak observer so any further construction of intrusive_weak_ptr's cannot fail.
	void reserve_weak() const volatile {
		__require_weak_view();
	}

	template<typename _U = _T>
	intrusive_ptr<const volatile _U> shared_from_this() const volatile noexcept {
		return __fork_strong<const volatile _U>(this);
	}
	template<typename _U = _T>
	intrusive_ptr<const _U> shared_from_this() const noexcept {
		return __fork_strong<const _U>(this);
	}
	template<typename _U = _T>
	intrusive_ptr<volatile _U> shared_from_this() volatile noexcept {
		return __fork_strong<volatile _U>(this);
	}
	template<typename _U = _T>
	intrusive_ptr<_U> shared_from_this() noexcept {
		return __fork_strong<_U>(this);
	}

	template<typename _U = _T>
	intrusive_weak_ptr<const volatile _U> weak_from_this() const volatile {
		return __fork_weak<const volatile _U>(this);
	}
	template<typename _U = _T>
	intrusive_weak_ptr<const _U> weak_from_this() const {
		return __fork_weak<const _U>(this);
	}
	template<typename _U = _T>
	intrusive_weak_ptr<volatile _U> weak_from_this() volatile {
		return __fork_weak<volatile _U>(this);
	}
	template<typename _U = _T>
	intrusive_weak_ptr<_U> weak_from_this() {
		return __fork_weak<_U>(this);
	}
};

namespace _Impl_intrusive_ptr {
	// Helper function template for casting pointers to user-defined classes to pointers to classes instantiated from intrusive_base.
	template<typename _T, class _D>
	const volatile intrusive_base<_T, _D> * __locate_intrusive_base(const volatile intrusive_base<_T, _D> * __p) noexcept {
		return __p;
	}
	template<typename _T, class _D>
	const intrusive_base<_T, _D> * __locate_intrusive_base(const intrusive_base<_T, _D> * __p) noexcept {
		return __p;
	}
	template<typename _T, class _D>
	volatile intrusive_base<_T, _D> * __locate_intrusive_base(volatile intrusive_base<_T, _D> * __p) noexcept {
		return __p;
	}
	template<typename _T, class _D>
	intrusive_base<_T, _D> * __locate_intrusive_base(intrusive_base<_T, _D> * __p) noexcept {
		return __p;
	}
}

template<typename _T>
class intrusive_ptr {
	template<typename>
	friend class intrusive_ptr;
	template<typename>
	friend class intrusive_weak_ptr;

public:
	using pointer      = _T *;
	using element_type = _T;
	using deleter_type = typename remove_pointer<decltype(_Impl_intrusive_ptr::__locate_intrusive_base(declval<_T *>()))>::type::deleter_type;

private:
	_T * __x_t;

private:
	_T * __fork() const noexcept {
		const auto __t = __x_t;
		if(__t){
			__t->_Impl_intrusive_ptr::_Ref_count_base::__add_ref();
		}
		return __t;
	}

public:
	constexpr intrusive_ptr(nullptr_t = nullptr) noexcept
		: __x_t(nullptr)
	{
	}
	explicit constexpr intrusive_ptr(_T * __t) noexcept
		: __x_t(__t)
	{
	}
	template<typename _U, typename _E,
		typename enable_if<
			is_convertible<typename unique_ptr<_U, _E>::pointer, pointer>::value && is_convertible<typename unique_ptr<_U, _E>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_ptr(unique_ptr<_U, _E> && __r) noexcept
		: __x_t(__r.release())
	{
	}
	template<typename _U,
		typename enable_if<
			is_convertible<typename intrusive_ptr<_U>::pointer, pointer>::value && is_convertible<typename intrusive_ptr<_U>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_ptr(const intrusive_ptr<_U> & __r) noexcept
		: __x_t(__r.__fork())
	{
	}
	template<typename _U,
		typename enable_if<
			is_convertible<typename intrusive_ptr<_U>::pointer, pointer>::value && is_convertible<typename intrusive_ptr<_U>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_ptr(intrusive_ptr<_U> && __r) noexcept
		: __x_t(__r.release())
	{
	}
	intrusive_ptr(const intrusive_ptr & __r) noexcept
		: __x_t(__r.__fork())
	{
	}
	intrusive_ptr(intrusive_ptr && __r) noexcept
		: __x_t(__r.release())
	{
	}
	intrusive_ptr & operator=(const intrusive_ptr & __r) noexcept {
		intrusive_ptr(__r).swap(*this);
		return *this;
	}
	intrusive_ptr & operator=(intrusive_ptr && __r) noexcept {
		intrusive_ptr(move(__r)).swap(*this);
		return *this;
	}
	~intrusive_ptr(){
		const auto __t = __x_t;
		if(__t){
			if(__t->_Impl_intrusive_ptr::_Ref_count_base::__drop_ref()){
				auto __d = move(_Impl_intrusive_ptr::__locate_intrusive_base(__t)->get_deleter());
				move(__d)(__t);
			}
		}
	}

public:
	constexpr _T * get() const noexcept {
		return __x_t;
	}
	_T * release() noexcept {
		const auto __t = __x_t;
		__x_t = nullptr;
		return __t;
	}
	bool unique() const noexcept {
		const auto __t = __x_t;
		if(!__t){
			return false;
		}
		return __t->_Impl_intrusive_ptr::_Ref_count_base::__unique();
	}
	long use_count() const noexcept {
		const auto __t = __x_t;
		if(!__t){
			return 0;
		}
		return __t->_Impl_intrusive_ptr::_Ref_count_base::__get_ref();
	}
	long weak_count() const noexcept {
		const auto __t = __x_t;
		if(!__t){
			return 0;
		}
		const auto __v = _Impl_intrusive_ptr::__locate_intrusive_base(__t)->__get_weak_view();
		if(!__v){
			return 0;
		}
		return __v->_Impl_intrusive_ptr::_Ref_count_base::__get_ref() - 1;
	}

	void reset(nullptr_t = nullptr) noexcept {
		intrusive_ptr().swap(*this);
	}
	void reset(_T * __t) noexcept {
		intrusive_ptr(__t).swap(*this);
	}

	void swap(intrusive_ptr & __r) noexcept {
		const auto __t = __x_t;
		__x_t = __r.__x_t;
		__r.__x_t = __t;
	}

public:
	explicit constexpr operator bool() const noexcept {
		return get() != nullptr;
	}

	constexpr _T & operator*() const {
		assert(get());

		return *get();
	}
	constexpr _T * operator->() const {
		assert(get());

		return get();
	}

	template<typename _U>
	constexpr bool operator==(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t == __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator==(_U * __u) const noexcept {
		return __x_t == __u;
	}
	template<typename _U>
	friend constexpr bool operator==(_U * __u, const intrusive_ptr<_T> & __r) noexcept {
		return __u == __r.__x_t;
	}

	template<typename _U>
	constexpr bool operator!=(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t != __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator!=(_U * __u) const noexcept {
		return __x_t != __u;
	}
	template<typename _U>
	friend constexpr bool operator!=(_U * __u, const intrusive_ptr<_T> & __r) noexcept {
		return __u != __r.__x_t;
	}

	template<typename _U>
	constexpr bool operator<(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t < __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator<(_U * __u) const noexcept {
		return __x_t < __u;
	}
	template<typename _U>
	friend constexpr bool operator<(_U * __u, const intrusive_ptr<_T> & __r)  noexcept {
		return __u < __r.__x_t;
	}

	template<typename _U>
	constexpr bool operator>(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t > __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator>(_U * __u) const noexcept {
		return __x_t > __u;
	}
	template<typename _U>
	friend constexpr bool operator>(_U * __u, const intrusive_ptr<_T> & __r) noexcept {
		return __u > __r.__x_t;
	}

	template<typename _U>
	constexpr bool operator<=(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t <= __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator<=(_U * __u) const noexcept {
		return __x_t <= __u;
	}
	template<typename _U>
	friend constexpr bool operator<=(_U * __u, const intrusive_ptr<_T> & __r) noexcept {
		return __u <= __r.__x_t;
	}

	template<typename _U>
	constexpr bool operator>=(const intrusive_ptr<_U> & __r) const noexcept {
		return __x_t >= __r.__x_t;
	}
	template<typename _U>
	constexpr bool operator>=(_U * __u) const noexcept {
		return __x_t >= __u;
	}
	template<typename _U>
	friend constexpr bool operator>=(_U * __u, const intrusive_ptr<_T> & __r) noexcept {
		return __u >= __r.__x_t;
	}

	friend void swap(intrusive_ptr<_T> & __l, intrusive_ptr<_T> & __r) noexcept {
		__l.swap(__r);
	}
};

template<typename _T, typename ..._Args>
intrusive_ptr<_T> make_intrusive(_Args &&... __args){
	static_assert(!is_array<_T>::value, "intrusive_ptr does not accept arrays.");
	static_assert(!is_reference<_T>::value, "intrusive_ptr does not accept references.");

	return intrusive_ptr<_T>(new _T(forward<_Args>(__args)...));
}

template<typename _U, typename _T>
intrusive_ptr<_U> static_pointer_cast(intrusive_ptr<_T> __r) noexcept {
	const auto __u = static_cast<_U>(__r.get());
	__r.release();
	return intrusive_ptr<_U>(__u);
}
template<typename _U, typename _T>
intrusive_ptr<_U> dynamic_pointer_cast(intrusive_ptr<_T> __r) noexcept {
	const auto __u = dynamic_cast<_U>(__r.get());
	if(__u){
		__r.release();
	}
	return intrusive_ptr<_U>(__u);
}
template<typename _U, typename _T>
intrusive_ptr<_U> const_pointer_cast(intrusive_ptr<_T> __r) noexcept {
	const auto __u = const_cast<_U>(__r.get());
	__r.release();
	return intrusive_ptr<_U>(__u);
}

template<typename _T>
class intrusive_weak_ptr {
	template<typename>
	friend class intrusive_ptr;
	template<typename>
	friend class intrusive_weak_ptr;

public:
	using pointer      = typename intrusive_ptr<_T>::pointer;
	using element_type = typename intrusive_ptr<_T>::element_type;
	using deleter_type = typename intrusive_ptr<_T>::deleter_type;

private:
	using _Weak_view = typename remove_pointer<decltype(_Impl_intrusive_ptr::__locate_intrusive_base(declval<pointer>()))>::type::_Weak_view;

private:
	static _Weak_view * _Create_view_from_element(const volatile _T * __t){
		if(!__t){
			return nullptr;
		}
		const auto __v = _Impl_intrusive_ptr::__locate_intrusive_base(__t)->__require_weak_view();
		__v->_Impl_intrusive_ptr::_Ref_count_base::__add_ref();
		return __v;
	}

private:
	_Weak_view * __x_v;

public:
	constexpr intrusive_weak_ptr(nullptr_t = nullptr) noexcept
		: __x_v(nullptr)
	{
	}
	explicit intrusive_weak_ptr(_T * __t)
		: __x_v(_Create_view_from_element(__t))
	{
	}
	template<typename _U,
		typename enable_if<
			is_convertible<typename intrusive_ptr<_U>::pointer, pointer>::value && is_convertible<typename intrusive_ptr<_U>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_weak_ptr(const intrusive_ptr<_U> & __r)
		: __x_v(_Create_view_from_element(__r.get()))
	{
	}
	template<typename _U,
		typename enable_if<
			is_convertible<typename intrusive_weak_ptr<_U>::pointer, pointer>::value && is_convertible<typename intrusive_weak_ptr<_U>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_weak_ptr(const intrusive_weak_ptr<_U> & __r) noexcept
		: __x_v(__r.__fork())
	{
	}
	template<typename _U,
		typename enable_if<
			is_convertible<typename intrusive_weak_ptr<_U>::pointer, pointer>::value && is_convertible<typename intrusive_weak_ptr<_U>::deleter_type, deleter_type>::value,
			int>::type = 0>
	intrusive_weak_ptr(intrusive_weak_ptr<_U> && __r) noexcept
		: __x_v(__r.__release())
	{
	}
	intrusive_weak_ptr(const intrusive_weak_ptr & __r) noexcept
		: __x_v(__r.__fork())
	{
	}
	intrusive_weak_ptr(intrusive_weak_ptr && __r) noexcept
		: __x_v(__r.__release())
	{
	}
	intrusive_weak_ptr & operator=(const intrusive_weak_ptr & __r) noexcept {
		intrusive_weak_ptr(__r).swap(*this);
		return *this;
	}
	intrusive_weak_ptr & operator=(intrusive_weak_ptr && __r) noexcept {
		intrusive_weak_ptr(move(__r)).swap(*this);
		return *this;
	}
	~intrusive_weak_ptr(){
		const auto __v = __x_v;
		if(__v){
			if(__v->_Impl_intrusive_ptr::_Ref_count_base::__drop_ref()){
				delete __v;
			}
		}
	}

private:
	_Weak_view * __fork() const noexcept {
		const auto __v = __x_v;
		if(!__v){
			return nullptr;
		}
		__v->_Impl_intrusive_ptr::_Ref_count_base::__add_ref();
		return __v;
	}
	_Weak_view * __release() noexcept {
		const auto __v = __x_v;
		__x_v = nullptr;
		return __v;
	}

public:
	bool expired() const noexcept {
		const auto __v = __x_v;
		if(!__v){
			return true;
		}
		return __v->__expired();
	}
	long weak_count() const noexcept {
		const auto __v = __x_v;
		if(!__v){
			return 0;
		}
		return __v->_Impl_intrusive_ptr::_Ref_count_base::__get_ref() - 1;
	}
	template<typename _U = _T>
	intrusive_ptr<_U> lock() const noexcept {
		const auto __v = __x_v;
		if(!__v){
			return nullptr;
		}
		return __v->template __lock<_U>();
	}

	void reset(nullptr_t = nullptr) noexcept {
		intrusive_weak_ptr().swap(*this);
	}
	void reset(_T * __t){
		intrusive_weak_ptr(__t).swap(*this);
	}

	void swap(intrusive_weak_ptr & __r) noexcept {
		const auto __v = __x_v;
		__x_v = __r.__x_v;
		__r.__x_v = __v;
	}

public:
	template<typename _U>
	constexpr bool operator==(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v == __r.__x_v;
	}
	template<typename _U>
	constexpr bool operator!=(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v != __r.__x_v;
	}
	template<typename _U>
	constexpr bool operator<(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v < __r.__x_v;
	}
	template<typename _U>
	constexpr bool operator>(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v > __r.__x_v;
	}
	template<typename _U>
	constexpr bool operator<=(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v <= __r.__x_v;
	}
	template<typename _U>
	constexpr bool operator>=(const intrusive_weak_ptr<_U> & __r) const noexcept {
		return __x_v >= __r.__x_v;
	}

	friend void swap(intrusive_weak_ptr<_T> & __l, intrusive_weak_ptr<_T> & __r) noexcept {
		__l.swap(__r);
	}
};

}

#endif

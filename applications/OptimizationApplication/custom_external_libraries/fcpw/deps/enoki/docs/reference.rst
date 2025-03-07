.. cpp:namespace:: enoki

Reference
=========

.. warning::

    This section is a work in progress.

How to read this section
------------------------

The function signatures provided in this section are simplified to aid
readability. For instance, the array addition operator

.. cpp:function:: template <typename Array> Array operator+(Array x, Array y)

hides the fact that

1. The function supports mixed arguments with builtin broadcasting. For
   instance, the following is legal

   .. code-block:: cpp

       Array<float, 4> x = ...;
       Array<Array<double, 8>, 4> y = ...;

       auto z = x + y;

   The entries of ``x`` are broadcasted and promoted in the above example, and
   the type of ``z`` is thus ``Array<Array<double, 8>, 4>``. See the section on
   :ref:`broadcasting rules <broadcasting>` for details.

2. The operator uses SFINAE (``std::enable_if``) so that it only becomes active
   when ``x`` or ``y`` are Enoki arrays.

3. The operator replicates the C++ typing rules. Recall that adding two
   ``float&`` references in standard C++ produces a ``float`` temporary (i.e.
   *not* a reference) holding the result of the operation.

   Analogously, ``z`` is of type ``Array<float, 4>`` in the example below.

   .. code-block:: cpp

       Array<float&, 4> x = ...;
       Array<float&, 4> y = ...;

       auto z = x + y;

Writing out these type transformation rules in every function definition would
make for tedious reading, hence the simplifications. If in doubt, please look
at the source code of Enoki.

Global macro definitions
------------------------

.. c:macro:: ENOKI_VERSION_MAJOR

    Integer value denoting the major version of the Enoki release.

.. c:macro:: ENOKI_VERSION_MINOR

    Integer value denoting the minor version of the Enoki release.

.. c:macro:: ENOKI_VERSION_PATCH

    Integer value denoting the patch version of the Enoki release.

.. c:macro:: ENOKI_VERSION

    Enoki version string (e.g. ``"0.1.2"``).

.. c:macro:: ENOKI_LIKELY(condition)

    Signals that the branch is almost always taken, which can be used for
    improved code layout if supported by the compiler. An example is shown
    below:

    .. code-block:: cpp

        if (ENOKI_LIKELY(x > 0)) {
            /// ....
         }

.. c:macro:: ENOKI_UNLIKELY(condition)

    Signals that the branch is rarely taken analogous to
    :cpp:func:`ENOKI_LIKELY`.

.. c:macro:: ENOKI_UNROLL

    Cross-platform mechanism for asking the compiler to unroll a loop. The
    macro should be placed before the ``for`` statement.

.. c:macro:: ENOKI_NOUNROLL

    Cross-platform mechanism for asking the compiler to *never* unroll a loop
    analogous to :cpp:func:`ENOKI_UNROLL`.

.. c:macro:: ENOKI_INLINE

    Cross-platform mechanism for asking the compiler to *always* inline a
    function. The macro should be placed in front of the function declaration.

    .. code-block:: cpp

        ENOKI_INLINE void foo() { ... }

.. c:macro:: ENOKI_NOINLINE

    Cross-platform mechanism for asking the compiler to *never* inline a
    function analogous to :cpp:func:`ENOKI_INLINE`.


Global variable definitions
---------------------------

.. cpp:var:: static constexpr bool has_avx512dq

    Specifies whether AVX512DQ instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512vl

    Specifies whether AVX512VL instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512bw

    Specifies whether AVX512BW instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512cd

    Specifies whether AVX512CD instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512pf

    Specifies whether AVX512PF instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512er

    Specifies whether AVX512ER instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512vpopcntdq

    Specifies whether AVX512VPOPCNTDQ instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx512f

    Specifies whether AVX512F instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx2

    Specifies whether AVX2 instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_avx

    Specifies whether AVX instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_fma

    Specifies whether fused multiply-add instructions are available on the
    target architecture (ARM & x86).

.. cpp:var:: static constexpr bool has_f16c

    Specifies whether F16C instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_sse42

    Specifies whether SSE 4.2 instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_x86_32

    Specifies whether the target architecture is x86, 32 bit.

.. cpp:var:: static constexpr bool has_x86_64

    Specifies whether the target architecture is x86, 64 bit.

.. cpp:var:: static constexpr bool has_arm_neon

    Specifies whether ARM NEON instructions are available on the target architecture.

.. cpp:var:: static constexpr bool has_arm_32

    Specifies whether the target architecture is a 32-bit ARM processor (armv7).

.. cpp:var:: static constexpr bool has_arm_64

    Specifies whether the target architecture is aarch64 (armv8+).

.. cpp:var:: static constexpr size_t max_packet_size

   Denotes the maximal packet size (in bytes) that can be mapped to native
   vector registers. It is equal to 64 if AVX512 is present, 32 if AVX is
   present, and 16 for machines with SSE 4.2 or ARM NEON.

.. cpp:var:: static constexpr size_t array_default_size

   Denotes the default size of Enoki arrays. Equal to ``max_packet_size / 4``.

Static arrays
-------------

.. cpp:class:: template <typename Value, size_t Size = array_default_size>
               Array : StaticArrayImpl<Value, Size, Array<Value, Size>>

    The default Enoki array class -- a generic container that stores a
    fixed-size array of an arbitrary data type similar to the standard template
    library class ``std::array``. The main distinction between the two is that
    :cpp:class:`enoki::Array` forwards all arithmetic operations (and other
    standard mathematical functions) to the contained elements.

    It has several template parameters:

    * ``typename Value``: the type of individual array entries.

    * ``size_t Size``: the number of packed array entries.

    This class is just a small wrapper that instantiates
    :cpp:class:`enoki::StaticArrayImpl` using the Curiously Recurring Template
    Pattern (CRTP). The latter provides the actual machinery that is needed to
    evaluate array expressions. See :ref:`custom-arrays` for details.

.. cpp:class:: template <typename Value, size_t Size = array_default_size>
               Packet : StaticArrayImpl<Value, Size, Array<Value, Size>>

    The ``Packet`` type is identical to :cpp:class:`enoki::Array` except for
    its :ref:`broadcasting behavior <broadcasting>`.

.. cpp:class:: template <typename Value, size_t Size, typename Derived> StaticArrayImpl

    This base class provides the core implementation of an Enoki array. It
    cannot be instantiated directly and is used via the Curiously Recurring
    Template Pattern (CRTP). See :cpp:class:`Array` and :ref:`custom-arrays`
    for details on how to create custom array types.

    .. cpp:function:: StaticArrayImpl()

        Create an unitialized array. Floating point arrays are initialized
        using ``std::numeric_limits<Value>::quiet_NaN()`` when the application
        is compiled in debug mode.

    .. cpp:function:: StaticArrayImpl(Value type)

        Broadcast a constant value to all entries of the array.

    .. cpp:function:: template<typename... Args> StaticArrayImpl(Args... args)

        Initialize the individual array entries with ``args`` (where
        ``sizeof...(args) == Size``).

    .. cpp:function:: template<typename Value2, typename Derived2> \
                      StaticArrayImpl(const StaticArrayImpl<Value2, Size, Derived2> &other)

        Initialize the array with the contents of another given array that
        potentially has a different underlying type. Enoki will perform a
        vectorized type conversion if this is supported by the target
        processor.

    .. cpp:function:: size_t size() const

        Returns the size of the array.

    .. cpp:function:: const Value& operator[](size_t index) const

        Return a reference to an array element (const version). When the
        application is compiled in debug mode, the function performs a range
        check and throws ``std::out_of_range`` in case of an out-of-range
        access. This behavior can be disabled by defining
        ``ENOKI_DISABLE_RANGE_CHECK``.

    .. cpp:function:: Value& operator[](size_t index)

        Return a reference to an array element. When the application is
        compiled in debug mode, the function performs a range check and throws
        ``std::out_of_range`` in case of an out-of-range access. This behavior
        can be disabled by defining ``ENOKI_DISABLE_RANGE_CHECK``.

    .. cpp:function:: const Value& coeff(size_t index) const

        Just like :cpp:func:`operator[]`, but without the range check (const
        version).

    .. cpp:function:: Value& coeff(size_t index)

        Just like :cpp:func:`operator[]`, but without the range check.

    .. cpp:function:: Value& x()

        Access the first component.

    .. cpp:function:: const Value& x() const

        Access the first component (const version).

    .. cpp:function:: Value& y()

        Access the second component.

    .. cpp:function:: const Value& y() const

        Access the second component (const version).

    .. cpp:function:: Value& z()

        Access the third component.

    .. cpp:function:: const Value& z() const

        Access the third component (const version).

    .. cpp:function:: Value& w()

        Access the fourth component.

    .. cpp:function:: const Value& w() const

        Access the fourth component (const version).

Memory allocation
-----------------

.. cpp:function:: void *alloc(size_t size)

    Allocates ``size`` bytes of memory that are sufficiently aligned so that
    any Enoki array can be safely stored at the returned address.

.. cpp:function:: template <typename T> T *alloc(size_t count)

    Typed convenience alias for :cpp:func:`alloc`. Allocates ``count *
    sizeof(T)`` bytes of memory that are sufficiently aligned so that any Enoki
    array can be safely stored at the returned address.

.. cpp:function:: void dealloc(void *ptr)

    Release the given memory region previously allocated by :cpp:func:`alloc`.

Memory operations
-----------------

.. cpp:function:: template <typename Array> Array load(const void *mem, mask_t<Array> mask = true)

    Loads an array of type ``Array`` from the memory address ``mem`` (which is
    assumed to be aligned on a multiple of ``alignof(Array)`` bytes). No loads
    are performed for entries whose mask bit is ``false``---instead, these
    entries are initialized with zero.

    .. warning::

        Performing an aligned load from an unaligned memory address will cause a
        general protection fault that immediately terminates the application.

    .. note::

        When the ``mask`` parameter is specified, the function implements a
        *masked load*, which is fairly slow on machines without the AVX512
        instruction set.

.. cpp:function:: template <typename Array> Array load_unaligned(const void *mem, mask_t<Array> mask = true)

    Loads an array of type ``Array`` from the memory address ``mem`` (which is
    not required to be aligned). No loads are performed for entries whose mask
    bit is ``false``---instead, these entries are initialized with zero.

    .. note::

        When the ``mask`` parameter is specified, the function implements a
        *masked load*, which is fairly slow on machines without the AVX512
        instruction set.

.. cpp:function:: template <typename Array> void store(const void *mem, Array array, mask_t<Array> mask = true)

    Stores an array of type ``Array`` at the memory address ``mem`` (which is
    assumed to be aligned on a multiple of ``alignof(Array)`` bytes). No stores
    are performed for entries whose mask bit is ``false``.

    .. warning::

        Performing an aligned storefrom an unaligned memory address will cause a
        general protection fault that immediately terminates the application.

    .. note::

        When the ``mask`` parameter is specified, the function implements a
        *masked store*, which is fairly slow on machines without the AVX512
        instruction set.

.. cpp:function:: template <typename Array> void store_unaligned(const void *mem, Array array, mask_t<Array> mask = true)

    Stores an array of type ``Array`` at the memory address ``mem`` (which is
    not required to be aligned). No stores are performed for entries whose mask
    bit is ``false``.

    .. note::

        When the ``mask`` parameter is specified, the function implements a
        *masked store*, which is fairly slow on machines without the AVX512
        instruction set.

.. cpp:function:: template <typename Array, size_t Stride = sizeof(scalar_t<Array>), \
                            typename Index> \
                  Array gather(const void *mem, Index index, mask_t<Array> mask = true)

    Loads an array of type ``Array`` using a masked gather operation. This is
    equivalent to the following scalar loop (which is mapped to efficient
    hardware instructions if supported by the target hardware).

    .. code-block:: cpp

        Array result;
        for (size_t i = 0; i < Array::Size; ++i)
            if (mask[i])
                result[i] = ((Value *) mem)[index[i]];
            else
                result[i] = Value(0);

    The ``index`` parameter must be a 32 or 64 bit integer array having the
    same number of entries. It will be interpreted as a signed array regardless
    of whether the provided array is signed or unsigned.

    The default value of the ``Stride`` parameter indicates that the data at
    ``mem`` uses a packed memory layout (i.e. a stride value of
    ``sizeof(Value)``); other values override this behavior.

.. cpp:function:: template <size_t Stride = 0, typename Array, typename Index> \
                  void scatter(const void *mem, Array array, Index index, mask_t<Array> mask = true)

    Stores an array of type ``Array`` using a scatter operation. This is
    equivalent to the following scalar loop (which is mapped to efficient
    hardware instructions if supported by the target hardware).

    .. code-block:: cpp

        for (size_t i = 0; i < Array::Size; ++i)
            if (mask[i])
                ((Value *) mem)[index[i]] = array[i];

    The ``index`` parameter must be a 32 or 64 bit integer array having the
    same number of entries. It will be interpreted as a signed array regardless
    of whether the provided array is signed or unsigned.

    The default value of the ``Stride`` parameter indicates that the data at
    ``mem`` uses a packed memory layout (i.e. a stride value of
    ``sizeof(Value)``); other values override this behavior.

.. cpp:function:: template <typename Array, bool Write = false, size_t Level = 2, \
                            size_t Stride = sizeof(scalar_t<Array>), typename Index> \
                  void prefetch(const void *mem, Index index, mask_t<Array> mask = true)

    Pre-fetches an array of type ``Array`` into the L1 or L2 cache (as
    indicated via the ``Level`` template parameter) to reduce the latency of a
    future gather or scatter operation. If ``Write = true``, the
    the associated cache line should be acquired for write access (i.e. a
    *scatter* rather than a *gather* operation).

    The ``index`` parameter must be a 32 or 64 bit integer array having the
    same number of entries. It will be interpreted as a signed array regardless
    of whether the provided array is signed or unsigned.

    If provided, the mask parameter specifies which of the pre-fetches should
    actually be performed.

    The default value of the ``Stride`` parameter indicates that the data at
    ``mem`` uses a packed memory layout (i.e. a stride value of
    ``sizeof(Value)``); other values override this behavior.

.. cpp:function:: template <typename Array, typename Index> \
                  void scatter_add(const void *mem, Array array, Index index, mask_t<Array> mask = true)

    Performs a scatter-add operation that is equivalent to the following scalar
    loop (mapped to efficient hardware instructions if supported by the target
    hardware).

    .. code-block:: cpp

        for (size_t i = 0; i < Array::Size; ++i)
            if (mask[i])
                ((Value *) mem)[index[i]] += array[i];

    The implementation avoids conflicts in case multiple indices refer to the
    same entry.

.. cpp:function:: template <typename Output, typename Input, typename Mask> \
                  size_t compress(Output output, Input input, Mask mask)

    Tightly packs the input values selected by a provided mask and writes them
    to ``output``, which must be a pointer or a structure of pointers. See the
    :ref:`advanced topics section <compression>` with regards to usage. The
    function returns ``count(mask)`` and also advances the pointer by this
    amount.

.. cpp:function:: template <typename Array, typename Index, typename Mask, typename Func, typename... Args> \
                  void transform(scalar_t<Array> *mem, Index index, Mask mask, Func func, Args&&... args)

    Transforms referenced entries at ``mem`` by the function ``func`` while
    avoiding potential conflicts. The variadic template arguments ``args`` are
    forwarded to the function. The pseudocode for this operation is

    .. code-block:: cpp

        for (size_t i = 0; i < Array::Size; ++i) {
            if (mask[i])
                func(mem[index], args...);
        }

    See the section on :ref:`the histogram problem and conflict detection
    <transform>` on how to use this function.

    .. note::

        To efficiently perform the transformation at the hardware level, the
        ``Index`` and ``Array`` types should occupy the same size. The
        implementation ensures that this is the case by performing an explicit
        cast of the index parameter to ``int_array_t<Array>``.


.. cpp:function:: template <typename Array, typename Index, typename Func, typename... Args> \
                  void transform(scalar_t<Array> *mem, Index index, Func func, Args&&... args)

    Unmasked version of :cpp:func:`transform`.

Miscellaneous initialization
----------------------------

.. cpp:function:: template <typename Array> Array empty()

    Returns an unitialized static array.

.. cpp:function:: template <typename DArray> DArray empty(size_t size)

    Allocates and returns a dynamic array of type ``DArray`` of size ``size``.
    The array contents are uninitialized.

.. cpp:function:: template <typename Array> Array zero()

    Returns a static array filled with zeros. This is analogous to writing
    ``Array(0)`` but makes it more explicit to the compiler that a specific
    efficient instruction sequence should be used for zero-initialization.

.. cpp:function:: template <typename DArray> DArray zero(size_t size)

    Allocates and returns a dynamic array of type ``DArray`` that is filled
    with zeros.

.. cpp:function:: template <typename Array> Array arange()

    Return an array initialized with an index sequence, i.e. ``0, 1, .., Array::Size-1``.

.. cpp:function:: template <typename DArray> DArray arange(size_t size)

    Allocates and returns a dynamic array of type ``DArray`` that is filled an
    index sequence ``0..size-1``.

.. cpp:function:: template <typename Array> Array linspace(scalar_t<Array> min, scalar_t<Array> max)

    Return an array initialized with linear linearly spaced entries including
    the endpoints ``min`` and ``max``.

.. cpp:function:: template <typename DArray> DArray linspace(scalar_t<DArray> min, scalar_t<DArray> max, size_t size)

    Allocates and returns a dynamic array initialized with ``size`` linear
    linearly spaced entries including the endpoints ``min`` and ``max``.

.. cpp:function:: template <typename DArray> Array<DArray, 2> meshgrid(const DArray &x, const DArray &y)

    Creates a 2D coordinate array containing all pairs of entries from the
    ``x`` and ``y`` arrays. Analogous to the ``meshgrid`` function in NumPy.

    .. code-block:: cpp

        using FloatP = Array<float>;
        using FloatX = DynamicArray<FloatP>;

        auto x = linspace<FloatX>(0.f, 1.f, 4);
        auto y = linspace<FloatX>(2.f, 3.f, 4);
        Array<FloatX, 2> grid = meshgrid(x, y);

        std::cout << grid << std::endl;

        /* Prints:

            [[0, 2],
             [0.333333, 2],
             [0.666667, 2],
             [1, 2],
             [0, 2.33333],
             [0.333333, 2.33333],
             [0.666667, 2.33333],
             [1, 2.33333],
             [0, 2.66667],
             [0.333333, 2.66667],
             [0.666667, 2.66667],
             [1, 2.66667],
             [0, 3],
             [0.333333, 3],
             [0.666667, 3],
             [1, 3]]
        */

.. cpp:function:: template <typename Predicate, typename Index> Index binary_search(scalar_t<Index> start, scalar_t<Index> end, Predicate pred)

    Perform binary search over a range given a predicate ``pred``, which
    monotonically decreases over this range (i.e. max one ``true`` -> ``false``
    transition).

    Given a (scalar) ``start`` and ``end`` index of a range, this function
    evaluates a predicate ``floor(log2(end-start) + 1)`` times with index
    values on the interval [start, end] (inclusive) to find the first index
    that no longer satisfies it. Note that the template parameter ``Index`` is
    automatically inferred from the supplied predicate. Specifically, the
    predicate takes an index or an index vector of type ``Index`` as input
    argument. In the vectorized case, each vector lane can thus evaluate a
    different predicate. When ``pred`` is ``false`` for all entries, the
    function returns ``start``, and when it is ``true`` for all cases, it
    returns ``end``.

    The following code example shows a typical use case: ``array`` contains a
    sorted list of floating point numbers, and the goal is to map floating
    point entries of ``x`` to the first index ``j`` such that ``array[j] >= x``
    (and all of this of course in parallel for each vector element).

    .. code-block:: cpp

        std::vector<float> array = ...;
        Float32P x = ...;

        UInt32P j = binary_search(
           0,
           array.size() - 1,
           [&](UInt32P index) {
               return gather(array.data(), index) < x;
           }
        );


Elementary Arithmetic Operators
-------------------------------

.. cpp:function:: template <typename Array> Array operator+(Array x, Array y)

    Binary addition operator.

.. cpp:function:: template <typename Array> Array operator-(Array x, Array y)

    Binary subtraction operator.

.. cpp:function:: template <typename Array> Array operator-(Array x)

    Unary minus operator.

.. cpp:function:: template <typename Array> Array operator*(Array x, Array y)

    Binary multiplication operator.

.. cpp:function:: template <typename Array> Array& operator*=(Array& x, Array y)

    Compound assignment operator for multiplication. Concerning ``Array`` types
    with non-commutative multiplication: the operation expands to ``x = x * y;``.

.. cpp:function:: template <typename Array> Array mulhi(Array x, Array y)

    Returns the high part of an integer multiplication. For 32-bit scalar
    input, this is e.g. equivalent to the following expression

    .. code-block:: cpp

        (int32_t) (((int64_t) x * (int64_t) y) >> 32);

.. cpp:function:: template <typename Array> Array operator/(Array x, Array y)

    Binary division operator. A special overload to multiply by the reciprocal
    when the second argument is a scalar.

    Integer division is handled specially, see :ref:`integer-division` for
    details.

.. cpp:function:: template <typename Array> Array operator|(Array x, Array y)

    Binary bitwise OR operator.

.. cpp:function:: template <typename Array> Array operator||(Array x, Array y)

    Binary logical OR operator (identical to ``operator|``, as no
    short-circuiting is supported in operator overloads).

.. cpp:function:: template <typename Array> Array operator&(Array x, Array y)

    Binary bitwise AND operator.

.. cpp:function:: template <typename Array> Array operator&&(Array x, Array y)

    Binary logical AND operator. (identical to ``operator&``, as no
    short-circuiting is supported in operator overloads).

.. cpp:function:: template <typename Array> Array andnot(Array x, Array y)

    Binary logical AND NOT operator. (identical to ``x & ~y``).

.. cpp:function:: template <typename Array> Array operator^(Array x, Array y)

    Binary bitwise XOR operator.

.. cpp:function:: template <typename Array> Array operator<<(Array x, Array y)

    Left shift operator. See also: :cpp:func:`sl` and :cpp:func:`rol`.

.. cpp:function:: template <typename Array> Array operator>>(Array x, Array y)

    Right shift operator. See also: :cpp:func:`sr` and :cpp:func:`ror`.

.. cpp:function:: template <typename Array> mask_t<Array> operator<(Array x, Array y)

    Less-than comparison operator.

.. cpp:function:: template <typename Array> mask_t<Array> operator<=(Array x, Array y)

    Less-than-or-equal comparison operator.

.. cpp:function:: template <typename Array> mask_t<Array> operator>(Array x, Array y)

    Greater-than comparison operator.

.. cpp:function:: template <typename Array> mask_t<Array> operator>=(Array x, Array y)

    Greater-than-or-equal comparison operator.

.. cpp:function:: template <typename Array> mask_t<Array> eq(Array x, Array y)

    Equality operator (vertical operation).

.. cpp:function:: template <typename Array> mask_t<Array> neq(Array x, Array y)

    Inequality operator (vertical operation).

.. cpp:function:: template <size_t Imm, typename Array> Array sl(Array x)

    Left shift by an immediate amount ``Imm``.

.. cpp:function:: template <size_t Imm, typename Array> Array sr(Array x)

    Right shift by an immediate amount ``Imm``.

.. cpp:function:: template <typename Array> Array rol(Array x, Array y)

    Left shift with rotation.

.. cpp:function:: template <typename Array> Array ror(Array x, Array y)

    Right shift with rotation.

.. cpp:function:: template <size_t Imm, typename Array> Array rol(Array x)

    Left shift with rotation by an immediate amount ``Imm``.

.. cpp:function:: template <size_t Imm, typename Array> Array ror(Array x)

    Right shift with rotation by an immediate amount ``Imm``.

.. cpp:function:: template <size_t Imm, typename Array> Array ror_array(Array x)

    Rotate the entire array by ``Imm`` entries towards the right, i.e.
    ``coeff[0]`` becomes ``coeff[Imm]``, etc.

.. cpp:function:: template <size_t Imm, typename Array> Array rol_array(Array x)

    Rotate the entire array by ``Imm`` entries towards the left, i.e.
    ``coeff[Imm]`` becomes ``coeff[0]``, etc.

.. cpp:function:: template <typename Target, typename Source> Target reinterpret_array(Source x)

    Reinterprets the bit-level representation of an array (e.g. from
    ``uint32_t`` to ``float``). See the section on :ref:`reinterpreting array
    contents <reinterpret>` for further details.

Elementary Arithmetic Functions
-------------------------------

.. cpp:function:: template <typename Array> Array rcp(Array x)

    Computes the reciprocal :math:`\frac{1}{x}`. A slightly less accurate (but
    more efficient) implementation is used when approximate mode is enabled for
    ``Array``. Relies on AVX512ER instructions if available.


.. cpp:function:: template <typename Array> Array rsqrt(Array x)

    Computes the reciprocal square root :math:`\frac{1}{\sqrt{x}}`. A slightly
    less accurate (but more efficient) implementation is used when approximate
    mode is enabled for ``Array``. Relies on AVX512ER instructions if available.

.. cpp:function:: template <typename Array> Array abs(Array x)

    Computes the absolute value :math:`|x|` (analogous to ``std::abs``).

.. cpp:function:: template <typename Array> Array max(Array x, Array y)

    Returns the maximum of :math:`x` and :math:`y` (analogous to ``std::max``).

.. cpp:function:: template <typename Array> Array min(Array x, Array y)

    Returns the minimum of :math:`x` and :math:`y` (analogous to ``std::min``).

.. cpp:function:: template <typename Array> Array sign(Array x)

    Computes the signum function :math:`\begin{cases}1,&\mathrm{if}\ x\ge 0\\-1,&\mathrm{otherwise}\end{cases}`

    Analogous to ``std::copysign(1.f, x)``.

.. cpp:function:: template <typename Array> Array copysign(Array x, Array y)

    Copies the sign of the array ``y`` to ``x`` (analogous to ``std::copysign``).

.. cpp:function:: template <typename Array> Array copysign_neg(Array x, Array y)

    Copies the sign of the array ``-y`` to ``x``.

.. cpp:function:: template <typename Array> Array mulsign(Array x, Array y)

    Efficiently multiplies ``x`` by the sign of ``y``.

.. cpp:function:: template <typename Array> Array mulsign_neg(Array x, Array y)

    Efficiently multiplies ``x`` by the sign of ``-y``.

.. cpp:function:: template <typename Array> Array sqr(Array x)

    Computes the square of :math:`x` (analogous to ``x*x``)

.. cpp:function:: template <typename Array> Array sqrt(Array x)

    Computes the square root of :math:`x` (analogous to ``std::sqrt``).

.. cpp:function:: template <typename Array> Array cbrt(Array x)

    Computes the cube root of :math:`x` (analogous to ``std::cbrt``).

.. cpp:function:: template <typename Array> Array hypot(Array x, Array y)

    Computes :math:`\sqrt{x^2+y^2}` while avoiding overflow and underflow.

.. cpp:function:: template <typename Array> Array ceil(Array x)

    Computes the ceiling of :math:`x` (analogous to ``std::ceil``).

.. cpp:function:: template <typename Array> Array floor(Array x)

    Computes the floor of :math:`x` (analogous to ``std::floor``).


.. cpp:function:: template <typename Target, typename Array> Array ceil2int(Array x)

    Computes the ceiling of :math:`x` and converts the result to an integer. If
    supported by the hardware, the combined operation is more efficient than
    the analogous expression ``Target(ceil(x))``.

.. cpp:function:: template <typename Target, typename Array> Array floor2int(Array x)

    Computes the floor of :math:`x` and converts the result to an integer. If
    supported by the hardware, the combined operation is more efficient than
    the analogous expression ``Target(floor(x))``.

.. cpp:function:: template <typename Array> Array round(Array x)

    Rounds :math:`x` to the nearest integer using Banker's rounding for
    half-way values.

    .. note::

        This is analogous to ``std::rint``, not ``std::round``.

.. cpp:function:: template <typename Array> Array trunc(Array x)

    Rounds :math:`x` towards zero (analogous to ``std::trunc``).

.. cpp:function:: template <typename Array> Array fmod(Array x, Array y)

    Computes the floating-point remainder of the division operation ``x/y``

.. cpp:function:: template <typename Array> Array fmadd(Array x, Array y, Array z)

    Performs a fused multiply-add operation if supported by the target
    hardware. Otherwise, the operation is emulated using conventional
    multiplication and addition (i.e. ``x * y + z``).

.. cpp:function:: template <typename Array> Array fnmadd(Array x, Array y, Array z)

    Performs a fused negative multiply-add operation if supported by the target
    hardware. Otherwise, the operation is emulated using conventional
    multiplication and addition (i.e. ``-x * y + z``).

.. cpp:function:: template <typename Array> Array fmsub(Array x, Array y, Array z)

    Performs a fused multiply-subtract operation if supported by the target
    hardware. Otherwise, the operation is emulated using conventional
    multiplication and subtraction (i.e. ``x * y - z``).

.. cpp:function:: template <typename Array> Array fnmsub(Array x, Array y, Array z)

    Performs a fused negative multiply-subtract operation if supported by the
    target hardware. Otherwise, the operation is emulated using conventional
    multiplication and subtraction (i.e. ``-x * y - z``).

.. cpp:function:: template <typename Array> Array fmaddsub(Array x, Array y, Array z)

    Performs a fused multiply-add and multiply-subtract operation for alternating elements.
    The pseudocode for this operation is

    .. code-block:: cpp

        Array result;
        for (size_t i = 0; i < Array::Size; ++i) {
            if (i % 2 == 0)
                result[i] = x[i] * y[i] - c[i];
            else
                result[i] = x[i] * y[i] + c[i];
        }

.. cpp:function:: template <typename Array> Array fmsubadd(Array x, Array y, Array z)

    Performs a fused multiply-add and multiply-subtract operation for alternating elements.
    The pseudocode for this operation is

    .. code-block:: cpp

        Array result;
        for (size_t i = 0; i < Array::Size; ++i) {
            if (i % 2 == 0)
                result[i] = x[i] * y[i] + c[i];
            else
                result[i] = x[i] * y[i] - c[i];
        }

.. cpp:function:: template <typename Array> Array ldexp(Array x, Array n)

    Multiplies :math:`x` by :math:`2^n`. Analogous to ``std::ldexp`` except
    that ``n`` is a floating point argument.

.. cpp:function:: template <typename Array> std::pair<Array, Array> frexp(Array x)

    Breaks the floating-point number :math:`x` into a normalized fraction and
    power of 2. Analogous to ``std::frexp`` except that both return values are
    floating point values.

.. cpp:function:: template <typename Array> Array lerp(Array a, Array b, Array t)

    Blends between the values :math:`a` and :math:`b` using the expression
    :math:`a(1-t) + t*b`.

.. cpp:function:: template <typename Array> Array clip(Array a, Array min, Array max)

    Clips :math:`a` to the specified interval :math:`[\texttt{min}, \texttt{max}]`.

Horizontal operations
---------------------

.. cpp:function:: template <typename Array> bool operator==(Array x, Array y)

    Equality operator.

    .. warning::

        Following the principle of least surprise,
        :cpp:func:`enoki::operator==` is a horizontal operations that returns a
        boolean value; a vertical alternatives named :cpp:func:`eq` is also
        available. The following pair of operations is equivalent:

        .. code-block:: cpp

            bool b1 = (f1 == f2);
            bool b2 = all(eq(f1, f2));

.. cpp:function:: template <typename Array> bool operator!=(Array x, Array y)

    .. warning::

        Following the principle of least surprise,
        :cpp:func:`enoki::operator!=` is a horizontal operations that returns a
        boolean value; a vertical alternatives named :cpp:func:`neq` is also
        available. The following pair of operations is equivalent:

        .. code-block:: cpp

            bool b1 = (f1 != f2);
            bool b2 = any(neq(f1, f2));

.. cpp:function:: template <typename Array> Array reverse(Array value)

    Returns the input array with all components reversed, i.e.

    .. code-block:: cpp


        value[Array::Size-1], .., value[0]

    For multidimensional arrays, the *outermost* dimension of ``Array`` are
    reversed. Note that this operation is currently not efficiently vectorized
    on 1D CPU arrays (though GPU and/or multi-dimensional arrays are fine).

.. cpp:function:: template <typename Array> Array psum(Array value)

    Computes the inclusive prefix sum of the components of ``value``, i.e.

    .. code-block:: cpp

        value[0], value[0] + value[1], .., value[0] + .. + value[Array::Size-1];

    For multidimensional arrays, the horizontal reduction is performed over the
    *outermost* dimension of ``Array``. Note that this operation is currently
    not efficiently vectorized on 1D CPU arrays (though GPU and/or
    multi-dimensional arrays are fine).

.. cpp:function:: template <typename Array> value_t<Array> hsum(Array value)

    Efficiently computes the horizontal sum of the components of ``value``, i.e.

    .. code-block:: cpp

        value[0] + .. + value[Array::Size-1];

    For 1D arrays, ``hsum()`` returns a scalar result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Array``, and the result is of type ``value_t<Array>``.

.. cpp:function:: template <typename Array> auto hsum_inner(Array value)

    Analogous to :cpp:func:`hsum`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Array`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Array> scalar_t<Array> hsum_nested(Array value)

    Recursive version of :cpp:func:`hsum`, which nests through all dimensions
    and always returns a scalar.

.. cpp:function:: template <typename Array> value_t<Array> hprod(Array value)

    Efficiently computes the horizontal product of the components of ``value``, i.e.

    .. code-block:: cpp

        value[0] * .. * value[Array::Size-1];

    For 1D arrays, ``hprod()`` returns a scalar result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Array``, and the result is of type ``value_t<Array>``.

.. cpp:function:: template <typename Array> auto hprod_inner(Array value)

    Analogous to :cpp:func:`hprod`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Array`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Array> scalar_t<Array> hprod_nested(Array value)

    Recursive version of :cpp:func:`hprod`, which nests through all dimensions
    and always returns a scalar.

.. cpp:function:: template <typename Array> value_t<Array> hmax(Array value)

    Efficiently computes the horizontal maximum of the components of ``value``, i.e.

    .. code-block:: cpp

        max(value[0], max(value[1], ...))

    For 1D arrays, ``hmax()`` returns a scalar result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Array``, and the result is of type ``value_t<Array>``.

.. cpp:function:: template <typename Array> auto hmax_inner(Array value)

    Analogous to :cpp:func:`hmax`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Array`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Array> scalar_t<Array> hmax_nested(Array value)

    Recursive version of :cpp:func:`hmax`, which nests through all dimensions
    and always returns a scalar.

.. cpp:function:: template <typename Array> value_t<Array> hmin(Array value)

    Efficiently computes the horizontal minimum of the components of ``value``, i.e.

    .. code-block:: cpp

        min(value[0], min(value[1], ...))

    For 1D arrays, ``hmin()`` returns a scalar result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Array``, and the result is of type ``value_t<Array>``.

.. cpp:function:: template <typename Array> auto hmin_inner(Array value)

    Analogous to :cpp:func:`hmin`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Array`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Array> scalar_t<Array> hmin_nested(Array value)

    Recursive version of :cpp:func:`hmin`, which nests through all dimensions
    and always returns a scalar.

.. cpp:function:: template <typename Mask> auto all(Mask value)

    Efficiently computes the horizontal AND (i.e. logical conjunction) of the
    components of the mask ``value``, i.e.

    .. code-block:: cpp

        value[0] & ... & value[Size-1]

    For 1D arrays, ``all()`` returns a boolean result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Mask``, and the result is of type ``mask_t<value_t<Mask>>``.

.. cpp:function:: template <typename Mask> auto all_inner(Mask value)

    Analogous to :cpp:func:`all`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Mask`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Mask> bool all_nested(Mask value)

    Recursive version of :cpp:func:`all`, which nests through all dimensions
    and always returns a boolean value.

.. cpp:function:: template <typename Mask, bool Default> bool all_or(Mask value)

    This function calls returns the `Default` template argument when `Mask` is
    a GPU array. Otherwise, it falls back :cpp:func:`all`. See the section on
    :ref:`horizontal operations on the GPU <horizontal_ops_on_gpu>` for
    details.

.. cpp:function:: template <typename Mask> auto any(Mask value)

    Efficiently computes the horizontal OR (i.e. logical disjunction) of the
    components of the mask ``value``, i.e.

    .. code-block:: cpp

        value[0] | ... | value[Size-1]

    For 1D arrays, ``any()`` returns a boolean result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Mask``, and the result is of type ``mask_t<value_t<Mask>>``.

.. cpp:function:: template <typename Mask> auto any_inner(Mask value)

    Analogous to :cpp:func:`any`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Mask`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Mask> bool any_nested(Mask value)

    Recursive version of :cpp:func:`any`, which nests through all dimensions
    and always returns a boolean value.

.. cpp:function:: template <typename Mask, bool Default> bool any_or(Mask value)

    This function calls returns the `Default` template argument when `Mask` is
    a GPU array. Otherwise, it falls back :cpp:func:`any`. See the section on
    :ref:`horizontal operations on the GPU <horizontal_ops_on_gpu>` for
    details.

.. cpp:function:: template <typename Mask> auto none(Mask value)

    Efficiently computes the negated horizontal OR of the components of the
    mask ``value``, i.e.

    .. code-block:: cpp

        ~(value[0] | ... | value[Size-1])

    For 1D arrays, ``none()`` returns a boolean result. For multidimensional
    arrays, the horizontal reduction is performed over the *outermost* dimension
    of ``Mask``, and the result is of type ``mask_t<value_t<Mask>>``.

.. cpp:function:: template <typename Mask> auto none_inner(Mask value)

    Analogous to :cpp:func:`none`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Mask`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Mask> bool none_nested(Mask value)

    Recursive version of :cpp:func:`none`, which nests through all dimensions
    and always returns a boolean value.

.. cpp:function:: template <typename Mask, bool Default> bool none_or(Mask value)

    This function calls returns the `Default` template argument when `Mask` is
    a GPU array. Otherwise, it falls back :cpp:func:`none`. See the section on
    :ref:`horizontal operations on the GPU <horizontal_ops_on_gpu>` for
    details.

.. cpp:function:: template <typename Mask> auto count(Mask value)

    Efficiently computes the number of components whose mask bits
    are turned on, i.e.

    .. code-block:: cpp

        (value[0] ? 1 : 0) + ... (value[Size - 1] ? 1 : 0)

    For 1D arrays, ``count()`` returns a result of type ``size_t``. For
    multidimensional arrays, the horizontal reduction is performed over the
    *outermost* dimension of ``Mask``, and the result is of type
    ``size_array_t<value_t<Mask>>``.

.. cpp:function:: template <typename Mask> auto count_inner(Mask value)

    Analogous to :cpp:func:`count`, exept that the horizontal reduction is
    performed over the *innermost* dimension of ``Mask`` (which is only
    relevant in the case of a multidimensional input array).

.. cpp:function:: template <typename Mask> size_t count_nested(Mask value)

    Recursive version of :cpp:func:`count`, which nests through all dimensions
    and always returns a boolean value.

.. cpp:function:: template <typename Array> value_t<Array> dot(Array value1, Array value2)

    Efficiently computes the dot products of ``value1`` and ``value2``, i.e.:

    .. code-block:: cpp

        value1[0]*value2[0] + .. + value1[Array::Size-1]*value2[Array::Size-1];

    The return value is of type ``value_t<Array>``, which is a scalar (e.g.
    ``float``) for ordinary inputs and an array for nested array inputs.

.. cpp:function:: template <typename Array> value_t<Array> norm(Array value)

    Computes the 2-norm of the input array. The return value is of type
    ``value_t<Array>``, which is a scalar (e.g. ``float``) for ordinary inputs
    and an array for nested array inputs.

.. cpp:function:: template <typename Array> value_t<Array> squared_norm(Array value)

    Computes the squared 2-norm of the input array. The return value is of type
    ``value_t<Array>``, which is a scalar (e.g. ``float``) for ordinary inputs
    and an array for nested array inputs.

.. cpp:function:: template <typename Array> Array normalize(Array value)

    Normalizes the input array by multiplying by the inverse of ``norm(value)``.

Transcendental functions
------------------------

.. _transcendental-accuracy:

Accuracy of transcendental function approximations
**************************************************

Most approximations of transcendental functions are based on routines in the
CEPHES math library. The table below provides some statistics on their absolute
and relative error.

The CEPHES approximations are only used when approximate mode is enabled;
otherwise, the functions below will invoke the corresponding non-vectorized
standard C library routines.

.. note::

    The forward trigonometric functions (*sin*, *cos*, *tan*) are optimized for
    low error on the domain :math:`|x| < 8192` and don't perform as well beyond
    this range.

Single precision
________________

.. list-table::
    :widths: 5 8 8 10 8 10
    :header-rows: 1
    :align: center

    * - Function
      - Tested domain
      - Abs. error (mean)
      - Abs. error (max)
      - Rel. error (mean)
      - Rel. error (max)
    * - :math:`\mathrm{sin}()`
      - :math:`-8192 < x < 8192`
      - :math:`1.2 \cdot 10^{-8}`
      - :math:`1.2 \cdot 10^{-7}`
      - :math:`1.9 \cdot 10^{-8}\,(0.25\,\mathrm{ulp})`
      - :math:`1.8 \cdot 10^{-6}\,(19\,\mathrm{ulp})`
    * - :math:`\mathrm{cos}()`
      - :math:`-8192 < x < 8192`
      - :math:`1.2 \cdot 10^{-8}`
      - :math:`1.2 \cdot 10^{-7}`
      - :math:`1.9 \cdot 10^{-8}\,(0.25\,\mathrm{ulp})`
      - :math:`3.1 \cdot 10^{-6}\,(47\,\mathrm{ulp})`
    * - :math:`\mathrm{tan}()`
      - :math:`-8192 < x < 8192`
      - :math:`4.7 \cdot 10^{-6}`
      - :math:`8.1 \cdot 10^{-1}`
      - :math:`3.4 \cdot 10^{-8}\,(0.42\,\mathrm{ulp})`
      - :math:`3.1 \cdot 10^{-6}\,(30\,\mathrm{ulp})`
    * - :math:`\mathrm{cot}()`
      - :math:`-8192 < x < 8192`
      - :math:`2.6 \cdot 10^{-6}`
      - :math:`0.11 \cdot 10^{1}`
      - :math:`3.5 \cdot 10^{-8}\,(0.42\,\mathrm{ulp})`
      - :math:`3.1 \cdot 10^{-6}\,(47\,\mathrm{ulp})`
    * - :math:`\mathrm{asin}()`
      - :math:`-1 < x < 1`
      - :math:`2.3 \cdot 10^{-8}`
      - :math:`1.2 \cdot 10^{-7}`
      - :math:`2.9 \cdot 10^{-8}\,(0.33\,\mathrm{ulp})`
      - :math:`2.3 \cdot 10^{-7}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{acos}()`
      - :math:`-1 < x < 1`
      - :math:`4.7 \cdot 10^{-8}`
      - :math:`2.4 \cdot 10^{-7}`
      - :math:`2.9 \cdot 10^{-8}\,(0.33\,\mathrm{ulp})`
      - :math:`1.2 \cdot 10^{-7}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{atan}()`
      - :math:`-1 < x < 1`
      - :math:`1.8 \cdot 10^{-7}`
      - :math:`6 \cdot 10^{-7}`
      - :math:`4.2 \cdot 10^{-7}\,(4.9\,\mathrm{ulp})`
      - :math:`8.2 \cdot 10^{-7}\,(12\,\mathrm{ulp})`
    * - :math:`\mathrm{sinh}()`
      - :math:`-10 < x < 10`
      - :math:`2.6 \cdot 10^{-5}`
      - :math:`2 \cdot 10^{-3}`
      - :math:`2.8 \cdot 10^{-8}\,(0.34\,\mathrm{ulp})`
      - :math:`2.7 \cdot 10^{-7}\,(3\,\mathrm{ulp})`
    * - :math:`\mathrm{cosh}()`
      - :math:`-10 < x < 10`
      - :math:`2.9 \cdot 10^{-5}`
      - :math:`2 \cdot 10^{-3}`
      - :math:`2.9 \cdot 10^{-8}\,(0.35\,\mathrm{ulp})`
      - :math:`2.5 \cdot 10^{-7}\,(4\,\mathrm{ulp})`
    * - :math:`\mathrm{tanh}()`
      - :math:`-10 < x < 10`
      - :math:`4.8 \cdot 10^{-8}`
      - :math:`4.2 \cdot 10^{-7}`
      - :math:`5 \cdot 10^{-8}\,(0.76\,\mathrm{ulp})`
      - :math:`5 \cdot 10^{-7}\,(7\,\mathrm{ulp})`
    * - :math:`\mathrm{csch}()`
      - :math:`-10 < x < 10`
      - :math:`5.7 \cdot 10^{-8}`
      - :math:`7.8 \cdot 10^{-3}`
      - :math:`4.4 \cdot 10^{-8}\,(0.54\,\mathrm{ulp})`
      - :math:`3.1 \cdot 10^{-7}\,(5\,\mathrm{ulp})`
    * - :math:`\mathrm{sech}()`
      - :math:`-10 < x < 10`
      - :math:`6.7 \cdot 10^{-9}`
      - :math:`1.8 \cdot 10^{-7}`
      - :math:`4.3 \cdot 10^{-8}\,(0.54\,\mathrm{ulp})`
      - :math:`3.2 \cdot 10^{-7}\,(4\,\mathrm{ulp})`
    * - :math:`\mathrm{coth}()`
      - :math:`-10 < x < 10`
      - :math:`1.2 \cdot 10^{-7}`
      - :math:`7.8 \cdot 10^{-3}`
      - :math:`6.9 \cdot 10^{-8}\,(0.61\,\mathrm{ulp})`
      - :math:`5.7 \cdot 10^{-7}\,(8\,\mathrm{ulp})`
    * - :math:`\mathrm{asinh}()`
      - :math:`-30 < x < 30`
      - :math:`2.8 \cdot 10^{-8}`
      - :math:`4.8 \cdot 10^{-7}`
      - :math:`1 \cdot 10^{-8}\,(0.13\,\mathrm{ulp})`
      - :math:`1.7 \cdot 10^{-7}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{acosh}()`
      - :math:`1 < x < 10`
      - :math:`2.9 \cdot 10^{-8}`
      - :math:`2.4 \cdot 10^{-7}`
      - :math:`1.5 \cdot 10^{-8}\,(0.18\,\mathrm{ulp})`
      - :math:`2.4 \cdot 10^{-7}\,(3\,\mathrm{ulp})`
    * - :math:`\mathrm{atanh}()`
      - :math:`-1 < x < 1`
      - :math:`9.9 \cdot 10^{-9}`
      - :math:`2.4 \cdot 10^{-7}`
      - :math:`1.5 \cdot 10^{-8}\,(0.18\,\mathrm{ulp})`
      - :math:`1.2 \cdot 10^{-7}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{exp}()`
      - :math:`-20 < x < 30`
      - :math:`0.72 \cdot 10^{4}`
      - :math:`0.1 \cdot 10^{7}`
      - :math:`2.4 \cdot 10^{-8}\,(0.27\,\mathrm{ulp})`
      - :math:`1.2 \cdot 10^{-7}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{log}()`
      - :math:`10^{-20} < x < 2\cdot 10^{30}`
      - :math:`9.6 \cdot 10^{-9}`
      - :math:`7.6 \cdot 10^{-6}`
      - :math:`1.4 \cdot 10^{-10}\,(0.0013\,\mathrm{ulp})`
      - :math:`1.2 \cdot 10^{-7}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{erf}()`
      - :math:`-1 < x < 1`
      - :math:`3.2 \cdot 10^{-8}`
      - :math:`1.8 \cdot 10^{-7}`
      - :math:`6.4 \cdot 10^{-8}\,(0.78\,\mathrm{ulp})`
      - :math:`3.3 \cdot 10^{-7}\,(4\,\mathrm{ulp})`
    * - :math:`\mathrm{erfc}()`
      - :math:`-1 < x < 1`
      - :math:`3.4 \cdot 10^{-8}`
      - :math:`2.4 \cdot 10^{-7}`
      - :math:`6.4 \cdot 10^{-8}\,(0.79\,\mathrm{ulp})`
      - :math:`1 \cdot 10^{-6}\,(11\,\mathrm{ulp})`

Double precision
________________

.. list-table::
    :widths: 5 8 8 10 8 10
    :header-rows: 1
    :align: center

    * - Function
      - Tested domain
      - Abs. error (mean)
      - Abs. error (max)
      - Rel. error (mean)
      - Rel. error (max)
    * - :math:`\mathrm{sin}()`
      - :math:`-8192 < x < 8192`
      - :math:`2.2 \cdot 10^{-17}`
      - :math:`2.2 \cdot 10^{-16}`
      - :math:`3.6 \cdot 10^{-17}\,(0.25\,\mathrm{ulp})`
      - :math:`3.1 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{cos}()`
      - :math:`-8192 < x < 8192`
      - :math:`2.2 \cdot 10^{-17}`
      - :math:`2.2 \cdot 10^{-16}`
      - :math:`3.6 \cdot 10^{-17}\,(0.25\,\mathrm{ulp})`
      - :math:`3 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{tan}()`
      - :math:`-8192 < x < 8192`
      - :math:`6.8 \cdot 10^{-16}`
      - :math:`1.2 \cdot 10^{-10}`
      - :math:`5.4 \cdot 10^{-17}\,(0.35\,\mathrm{ulp})`
      - :math:`4.1 \cdot 10^{-16}\,(3\,\mathrm{ulp})`
    * - :math:`\mathrm{cot}()`
      - :math:`-8192 < x < 8192`
      - :math:`4.9 \cdot 10^{-16}`
      - :math:`1.2 \cdot 10^{-10}`
      - :math:`5.5 \cdot 10^{-17}\,(0.36\,\mathrm{ulp})`
      - :math:`4.4 \cdot 10^{-16}\,(3\,\mathrm{ulp})`
    * - :math:`\mathrm{asin}()`
      - :math:`-1 < x < 1`
      - :math:`1.3 \cdot 10^{-17}`
      - :math:`2.2 \cdot 10^{-16}`
      - :math:`1.5 \cdot 10^{-17}\,(0.098\,\mathrm{ulp})`
      - :math:`2.2 \cdot 10^{-16}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{acos}()`
      - :math:`-1 < x < 1`
      - :math:`5.4 \cdot 10^{-17}`
      - :math:`4.4 \cdot 10^{-16}`
      - :math:`3.5 \cdot 10^{-17}\,(0.23\,\mathrm{ulp})`
      - :math:`2.2 \cdot 10^{-16}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{atan}()`
      - :math:`-1 < x < 1`
      - :math:`4.3 \cdot 10^{-17}`
      - :math:`3.3 \cdot 10^{-16}`
      - :math:`1 \cdot 10^{-16}\,(0.65\,\mathrm{ulp})`
      - :math:`7.1 \cdot 10^{-16}\,(5\,\mathrm{ulp})`
    * - :math:`\mathrm{sinh}()`
      - :math:`-10 < x < 10`
      - :math:`3.1 \cdot 10^{-14}`
      - :math:`1.8 \cdot 10^{-12}`
      - :math:`3.3 \cdot 10^{-17}\,(0.22\,\mathrm{ulp})`
      - :math:`4.3 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{cosh}()`
      - :math:`-10 < x < 10`
      - :math:`2.2 \cdot 10^{-14}`
      - :math:`1.8 \cdot 10^{-12}`
      - :math:`2 \cdot 10^{-17}\,(0.13\,\mathrm{ulp})`
      - :math:`2.9 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{tanh}()`
      - :math:`-10 < x < 10`
      - :math:`5.6 \cdot 10^{-17}`
      - :math:`3.3 \cdot 10^{-16}`
      - :math:`6.1 \cdot 10^{-17}\,(0.52\,\mathrm{ulp})`
      - :math:`5.5 \cdot 10^{-16}\,(3\,\mathrm{ulp})`
    * - :math:`\mathrm{csch}()`
      - :math:`-10 < x < 10`
      - :math:`4.5 \cdot 10^{-17}`
      - :math:`1.8 \cdot 10^{-12}`
      - :math:`3.3 \cdot 10^{-17}\,(0.21\,\mathrm{ulp})`
      - :math:`5.1 \cdot 10^{-16}\,(4\,\mathrm{ulp})`
    * - :math:`\mathrm{sech}()`
      - :math:`-10 < x < 10`
      - :math:`3 \cdot 10^{-18}`
      - :math:`2.2 \cdot 10^{-16}`
      - :math:`2 \cdot 10^{-17}\,(0.13\,\mathrm{ulp})`
      - :math:`4.3 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{coth}()`
      - :math:`-10 < x < 10`
      - :math:`1.2 \cdot 10^{-16}`
      - :math:`3.6 \cdot 10^{-12}`
      - :math:`6.2 \cdot 10^{-17}\,(0.3\,\mathrm{ulp})`
      - :math:`6.7 \cdot 10^{-16}\,(5\,\mathrm{ulp})`
    * - :math:`\mathrm{asinh}()`
      - :math:`-30 < x < 30`
      - :math:`5.1 \cdot 10^{-17}`
      - :math:`8.9 \cdot 10^{-16}`
      - :math:`1.9 \cdot 10^{-17}\,(0.13\,\mathrm{ulp})`
      - :math:`4.4 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{acosh}()`
      - :math:`1 < x < 10`
      - :math:`4.9 \cdot 10^{-17}`
      - :math:`4.4 \cdot 10^{-16}`
      - :math:`2.6 \cdot 10^{-17}\,(0.17\,\mathrm{ulp})`
      - :math:`6.6 \cdot 10^{-16}\,(5\,\mathrm{ulp})`
    * - :math:`\mathrm{atanh}()`
      - :math:`-1 < x < 1`
      - :math:`1.8 \cdot 10^{-17}`
      - :math:`4.4 \cdot 10^{-16}`
      - :math:`3.2 \cdot 10^{-17}\,(0.21\,\mathrm{ulp})`
      - :math:`3 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{exp}()`
      - :math:`-20 < x < 30`
      - :math:`4.7 \cdot 10^{-6}`
      - :math:`2 \cdot 10^{-3}`
      - :math:`2.5 \cdot 10^{-17}\,(0.16\,\mathrm{ulp})`
      - :math:`3.3 \cdot 10^{-16}\,(2\,\mathrm{ulp})`
    * - :math:`\mathrm{log}()`
      - :math:`10^{-20} < x < 2\cdot 10^{30}`
      - :math:`1.9 \cdot 10^{-17}`
      - :math:`1.4 \cdot 10^{-14}`
      - :math:`2.7 \cdot 10^{-19}\,(0.0013\,\mathrm{ulp})`
      - :math:`2.2 \cdot 10^{-16}\,(1\,\mathrm{ulp})`
    * - :math:`\mathrm{erf}()`
      - :math:`-1 < x < 1`
      - :math:`4.7 \cdot 10^{-17}`
      - :math:`4.4 \cdot 10^{-16}`
      - :math:`9.6 \cdot 10^{-17}\,(0.63\,\mathrm{ulp})`
      - :math:`5.9 \cdot 10^{-16}\,(5\,\mathrm{ulp})`
    * - :math:`\mathrm{erfc}()`
      - :math:`-1 < x < 1`
      - :math:`4.8 \cdot 10^{-17}`
      - :math:`4.4 \cdot 10^{-16}`
      - :math:`9.6 \cdot 10^{-17}\,(0.64\,\mathrm{ulp})`
      - :math:`2.5 \cdot 10^{-15}\,(16\,\mathrm{ulp})`

Trigonometric functions
***********************

.. cpp:function:: template <typename Array> Array sin(Array x)

    Sine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array cos(Array x)

    Cosine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> std::pair<Array, Array> sincos(Array x)

    Simultaneous sine and cosine function approximation based on the CEPHES
    library.

.. cpp:function:: template <typename Array> Array tan(Array x)

    Tangent function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array csc(Array x)

    Cosecant convenience function implemented as ``rcp(sin(x))``.

.. cpp:function:: template <typename Array> Array sec(Array x)

    Cosecant convenience function implemented as ``rcp(cos(x))``.

.. cpp:function:: template <typename Array> Array cot(Array x)

    Cotangent convenience function implemented as ``rcp(tan(x))``.

.. cpp:function:: template <typename Array> Array asin(Array x)

    Arcsine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array acos(Array x)

    Arccosine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array atan(Array x)

    Arctangent function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array atan2(Array y, Array x)

    Arctangent function of two variables.

Hyperbolic functions
********************

.. cpp:function:: template <typename Array> Array sinh(Array x)

    Hyperbolic sine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array cosh(Array x)

    Hyperbolic cosine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> std::pair<Array, Array> sincosh(Array x)

    Simultaneous hyperbolic sine and cosine function approximation based on the
    CEPHES library.

.. cpp:function:: template <typename Array> Array tanh(Array x)

    Hyperbolic tangent function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array csch(Array x)

    Hyperbolic cosecant convenience function implemented as ``rcp(sinh(x))``.

.. cpp:function:: template <typename Array> Array sech(Array x)

    Hyperbolic secant convenience function.

.. cpp:function:: template <typename Array> Array coth(Array x)

    Hyperbolic cotangent convenience function implemented as ``rcp(tanh(x))``.

.. cpp:function:: template <typename Array> Array asinh(Array x)

    Hyperbolic arcsine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array acosh(Array x)

    Hyperbolic arccosine function approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array atanh(Array x)

    Hyperbolic arctangent function approximation based on the CEPHES library.

Exponential, logarithm, and others
**********************************

.. cpp:function:: template <typename Array> Array exp(Array x)

   Natural exponential function approximation based on the CEPHES library.
   Relies on AVX512ER instructions if available.

.. cpp:function:: template <typename Array> Array log(Array x)

    Natural logarithm approximation based on the CEPHES library.

.. cpp:function:: template <typename Array> Array pow(Array x, Array y)

    Computes the power function :math:`x^y`.

"Safe" versions of mathematical functions
-----------------------------------------

.. cpp:function:: template <typename Array> Array safe_sqrt(Array x)

    Computes ``sqrt(max(0, x))`` to avoid issues with negative inputs
    (e.g. due to roundoff error in a prior calculation).

.. cpp:function:: template <typename Array> Array safe_rsqrt(Array x)

    Computes ``rsqrt(max(0, x))`` to avoid issues with negative inputs
    (e.g. due to roundoff error in a prior calculation).

.. cpp:function:: template <typename Array> Array safe_asin(Array x)

    Computes ``asin(min(1, max(-1, x)))`` to avoid issues with
    out-of-range inputs (e.g. due to roundoff error in a prior calculation).

.. cpp:function:: template <typename Array> Array safe_acos(Array x)

    Computes ``acos(min(1, max(-1, x)))`` to avoid issues with
    out-of-range inputs (e.g. due to roundoff error in a prior calculation).

Special functions
-----------------

The following special functions require including the header
:file:`enoki/special.h`.

.. cpp:function:: template <typename Array> Array erf(Array x)

    Evaluates the error function defined as

    .. math::

        \mathrm{erf}(x)=\frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}\,\mathrm{d}t.

    Requires a real-valued input array ``x``.

.. cpp:function:: template <typename Array> Array erfi(Array x)

    Evaluates the imaginary error function defined as

    .. math::

        \mathrm{erfi}(x)=-i\,\mathrm{erf}(ix).

    Requires a real-valued input array ``x``.

.. cpp:function:: template <typename Array> Array erfinv(Array x)

    Evaluates the inverse of the error function :cpp:func:`erf`.

.. cpp:function:: template <typename Array> Array i0e(Array x)

    Evaluates the exponentially scaled modified Bessel function of order zero
    defined as

    .. math::

        I_0^{(e)}(x) = e^{-|x|} I_0(x),

    where

    .. math::

        I_0(x) = \frac{1}{\pi} \int_{0}^\pi e^{x\cos \theta}\mathrm{d}\theta.

.. cpp:function:: template <typename Array> Array gamma(Array x)

    Evaluates the Gamma function defined as

    .. math::

        \Gamma(x)=\int_0^\infty t^{x-1} e^{-t}\,\mathrm{d}t.

.. cpp:function:: template <typename Array> Array lgamma(Array x)

    Evaluates the natural logarithm of the Gamma function.

.. cpp:function:: template <typename Array> Array dawson(Array x)

    Evaluates Dawson's integral defined as

    .. math::

        D(x)=e^{-x^2}\int_0^x e^{t^2}\,\mathrm{d}t.

.. cpp:function:: template <typename Array> Array ellint_1(Array phi, Array k)

    Evaluates the incomplete elliptic integral of the first kind

    .. math::

        F(\phi, k)=\int_0^\phi (1-k^2\sin^2\theta)^{-\frac{1}{2}}\,\mathrm{d}\theta

.. cpp:function:: template <typename Array> Array comp_ellint_1(Array k)

    Evaluates the complete elliptic integral of the first kind

    .. math::

        F(k)=\int_0^\frac{\pi}{2} (1-k^2\sin^2\theta)^{-\frac{1}{2}}\,\mathrm{d}\theta

.. cpp:function:: template <typename Array> Array ellint_2(Array phi, Array k)

    Evaluates the incomplete elliptic integral of the second kind

    .. math::

        E(\phi, k)=\int_0^\phi (1-k^2\sin^2\theta)^{\frac{1}{2}}\,\mathrm{d}\theta

.. cpp:function:: template <typename Array> Array comp_ellint_2(Array k)

    Evaluates the complete elliptic integral of the second kind

    .. math::

        E(k)=\int_0^\frac{\pi}{2} (1-k^2\sin^2\theta)^{\frac{1}{2}}\,\mathrm{d}\theta

.. cpp:function:: template <typename Array> Array ellint_3(Array phi, Array k, Array nu)

    Evaluates the incomplete elliptic integral of the third kind

    .. math::

        \Pi(\phi, k, \nu)=\int_0^\phi (1+\nu\sin^2\theta)^{-1}(1-k^2\sin^2\theta)^{-\frac{1}{2}}\,\mathrm{d}\theta

.. cpp:function:: template <typename Array> Array comp_ellint_3(Array k, Array nu)

    Evaluates the complete elliptic integral of the third kind

    .. math::

        \Pi(k, \nu)=\int_0^\frac{\pi}{2} (1+\nu\sin^2\theta)^{-1}(1-k^2\sin^2\theta)^{-\frac{1}{2}}\,\mathrm{d}\theta


Miscellaneous operations
------------------------

.. cpp:function:: template <typename Array> mask_t<Array> isnan(Array x)

    Checks for NaN values and returns a mask, analogous to ``std::isnan``.

.. cpp:function:: template <typename Array> mask_t<Array> isinf(Array x)

    Checks for infinite values and returns a mask, analogous to ``std::isinf``.

.. cpp:function:: template <typename Array> mask_t<Array> isfinite(Array x)

    Checks for finite values and returns a mask, analogous to ``std::isfinite``.

.. cpp:function:: template <typename Array> mask_t<Array> isdenormal(Array x)

    Checks for denormalized values and returns a mask.

.. cpp:function:: template <typename Array> Array deg_to_rad(Array array)

    Convenience function which multiplies the input array by :math:`\frac{\pi}{180}`.

.. cpp:function:: template <typename Array> Array rad_to_deg(Array array)

    Convenience function which multiplies the input array by :math:`\frac{180}{\pi}`.

.. cpp:function:: template <typename Array> Array prev_float(Array array)

    Return the prev representable floating point value for each element of
    ``array`` analogous to ``std::nextafter(array, -∞)``. Special values
    (infinities & not-a-number values) are returned unchanged.

.. cpp:function:: template <typename Array> Array next_float(Array array)

    Return the next representable floating point value for each element of
    ``array`` analogous to ``std::nextafter(array, ∞)``. Special values
    (infinities & not-a-number values) are returned unchanged.

.. cpp:function:: template <typename Array> Array tzcnt(Array array)

    Return the number of trailing zero bits (assumes that ``Array`` is an integer array).

.. cpp:function:: template <typename Array> Array lzcnt(Array array)

    Return the number of leading zero bits (assumes that ``Array`` is an integer array).

.. cpp:function:: template <typename Array> Array popcnt(Array array)

    Return the number nonzero bits (assumes that ``Array`` is an integer array).

.. cpp:function:: template <typename Array> Array log2i(Array array)

    Return the floor of the base-two logarithm (assumes that ``Array`` is an integer array).

.. cpp:function:: template <typename Index> std::pair<Index, mask_t<Index>> range(scalar_t<Index> size)

    Returns an iterable, which generates linearly increasing index vectors from
    ``0`` to ``size-1``. This function is meant to be used with the C++11
    range-based for loop:

    .. code-block:: cpp

        for (auto pair : range<Index>(1000)) {
            Index index = pair.first;
            mask_t<Index> mask = pair.second;

            // ...
        }

    The mask specifies which index vector entries are active: unless the number
    of interations is exactly divisible by the packet size, the last loop
    iteration will generally have several disabled entries.

    The implementation also supports iteration of multidimensional arrays

    .. code-block:: cpp

        using UInt32P = Packet<uint32_t>;
        using Vector3uP = Array<UInt32P, 3>;

        for (auto pair : range<Vector3uP>(10, 20, 30)) {
            Vector3uP index = pair.first;
            mask_t<Index> mask = pair.second;

            // ...
        }

    which will generate indices ``(0, 0, 0)``, ``(0, 0, 1)``, etc. As before, the
    last loop iteration will generally have several disabled entries.


.. cpp:function:: void set_flush_denormals(bool value)

    Arithmetic involving denormalized floating point numbers triggers a `slow
    microcode handler <https://en.wikipedia.org/wiki/Denormal_number#Performance_issues>`_
    on most current architectures, which leads to severe performance penalties.
    This function can be used to specify whether denormalized floating point
    values are simply flushed to zero, which sidesteps the performance issues.

    Enoki also provides a tiny a RAII wrapper named `scoped_flush_denormals`
    which sets (and later resets) this parameter.

.. cpp:function:: bool flush_denormals()

    Returns the denormals are flushed to zero (see :cpp:func:`set_flush_denormals`).

.. cpp:function:: template <typename Array> auto unit_angle(const Array &a, const Array &b)

    Numerically well-behaved routine for computing the angle between two unit
    direction vectors. This should be used wherever one is tempted to compute
    the arc cosine of a dot product, i.e. ``acos(dot(a, b))``.

    Proposed by `Don Hatch <http://www.plunk.org/~hatch/rightway.php>`_.

.. cpp:function:: template <typename Array> auto unit_angle_z(const Array &v)

    Numerically well-behaved routine for computing the angle between the unit
    vector ``v`` and the Z axis ``[0, 0, 1]``. This should be used wherever one
    is tempted to compute the arc cosine, i.e. ``acos(v.z())``.

    Proposed by `Don Hatch <http://www.plunk.org/~hatch/rightway.php>`_.


Rearranging contents of arrays
------------------------------

.. cpp:function:: template <size_t... Index, typename Array> shuffle(Array a)

    Shuffles the contents of an array given indices specified as a template
    parameter. The pseudocode for this operation is

    .. code-block:: cpp

        Array out;
        for (size_t i = 0; i<Array::Size; ++i)
            out[i] = a[Index[i]];
        return out;

.. cpp:function:: template <typename Array, typename Indices> shuffle(Array a, Indices ind)

    Shuffles the contents of an array given indices specified via an integer
    array. The pseudocode for this operation is

    .. code-block:: cpp

        Array out;
        for (size_t i = 0; i<Array::Size; ++i)
            out[i] = a[ind[i]];
        return out;

.. cpp:function:: template <typename Array1, typename Array2> auto concat(Array1 a1, Array2 a2)

    Concatenates the contents of two arrays ``a1`` and ``a2``. The pseudocode
    for this operation is

    .. code-block:: cpp

        Array<value_t<Array1>, Array1::Size + Array2::Size> out;
        for (size_t i = 0; i<Array1::Size; ++i)
            out[i] = a1[i];
        for (size_t i = 0; i<Array2::Size; ++i)
            out[i + Array1::Size] = a2[i];
        return out;

.. cpp:function:: template <typename Array> auto low(Array a)

    Returns the low part of the input array ``a``. The length of the low part
    is defined as the largest power of two that is smaller than
    ``Array::Size``. For power-of-two sized input, this function simply returns
    the low half.

.. cpp:function:: template <typename Array> auto high(Array a)

    Returns the high part of the input array ``a``. The length of the high part
    is equal to ``Array::Size`` minus the size of the low part. For
    power-of-two sized input, this function simply returns the high half.

.. cpp:function:: template <size_t Size, typename Array> auto head(Array a)

    Returns a new array containing the leading ``Size`` elements of ``a``.

.. cpp:function:: template <size_t Size, typename Array> auto tail(Array a)

    Returns a new array containing the trailing ``Size`` elements of ``a``.

.. cpp:function:: template <typename Outer, typename Inner> replace_scalar_t<Outer, Inner> full(const Inner &inner)

    Given an array type ``Outer`` and a value of type ``Inner``, this function
    returns a new composite type of shape ``[<shape of Outer>, <shape of
    Inner>]`` that is filled with the value of the ``inner`` argument.

    In the simplest case, this can be used to create a (potentially nested)
    array that is filled with constant values.

    .. code-block:: cpp

        using Vector4f = Array<float, 4>;
        using MyMatrix = Array<Vector4f, 4>;
        MyMatrix result = full<MyMatrix>(10.f);
        std::cout << result << std::endl;

        /* Prints:

             [[10, 10, 10, 10],
              [10, 10, 10, 10],
              [10, 10, 10, 10],
              [10, 10, 10, 10]]
        */

    Another use case entails replicating an array over the trailing dimensions
    of a new array:

    .. code-block:: cpp

        result = full<Vector4f>(Vector4f(1, 2, 3, 4))
        std::cout << result << std::endl;

        /* Prints:

            [[1, 1, 1, 1],
             [2, 2, 2, 2],
             [3, 3, 3, 3],
             [4, 4, 4, 4]]
        */

    Note how this is different from the default broadcasting behavior of
    arrays. In this case, ``Vector4f`` and ``MyMatrix`` have the same size
    in the leading dimension, which would replicate the vector over that
    axis instead:

    .. code-block:: cpp

        result = MyMatrix(Vector4f(1, 2, 3, 4));
        std::cout << result << std::endl;

        /* Prints:

            [[1, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 4]]
        */

.. cpp:function:: template <typename Array> std::array<size_t, array_depth_v<Array>> shape(const Array &a)

    Returns a ``std::array``, whose entries describe the shape of the
    (potentially multi-dimensional) input array along each dimension. It works
    for both static and dynamic arrays.


Operations for dynamic arrays
-----------------------------

.. cpp:function:: template <typename DArray> auto packet(DArray &&a, size_t i)

    Extracts the :math:`i`-th packet from a dynamic array or data structure. See
    the chapter on :ref:`dynamic arrays <dynamic>` on how to use this function.

.. cpp:function:: template <typename DArray> size_t packets(const DArray &a)

    Return the number of packets stored in the given dynamic array or data structure.

.. cpp:function:: template <typename DArray> auto slice(DArray &&a, size_t i)

    Extracts the :math:`i`-th slice from a dynamic array or data structure. See
    the chapter on :ref:`dynamic arrays <dynamic>` on how to use this function.

.. cpp:function:: template <typename DArray> size_t slices(const DArray &a)

    Return the number of packets stored in the given dynamic array or data structure.

.. cpp:function:: template <typename DArray> void set_slices(DArray &a, size_t size)

    Resize the given dynamic array or data structure so that there is space for
    ``size`` slices. When reducing the size of the array, any memory allocated
    so far is kept. Since there's no exponential allocation mechanism, it is not
    recommended to call ``set_slices`` repeatedly to append elements.

    Unlike ``std::vector::resize()``, previous values are *not* preserved when
    enlarging the array.

.. _type-traits:

Type traits
-----------

The following type traits are available to query the properties of arrays at
compile time.

Accessing types related to Enoki arrays
***************************************

.. cpp:type:: template <typename T> value_t

    Given an Enoki array ``T``, :cpp:type:`value_t\<T>` provides access to the
    type of the individual array entries. For non-array types ``T``,
    :cpp:type:`value_t\<T>` equals to the input template parameter ``T``.

    A few examples are shown below:

    .. code-block:: cpp

        // Non-array input:
        // value_t<float> yields 'float'

        // Array input:
        using FloatP     = Array<float>;
        using Vector4f   = Array<float, 4>;
        // value_t<Vector4f> yields 'float'

        using Vector4fr  = Array<float&, 4>;
        // value_t<Vector4fr> yields 'float&'

        using Vector4fP  = Array<FloatP, 4>;
        // value_t<Vector4P> yields 'FloatP'

        using Vector4fPr = Array<FloatP&, 4>;
        // value_t<Vector4Pr> yields 'FloatP&'


.. cpp:type:: template <typename... Args> expr_t

    Given arrays ``a1``, ..., ``an`` of type ``T1``, ..., ``Tn``, ``expr_t<T1,
    .., Tn>`` returns the type of an arithmetic expression such as ``a1 + ... +
    an``. The type trait applies all of the standard C++ type promotion rules
    and strips away references occurring anywhere within the definition of the
    input types. In addition to array input, ``expr_t`` also works for
    recursively defined arrays and non-array inputs or mixed array & non-array
    input.

    A few examples are shown below:

    .. code-block:: cpp

        using FloatP     = Array<float>;
        using DoubleP    = Array<double, FloatP::Size>;

        using Vector4f   = Array<float, 4>;
        using Vector4fr  = Array<float&, 4>;

        using Vector4fP  = Array<FloatP, 4>;
        using Vector4fPr = Array<FloatP&, 4>;

        using Vector4d   = Array<double, 4>;
        using Vector4dP  = Array<DoubleP, 4>;

        /* Non-array input */
        static_assert(std::is_same_v<expr_t<float>,               float>);
        static_assert(std::is_same_v<expr_t<float&>,              float>);
        static_assert(std::is_same_v<expr_t<float, double>,       double>);
        static_assert(std::is_same_v<expr_t<float&, double&>,     double>);

        /* Array input */
        static_assert(std::is_same_v<expr_t<Vector4f>,            Vector4f>);
        static_assert(std::is_same_v<expr_t<Vector4fr>,           Vector4f>);
        static_assert(std::is_same_v<expr_t<Vector4fP>,           Vector4fP>);
        static_assert(std::is_same_v<expr_t<Vector4fPr>,          Vector4fP>);

        static_assert(std::is_same_v<expr_t<Vector4f, double>,    Vector4d>);
        static_assert(std::is_same_v<expr_t<Vector4fPr, double&>, Vector4dP>);

.. cpp:type:: template <typename T> scalar_t

    Given a (potentially nested) Enoki array ``T``, this trait class provides
    access to the scalar type underlying the array.
    For non-array
    types ``T``, :cpp:type:`scalar_t\<T>` is simply set to the template parameter ``T``.

    A few examples are shown below:

    .. code-block:: cpp

        using FloatP     = Array<float>;

        using Vector4f   = Array<float, 4>;
        using Vector4fr  = Array<float&, 4>;

        using Vector4fP  = Array<FloatP, 4>;
        using Vector4fPr = Array<FloatP&, 4>;

        /* Non-array input */
        static_assert(std::is_same_v<scalar_t<float>,               float>);
        static_assert(std::is_same_v<scalar_t<float&>,              float>);

        /* Array input */
        static_assert(std::is_same_v<scalar_t<Vector4f>,            float>);
        static_assert(std::is_same_v<scalar_t<Vector4fr>,           float>);
        static_assert(std::is_same_v<scalar_t<Vector4fP>,           float>);
        static_assert(std::is_same_v<scalar_t<Vector4fPr>,          float>);

.. cpp:type:: template <typename T> mask_t

    Given an Enoki array ``T``, :cpp:type:`mask_t\<T>` provides access to the
    underlying mask type (i.e. the type that would result from a comparison
    operation such as ``array < 0``). For non-array types ``T``,
    :cpp:type:`mask_t\<T>` is set to ``bool``.


.. cpp:class:: template <typename T> array_depth

    .. cpp:member:: static constexpr size_t value

        Given a type :cpp:any:`T` (which could be a nested Enoki array),
        :cpp:member:`value` specifies the nesting level and stores it in the
        :cpp:var:`value` member. Non-array types (e.g. ``int32_t``) have a
        nesting level of 0, a type such as ``Array<float>`` has nesting level
        1, and so on.

Replacing the scalar type of an array
*************************************

The :cpp:type:`enoki::replace_scalar_t` type trait and various aliases construct arrays
matching a certain layout, but with different-flavored data. This is often
helpful when defining custom data structures or function inputs. See the
section on :ref:`custom data structures <custom-structures>` for an example
usage.

.. cpp:type:: template <typename Array, typename Scalar> replace_scalar_t

    Replaces the scalar type underlying an array. For instance,
    ``replace_scalar_t<Array<Array<float, 16>, 32>, int>`` is equal to ``Array<Array<int,
    16>, 32>``.

    The type trait also works for scalar arguments. Pointers and reference
    arguments are copied---for instance, ``replace_scalar_t<const float *, int>`` is
    equal to ``const int *``.

.. cpp:type:: template <typename Array> uint32_array_t = replace_scalar_t<Array, uint32_t>

    Create a 32-bit unsigned integer array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> int32_array_t = replace_scalar_t<Array, int32_t>

    Create a 32-bit signed integer array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> uint64_array_t = replace_scalar_t<Array, uint64_t>

    Create a 64-bit unsigned integer array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> int64_array_t = replace_scalar_t<Array, int64_t>

    Create a 64-bit signed integer array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> int_array_t

    Create a signed integer array (with the same number of bits per entry as
    the input) matching the layout of ``Array``.

.. cpp:type:: template <typename Array> uint_array_t

    Create an unsigned integer array (with the same number of bits per entry as
    the input) matching the layout of ``Array``.

.. cpp:type:: template <typename Array> float16_array_t = replace_scalar_t<Array, half>

    Create a half precision array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> float32_array_t = replace_scalar_t<Array, float>

    Create a single precision array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> float64_array_t = replace_scalar_t<Array, double>

    Create a double precision array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> float_array_t

    Create a floating point array (with the same number of bits per entry as
    the input) matching the layout of ``Array``.

.. cpp:type:: template <typename Array> bool_array_t = replace_scalar_t<Array, bool>

    Create a boolean array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> size_array_t = replace_scalar_t<Array, size_t>

    Create a ``size_t``-valued array matching the layout of ``Array``.

.. cpp:type:: template <typename Array> ssize_array_t = replace_scalar_t<Array, ssize_t>

    Create a ``ssize_t``-valued array matching the layout of ``Array``.

SFINAE helper types
-------------------

The following section discusses helper types that can be used to selectively
enable or disable template functions for Enoki arrays, e.g. like so:

.. code-block:: cpp

    template <typename Value, enable_if_array_t<Value> = 0>
    void f(Value value) {
        /* Invoked if 'Value' is an Enoki array */
    }

    template <typename Value, enable_if_not_array_t<Value> = 0>
    void f(Value value) {
        /* Invoked if 'Value' is *not* an Enoki array */
    }


Detecting Enoki arrays
**********************

.. cpp:class:: template <typename T> is_array

    .. cpp:member:: static constexpr bool value

        Equal to ``true`` iff ``T`` is a static or dynamic Enoki array type.

.. cpp:type:: template <typename T> enable_if_array_t = std::enable_if_t<is_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for Enoki
    array types.

.. cpp:type:: template <typename T> enable_if_not_array_t = std::enable_if_t<!is_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for types
    that are not Enoki arrays.


Detecting Enoki masks
*********************

.. cpp:class:: template <typename T> is_mask

    .. cpp:member:: static constexpr bool value

        Equal to ``true`` iff ``T`` is a static or dynamic Enoki mask type.

.. cpp:type:: template <typename T> enable_if_mask_t = std::enable_if_t<is_mask_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for Enoki
    mask types.

.. cpp:type:: template <typename T> enable_if_not_mask_t = std::enable_if_t<!is_mask_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for types
    that are not Enoki masks.

Detecting static Enoki arrays
*****************************

.. cpp:class:: template <typename T> is_static_array

    .. cpp:member:: static constexpr bool value

        Equal to ``true`` iff ``T`` is a static Enoki array type.

.. cpp:type:: template <typename T> enable_if_static_array_t = std::enable_if_t<is_static_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for
    static Enoki array types.

.. cpp:type:: template <typename T> enable_if_not_static_array_t = std::enable_if_t<!is_static_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for
    static Enoki array types.

Detecting dynamic Enoki arrays
******************************

.. cpp:class:: template <typename T> is_dynamic_array

    .. cpp:member:: static constexpr bool value

        Equal to ``true`` iff ``T`` is a dynamic Enoki array type.

.. cpp:type:: template <typename T> enable_if_dynamic_array_t = std::enable_if_t<is_dynamic_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for
    dynamic Enoki array types.

.. cpp:type:: template <typename T> enable_if_not_dynamic_array_t = std::enable_if_t<!is_dynamic_array_v<T>, int>

    SFINAE alias to selectively enable a class or function definition for
    dynamic Enoki array types.

.. cpp:class:: template <typename T> is_dynamic

    .. cpp:member:: static constexpr bool value

        Equal to ``true`` iff ``T`` (which could be a nested Enoki array) contains
        a dynamic array at *any* level.

        This is different from :cpp:class:`is_dynamic_array`, which only cares
        about the outermost level -- for instance, given static array ``T``
        containing a nested dynamic array, ``is_dynamic_array_v<T> ==
        false``, while ``is_dynamic_v<T> == true``.



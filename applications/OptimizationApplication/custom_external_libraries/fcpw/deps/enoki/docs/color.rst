.. cpp:namespace:: enoki

Color space transformations
===========================

Enoki provides a set of helper functions for color space transformations. For
now, only sRGB and inverse sRGB gamma correction are available. To use them,
include the following header file:

.. code-block:: cpp

    #include <enoki/color.h>


Functions
*********

.. cpp:function:: template <typename Value> Value linear_to_srgb(Value value)

    Efficiently applies the sRGB gamma correction

    .. math ::

        x\mapsto\begin{cases}12.92x,&x\leq 0.0031308\\1.055x^{1/2.4}-0.055,&x>0.0031308\end{cases}

    to an input value in the interval :math:`(0, 1)`.

.. cpp:function:: template <typename Value> Value srgb_to_linear(Value value)

    Efficiently applies the inverse sRGB gamma correction

    .. math ::

        x\mapsto{\begin{cases}{\frac {x}{12.92}},&x\leq 0.04045\\\left({\frac {x+0.055}{1.055}}\right)^{2.4},&x>0.04045\end{cases}}

    to an input value in the interval :math:`(0, 1)`.


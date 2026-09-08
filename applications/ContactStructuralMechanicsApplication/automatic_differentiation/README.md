Define dependencies into symbols, is not longer available since Sympy 1.3:

https://github.com/sympy/sympy/wiki/Release-Notes-for-1.3

> Symbols no longer automatically convert to functions when called, e.g., if f = Symbol('f'), f(t) is now a TypeError. To create a function, use f = Function('f') or f = symbols('f', cls=Function).

To solve that temporally you can install the 1.2 version of Sympy: (in order to run the AD scripts contained on this folder)

~~~sh
pip install sympy==1.2 
~~~

or 

~~~sh
python3 -m pip install sympy==1.2
~~~

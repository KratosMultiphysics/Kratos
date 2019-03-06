Dependencies
------------

VexCL depends on [Boost](http://www.boost.org). In particular, you may need to
link to Boost.chrono, Boost.date_time, Boost.filesystem, Boost.system,
Boost.thread.

Boost versions 1.48 and above should work. However, there is a known issue with
Boost v1.49 through 1.51 resulting in a segmentation fault in code that uses
sparse matrix vector products or stencil convolutions.

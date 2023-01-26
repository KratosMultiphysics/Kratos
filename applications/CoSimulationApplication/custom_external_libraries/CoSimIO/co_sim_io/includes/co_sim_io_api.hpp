
#ifndef CO_SIM_IO_API_H
#define CO_SIM_IO_API_H

#ifdef CO_SIM_IO_STATIC_DEFINE
#  define CO_SIM_IO_API
#  define CO_SIM_IO_NO_EXPORT
#else
#  ifndef CO_SIM_IO_API
#    ifdef co_sim_io_EXPORTS
        /* We are building this library */
#      define CO_SIM_IO_API __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define CO_SIM_IO_API __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef CO_SIM_IO_NO_EXPORT
#    define CO_SIM_IO_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef CO_SIM_IO_DEPRECATED
#  define CO_SIM_IO_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CO_SIM_IO_DEPRECATED_EXPORT
#  define CO_SIM_IO_DEPRECATED_EXPORT CO_SIM_IO_API CO_SIM_IO_DEPRECATED
#endif

#ifndef CO_SIM_IO_DEPRECATED_NO_EXPORT
#  define CO_SIM_IO_DEPRECATED_NO_EXPORT CO_SIM_IO_NO_EXPORT CO_SIM_IO_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CO_SIM_IO_NO_DEPRECATED
#    define CO_SIM_IO_NO_DEPRECATED
#  endif
#endif

#endif /* CO_SIM_IO_API_H */

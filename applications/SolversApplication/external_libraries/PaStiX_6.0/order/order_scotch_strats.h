/**
 *
 * @file order_scotch_strats.h
 *
 * PaStiX order (PT-)Scotch strategy strings
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _order_scotch_strats_h_
#define _order_scotch_strats_h_

#define SCOTCH_STRAT_DIRECT                             \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>120)?m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"      \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>120)?(m{rat=0.8,"               \
    ""                         "vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}})|"      \
    ""                        "m{rat=0.8,"              \
    ""                          "vert=100,"             \
    ""                          "low=h{pass=10},"       \
    ""                          "asc=f{bal=0.2}};,"     \
    ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"     \
    ""      "ose=g}}"

#define SCOTCH_STRAT_INCOMP                             \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>120)?m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"      \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>120)?(m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}})|"      \
    ""                       "m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}};,"      \
    ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"     \
    ""            "ose=g}}"

#define SCOTCH_STRAT_PERSO                              \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>%ld)?m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"         \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>%ld)?(m{rat=0.8,"               \
    ""                         "vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}})|"      \
    ""                        "m{rat=0.8,"              \
    ""                          "vert=100,"             \
    ""                          "low=h{pass=10},"       \
    ""                          "asc=f{bal=0.2}};,"     \
    ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"         \
    ""      "ose=g}}"

#define PTSCOTCH_STRAT_DIRECT                           \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>120)?m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{rat=0.8,"                \
    ""                        "vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"       \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>120)?(m{type=h,"                \
    ""                         "rat=0.8,"               \
    ""                         "vert=100000,"           \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=08.2}})|"     \
    ""                       "m{type=h,"                \
    ""                         "rat=0.8,"               \
    ""                         "vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}};,"      \
    ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"     \
    ""      "ose=g}}"

#define PTSCOTCH_STRAT_INCOMP                           \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>120)?m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"      \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>120)?(m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}})|"      \
    ""                       "m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}};,"      \
    ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"     \
    ""      "ose=g}}"

#define PTSCOTCH_STRAT_PERSO                            \
    "c{rat=0.7,"                                        \
    """cpr=n{sep=/(vert>%ld)?m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}}|"        \
    ""                      "m{vert=100,"               \
    ""                        "low=h{pass=10},"         \
    ""                        "asc=f{bal=0.2}};,"       \
    ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"         \
    ""      "ose=g},"                                   \
    """unc=n{sep=/(vert>%ld)?(m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}})|"      \
    ""                       "m{vert=100,"              \
    ""                         "low=h{pass=10},"        \
    ""                         "asc=f{bal=0.2}};,"      \
    ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"         \
    ""      "ose=g}}"


#endif /* _order_scotch_strats_h_ */

/**
 *
 * @file simu_timer.h
 *
 * PaStiX simulation timer.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_simu
 * @{
 *
 **/
#ifndef _simu_timer_h_
#define _simu_timer_h_

/**
 * @brief Timer for the simulation.
 */
typedef struct simu_timer_s {
    double s; /**< Second in the timer */
    /*  double ms;*/
} SimuTimer;

/**
 * @brief Compare two timings
 * @param[in] t1
 *            The first timer.
 * @param[in] t2
 *            The second timer.
 * @return True if the t1 is smaller than t2, false otherwise.
 */
static inline int
timerComp(const SimuTimer *t1,
          const SimuTimer *t2)
{
    /* Return (t1 < t2) */
    if(t1->s < t2->s) {
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * @brief Increment the timer
 * @param[inout] timer
 *               The timer to update.
 * @param[in]    t
 *               The time to add to the timer.
 */
static inline void
timerAdd(SimuTimer *timer, double t)
{
    timer->s += t;
}

/**
 * @brief Get the timer value
 * @param[in] timer
 *            The timer to read.
 * @return The timer value in second.
 */
static inline double
timerVal(const SimuTimer *timer)
{
    return timer->s;
}

/**
 * @brief Set the timer value
 * @param[inout] timer
 *               The timer to set
 * @param[in]    t
 *               The value to set
 */
static inline void
timerSet(SimuTimer *timer, double t)
{
    timer->s = t;
}

/**
 * @brief Set the timer value if the value is greater than the actual one.
 * @param[inout] timer
 *               The timer to update
 * @param[in]    t
 *               The time to compare with and set if larger than the timer.
 */
static inline void
timerSetMax(SimuTimer *timer, double t)
{
    if ( t > timer->s ) {
        timer->s = t;
    }
}

#endif /* _simu_timer_h_ */

/**
 *@}
 */

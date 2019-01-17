#include "f2c.h"

/* Subroutine */ int dlaruv_slu(integer *iseed, integer *n, doublereal *x)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLARUV returns a vector of n random real numbers from a uniform (0,1) 
  
    distribution (n <= 128).   

    This is an auxiliary routine called by DLARNV and ZLARNV.   

    Arguments   
    =========   

    ISEED   (input/output) INTEGER array, dimension (4)   
            On entry, the seed of the random number generator; the array 
  
            elements must be between 0 and 4095, and ISEED(4) must be   
            odd.   
            On exit, the seed is updated.   

    N       (input) INTEGER   
            The number of random numbers to be generated. N <= 128.   

    X       (output) DOUBLE PRECISION array, dimension (N)   
            The generated random numbers.   

    Further Details   
    ===============   

    This routine uses a multiplicative congruential method with modulus   
    2**48 and multiplier 33952834046453 (see G.S.Fishman,   
    'Multiplicative congruential random number generators with modulus   
    2**b: an exhaustive analysis for b = 32 and a partial analysis for   
    b = 48', Math. Comp. 189, pp 331-344, 1990).   

    48-bit integers are stored in 4 integer array elements with 12 bits   
    per element. Hence the routine is portable across machines with   
    integers of 32 bits or more.   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Initialized data */
    static integer mm[512]	/* was [128][4] */ = { 494,2637,255,2008,1253,
	    3344,4084,1739,3143,3468,688,1657,1238,3166,1292,3422,1270,2016,
	    154,2862,697,1706,491,931,1444,444,3577,3944,2184,1661,3482,657,
	    3023,3618,1267,1828,164,3798,3087,2400,2870,3876,1905,1593,1797,
	    1234,3460,328,2861,1950,617,2070,3331,769,1558,2412,2800,189,287,
	    2045,1227,2838,209,2770,3654,3993,192,2253,3491,2889,2857,2094,
	    1818,688,1407,634,3231,815,3524,1914,516,164,303,2144,3480,119,
	    3357,837,2826,2332,2089,3780,1700,3712,150,2000,3375,1621,3090,
	    3765,1149,3146,33,3082,2741,359,3316,1749,185,2784,2202,2199,1364,
	    1244,2020,3160,2785,2772,1217,1822,1245,2252,3904,2774,997,2573,
	    1148,545,322,789,1440,752,2859,123,1848,643,2405,2638,2344,46,
	    3814,913,3649,339,3808,822,2832,3078,3633,2970,637,2249,2081,4019,
	    1478,242,481,2075,4058,622,3376,812,234,641,4005,1122,3135,2640,
	    2302,40,1832,2247,2034,2637,1287,1691,496,1597,2394,2584,1843,336,
	    1472,2407,433,2096,1761,2810,566,442,41,1238,1086,603,840,3168,
	    1499,1084,3438,2408,1589,2391,288,26,512,1456,171,1677,2657,2270,
	    2587,2961,1970,1817,676,1410,3723,2803,3185,184,663,499,3784,1631,
	    1925,3912,1398,1349,1441,2224,2411,1907,3192,2786,382,37,759,2948,
	    1862,3802,2423,2051,2295,1332,1832,2405,3638,3661,327,3660,716,
	    1842,3987,1368,1848,2366,2508,3754,1766,3572,2893,307,1297,3966,
	    758,2598,3406,2922,1038,2934,2091,2451,1580,1958,2055,1507,1078,
	    3273,17,854,2916,3971,2889,3831,2621,1541,893,736,3992,787,2125,
	    2364,2460,257,1574,3912,1216,3248,3401,2124,2762,149,2245,166,466,
	    4018,1399,190,2879,153,2320,18,712,2159,2318,2091,3443,1510,449,
	    1956,2201,3137,3399,1321,2271,3667,2703,629,2365,2431,1113,3922,
	    2554,184,2099,3228,4012,1921,3452,3901,572,3309,3171,817,3039,
	    1696,1256,3715,2077,3019,1497,1101,717,51,981,1978,1813,3881,76,
	    3846,3694,1682,124,1660,3997,479,1141,886,3514,1301,3604,1888,
	    1836,1990,2058,692,1194,20,3285,2046,2107,3508,3525,3801,2549,
	    1145,2253,305,3301,1065,3133,2913,3285,1241,1197,3729,2501,1673,
	    541,2753,949,2361,1165,4081,2725,3305,3069,3617,3733,409,2157,
	    1361,3973,1865,2525,1409,3445,3577,77,3761,2149,1449,3005,225,85,
	    3673,3117,3089,1349,2057,413,65,1845,697,3085,3441,1573,3689,2941,
	    929,533,2841,4077,721,2821,2249,2397,2817,245,1913,1997,3121,997,
	    1833,2877,1633,981,2009,941,2449,197,2441,285,1473,2741,3129,909,
	    2801,421,4073,2813,2337,1429,1177,1901,81,1669,2633,2269,129,1141,
	    249,3917,2481,3941,2217,2749,3041,1877,345,2861,1809,3141,2825,
	    157,2881,3637,1465,2829,2161,3365,361,2685,3745,2325,3609,3821,
	    3537,517,3017,2141,1537 };
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i, i1, i2, i3, i4, it1, it2, it3, it4;


#define MM(I) mm[(I)]
#define WAS(I) was[(I)]
#define ISEED(I) iseed[(I)-1]
#define X(I) x[(I)-1]



    i1 = ISEED(1);
    i2 = ISEED(2);
    i3 = ISEED(3);
    i4 = ISEED(4);

    i__1 = min(*n,128);
    for (i = 1; i <= min(*n,128); ++i) {

/*        Multiply the seed by i-th power of the multiplier modulo 2**
48 */

	it4 = i4 * MM(i + 383);
	it3 = it4 / 4096;
	it4 -= it3 << 12;
	it3 = it3 + i3 * MM(i + 383) + i4 * MM(i + 255);
	it2 = it3 / 4096;
	it3 -= it2 << 12;
	it2 = it2 + i2 * MM(i + 383) + i3 * MM(i + 255) + i4 * MM(i + 127);
	it1 = it2 / 4096;
	it2 -= it1 << 12;
	it1 = it1 + i1 * MM(i + 383) + i2 * MM(i + 255) + i3 * MM(i + 127) + 
		i4 * MM(i - 1);
	it1 %= 4096;

/*        Convert 48-bit integer to a real number in the interval (0,1
) */

	X(i) = ((doublereal) it1 + ((doublereal) it2 + ((doublereal) it3 + (
		doublereal) it4 * 2.44140625e-4) * 2.44140625e-4) * 
		2.44140625e-4) * 2.44140625e-4;
/* L10: */
    }

/*     Return final value of seed */

    ISEED(1) = it1;
    ISEED(2) = it2;
    ISEED(3) = it3;
    ISEED(4) = it4;
    return 0;

/*     End of DLARUV */

} /* dlaruv_slu */


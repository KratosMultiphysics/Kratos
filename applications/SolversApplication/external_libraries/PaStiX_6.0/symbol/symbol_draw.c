/**
 *
 * @file symbol_draw.c
 *
 * PaStiX symbol structure drawing function.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 **/
#include <time.h>
#include "common.h"
#include "symbol.h"

/**
 * @addtogroup symbol_dev_draw
 * @{
 *
 */

/**
 * @name PostScript parameter
 * @{
 *   PostScript (tm) output definitions.
 */

/**
 * @brief PostScript dots-per-inch
 */
#define SYMBOL_PSDPI      72
/**
 * @brief PostScript picture size (in inches)
 */
#define SYMBOL_PSPICTSIZE 6.6

/**
 * @}
 */

/**
 * @brief Predefined set of colors
 */
static float symbolDrawColorTab[16][3] = {
    { 1.00, 0.00, 0.00 }, /* Red          */
    { 0.00, 1.00, 0.00 }, /* Green        */
    { 1.00, 1.00, 0.00 }, /* Yellow       */
    { 0.00, 0.00, 1.00 }, /* Blue         */
    { 1.00, 0.00, 1.00 }, /* Magenta      */
    { 0.00, 1.00, 1.00 }, /* Cyan         */
    { 1.00, 0.50, 0.20 }, /* Orange       */
    { 0.30, 0.55, 0.00 }, /* Olive        */
    { 0.72, 0.47, 0.47 }, /* Dark pink    */
    { 0.33, 0.33, 0.81 }, /* Sea blue     */
    { 1.00, 0.63, 0.63 }, /* Pink         */
    { 0.62, 0.44, 0.65 }, /* Violet       */
    { 0.60, 0.80, 0.70 }, /* Pale green   */
    { 0.47, 0.20, 0.00 }, /* Brown        */
    { 0.00, 0.68, 0.68 }, /* Turquoise    */
    { 0.81, 0.00, 0.40 } }; /* Purple     */

/**
 *******************************************************************************
 *
 * @brief Return one of 16 predefined, visually distinct distinct colors.
 *
 *******************************************************************************
 *
 * @param[in]    labl
 *          The label of the color
 *
 * @param[inout] color
 *          The color array in which is returned the selected color.
 *
 *******************************************************************************/
void
pastixSymbolDrawColor ( const pastix_int_t labl,
                        float              color[] )
{
    color[0] = (float) symbolDrawColorTab[(labl - 1) % 16][0];
    color[1] = (float) symbolDrawColorTab[(labl - 1) % 16][1];
    color[2] = (float) symbolDrawColorTab[(labl - 1) % 16][2];
}

/**
 *******************************************************************************
 *
 * @brief Export the symbol structure in a PostScript format.
 *
 * This routine writes to the given stream a PostScript (tm) picture of the
 * symbolic block matrix.
 *
 *******************************************************************************
 *
 * @param[in]    symbptr
 *          The pointer to the symbolic structure to draw.
 *
 * @param[in]    diagfunc
 *          Personal function to draw the diagonal block. NULL by default
 *
 * @param[in]    offdfunc
 *          Personal function to draw the off-diagonal blocks. NULL by default
 *
 * @param[in]    dataptr
 *          Data structure for block coloring
 *
 * @param[inout] stream
 *          The stream where to write the output.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval -1 on error
 *
 *******************************************************************************/
int
pastixSymbolDrawFunc (
    const symbol_matrix_t * const  symbptr,
    int                         (* diagfunc) (const symbol_matrix_t * const, const symbol_blok_t * const, void * const, float * const),
    int                         (* offdfunc) (const symbol_matrix_t * const, const symbol_blok_t * const, void * const, float * const),
    void * const                   dataptr,              /* Data structure for block coloring */
    FILE * const                   stream)
{
    pastix_int_t cblknum;                    /* Number of current column block */
    pastix_int_t bloknum;                    /* Number of current block        */
    time_t       picttime;                   /* Creation time                  */
    double       pictsize;                   /* Number of distinct coordinates */
    int          o;

    time (&picttime);                               /* Get current time */
    pictsize = (double) (symbptr->nodenbr + 1);     /* Get matrix size  */

    fprintf (stream, "%%!PS-Adobe-2.0 EPSF-2.0\n"); /* Write header */
    fprintf (stream, "%%%%Title: pastixSymbolmatrix (%ld,%ld,%ld)\n",
             (long) symbptr->cblknbr, (long) symbptr->bloknbr, (long)symbptr->nodenbr);
    fprintf (stream, "%%%%Creator: pastixSymbolDraw (LaBRI, Universite Bordeaux I)\n");
    fprintf (stream, "%%%%CreationDate: %s", ctime (&picttime));
    fprintf (stream, "%%%%BoundingBox: 0 0 %ld %ld\n",
             (long) (SYMBOL_PSPICTSIZE * SYMBOL_PSDPI),
             (long) (SYMBOL_PSPICTSIZE * SYMBOL_PSDPI));
    fprintf (stream, "%%%%Pages: 0\n");
    fprintf (stream, "%%%%EndComments\n");          /* Write shortcuts */
    fprintf (stream, "/c { 4 2 roll pop pop newpath 2 copy 2 copy moveto dup lineto dup lineto closepath fill } bind def\n");
    fprintf (stream, "/b { 4 copy 2 index exch moveto lineto dup 3 index lineto exch lineto closepath fill pop } bind def\n");
    fprintf (stream, "/r { setrgbcolor } bind def\n");
    fprintf (stream, "/g { setgray } bind def\n");

    fprintf (stream, "gsave\n");                    /* Save context         */
    fprintf (stream, "0 setlinecap\n");             /* Use miter caps       */
    fprintf (stream, "%f dup scale\n",              /* Print scaling factor */
             (double) SYMBOL_PSDPI * SYMBOL_PSPICTSIZE / pictsize);
    fprintf (stream, "[ 1 0 0 -1 0 %d ] concat\n",  /* Reverse Y coordinate */
             (int) (symbptr->nodenbr + 1));

    fprintf (stream, "0 0\n");                      /* Output fake column block */
    for (cblknum = 0, bloknum = 0; cblknum < symbptr->cblknbr; cblknum ++) {
        float               coloval[3];               /* Color of diagonal block and previous color */
        pastix_int_t                 blokend;                  /* Number of end block for column             */

        coloval[0] =
            coloval[1] =
            coloval[2] = 0.5;
        if (diagfunc != NULL)                         /* Always display diagonal blocks */
            diagfunc (symbptr, &symbptr->bloktab[bloknum], dataptr, coloval);
        if ((coloval[0] == coloval[1]) &&
            (coloval[1] == coloval[2]))
            fprintf (stream, "%.2g g ",
                     (float) coloval[0]);
        else
            fprintf (stream, "%.2g %.2g %.2g r \n",
                     (float) coloval[0], (float) coloval[1], (float) coloval[2]);

        fprintf (stream, "%ld\t%ld\tc\n",             /* Begin new column block */
                 (long) (symbptr->cblktab[cblknum].fcolnum - symbptr->baseval),
                 (long) (symbptr->cblktab[cblknum].lcolnum - symbptr->baseval + 1));

        for (bloknum ++, blokend = symbptr->cblktab[cblknum + 1].bloknum; /* Skip diagonal block */
             bloknum < blokend; bloknum ++) {
            float               colbval[3];             /* Color of off-diagonal block */

            colbval[0] =
                colbval[1] =
                colbval[2] = 0.0;
            if ((offdfunc == NULL) || (offdfunc (symbptr, &symbptr->bloktab[bloknum], dataptr, colbval) != 0)) { /* If block is kept */
                if ((coloval[0] != colbval[0]) ||         /* If change with respect to previous  color */
                    (coloval[1] != colbval[1]) ||
                    (coloval[2] != colbval[2])) {
                    coloval[0] = colbval[0];                /* Save new color data */
                    coloval[1] = colbval[1];
                    coloval[2] = colbval[2];

                    if ((coloval[0] == coloval[1]) &&
                        (coloval[1] == coloval[2]))
                        fprintf (stream, "%.2g g ",
                                 (float) coloval[0]);
                    else
                        fprintf (stream, "%.2g %.2g %.2g r \n",
                                 (float) coloval[0], (float) coloval[1], (float) coloval[2]);
                }

                fprintf (stream, "%ld\t%ld\tb\n",         /* Write block in column block */
                         (long) (symbptr->bloktab[bloknum].frownum - symbptr->baseval),
                         (long) (symbptr->bloktab[bloknum].lrownum - symbptr->baseval + 1));
            }
        }
    }
    fprintf (stream, "pop pop\n");                  /* Purge last column block indexes */
    o = fprintf (stream, "grestore\nshowpage\n");   /* Restore context                 */

    return ((o != EOF) ? 0 : 1);
}

/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Export the symbol structure in a PostScript format.
 *
 * This routine writes to the given stream a PostScript (tm) picture of the
 * symbolic block matrix, with diagonal blocks in black and off-diagonal blocks
 * in dark gray.
 *
 *******************************************************************************
 *
 * @param[in]    symbptr
 *          The pointer to the symbolic structure to draw.
 *
 * @param[inout] stream
 *          The stream where to write the output.
 *
 *******************************************************************************
 *
 * @retval 0 on success
 * @retval -1 on error
 *
 *******************************************************************************/
int
pastixSymbolDraw ( const symbol_matrix_t * const  symbptr,
                   FILE * const                stream )
{
    return (pastixSymbolDrawFunc (symbptr, NULL, NULL, NULL, stream));
}

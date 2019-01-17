/**
 *
 * @file pastix_starpu_model.h
 *
 * Pastix StarPU model function
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2018-07-16
 *
 **/

#ifndef _pastix_starpu_model_h_
#define _pastix_starpu_model_h_

double blok_getrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_hetrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_potrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_pxtrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_sytrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_trsmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double blok_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );

double cblk_getrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double cblk_hetrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double cblk_potrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double cblk_pxtrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double cblk_sytrf_cost ( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );
double cblk_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch, unsigned nimpl );

#endif /* _pastix_starpu_model_h_ */

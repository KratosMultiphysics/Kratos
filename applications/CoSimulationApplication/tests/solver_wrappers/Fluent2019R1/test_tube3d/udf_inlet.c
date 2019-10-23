#include "udf.h"

DEFINE_PROFILE(pressure_inlet, face_thread, nv) {
#if !RP_HOST
	face_t face;
	real time, pressure;
	
	time = CURRENT_TIME;
	pressure = ((time <= 0.003) ? 1333.2 : 0.0);
	
	begin_f_loop(face,face_thread) {
		F_PROFILE(face,face_thread,nv) = pressure;
	} end_f_loop(face,face_thread)
#endif /* !RP_HOST */
}

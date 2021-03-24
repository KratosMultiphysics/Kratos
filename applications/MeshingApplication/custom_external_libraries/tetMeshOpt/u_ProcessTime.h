#if !defined(U_PROCESSTIME_H )
#define  U_PROCESSTIME_H

#include <vector>
#include "u_delphiClasses.h"
#include <time.h>
#include <stdio.h>
/*
struct timeval 
{
time_t tv_sec;
suseconds_t tv_usec;
}
*/
class TTimeSignal 
{
public : 
	clock_t _start, _end ;
	double minTime, maxTime,acum ,last ,prop;
	int counter;
	TTimeSignal* parent;
	char* name;
	TTimeSignal(char* name){ this->name = name; acum = 0;}
	void startTime();
	void endTime();
	double elapsedTime(); // in ms

};

void startProcess(char* procName);
void stopTimers();
void startTimers();
double endProcessR(char* procName);
void endProcess(char* procName);
TTimeSignal* findProcess(char* procName);
void showProcessTime();

#endif
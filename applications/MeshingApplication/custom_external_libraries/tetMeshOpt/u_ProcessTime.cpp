#include "u_ProcessTime.h"
#include <string>
#include <iostream>

#include <time.h>
#include <stdio.h>

TList<TTimeSignal*>* processNames ; 
bool TimingIsActive = true ;

double TTimeSignal::elapsedTime()
{
	double diffticks=_end-_start;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
}

void TTimeSignal::startTime()
{
	_start=clock();
}
void TTimeSignal::endTime()
{
	_end=clock();
	double diffticks=_end-_start;
	this->acum += (diffticks*1000)/CLOCKS_PER_SEC;


}

void stopTimers()
{
	TimingIsActive = false;
}
void startTimers()
{
	TimingIsActive = true;
}


void startProcess(char* procName)
{
	if (!TimingIsActive ) return ;
	if (processNames == NULL ) processNames = new TList<TTimeSignal*>();
	TTimeSignal *newT = findProcess(procName);
	if (newT == NULL) 
	{
		newT = new TTimeSignal(procName);
		processNames->Add(newT);
	}
	newT->startTime();
}
void endProcess(char* procName)
{
	if (!TimingIsActive ) return ;
	if (processNames == NULL ) return;
	TTimeSignal *newT = findProcess(procName);
	if (newT == NULL)  return;
	newT->endTime();
}

void showProcessTime()
{
	if (!TimingIsActive ) return ;
	for (int i=0; i<processNames->Count();i++)
	{
		TTimeSignal *t = (TTimeSignal*)(processNames->elementAt(i));
		//printf("%s%f", t->name, t->elapsedTime());
		std::cout << "........................................"<<"\n";
		std::cout <<t->name<<" lastTime " << t->elapsedTime()<< " totalTime "<< t->acum<<"\n";
	}

}



TTimeSignal* findProcess(char* procName)
{
	std::string str1 (procName);
	for (int i = 0 ; i<processNames->Count(); i++)
	{	
		TTimeSignal *t = (TTimeSignal*)(processNames->elementAt(i));
		std::string str2 (t->name);
		if ( str1.compare( str2) == 0)
			return t;
	}
	return NULL;
}

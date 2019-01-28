#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <cstring>
#include "DPC.h"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 7)
    {
        cout<<"Error! Usage: ./DPC [filename] [number of data points] [number of clusters] [number of dc] [output clustering results: True or False] [output decision graphs: True or False]"<<endl;
        return 0;
    }
    //parse input parameters
    string file(argv[1]);
    int N = atoi(argv[2]);
    int K = atoi(argv[3]);
    int P = atoi(argv[4]);

    DPC dpc(N);
    dpc.readDist(file);
    dpc.setK(K); //set the number of clusters

    //set percentages of average number of neighbors compared with all data points
    //the range is from 1% to 5%, so the formula we use is: (1 + i*4/n)%
    vector<double> percents; 
    for(double i = 0; i < P; i++)
        percents.push_back(1 + i / P * 5);
    //vector<double> percents = {0.5};


    timespec t_start, t_end;                       //for the measurement of total running time
    timespec t_startDc, t_endDc;                   //for the measurement of time spent on getting Dc
    timespec t_startRhoDelta, t_endRhoDelta;       //for the measurement of time spent on updating Rho and Delta
    timespec t_startHalo, t_endHalo;               //for the measurement of time spent on identifying cluster halo
    
    long timeTotal = 0, timeDc = 0, timeRhoDelta = 0, timeHalo = 0;
	
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    //measure running time for getting Dc and change the unit to millisecond
    clock_gettime(CLOCK_MONOTONIC, &t_startDc);
    dpc.findDcPos(percents);
    clock_gettime(CLOCK_MONOTONIC, &t_endDc);
    timeDc = 1000000*(t_endDc.tv_sec-t_startDc.tv_sec)+(t_endDc.tv_nsec-t_startDc.tv_nsec)/1000; 
    
    //iterating through all percentage choices
    for(int i = 0; i < percents.size(); i++)
    {
        dpc.setPercent(percents[i]);
        //start point of measuring running time for updating rho and delta
        clock_gettime(CLOCK_MONOTONIC, &t_startRhoDelta);
        dpc.calRho(false);

        //for the first dc, DPC Express use the original method to get delta
        //starting from the second dc, DPC Express applies optimization strategies 
        if(i == 0)
        {
            dpc.sortRho();   
            dpc.calDelta();  
        }
        else
            dpc.updateDelta();

        //end point of measuring running time for updating rho and delta

        clock_gettime(CLOCK_MONOTONIC, &t_endRhoDelta);
        timeRhoDelta += 1000000*(t_endRhoDelta.tv_sec-t_startRhoDelta.tv_sec)+(t_endRhoDelta.tv_nsec-t_startRhoDelta.tv_nsec)/1000;

        if(strcmp(argv[6], "True") == 0)
            dpc.outputDecisionGraph();
        dpc.calGamma();
        
        dpc.pickCenter();

        dpc.assignation();
        clock_gettime(CLOCK_MONOTONIC, &t_startHalo);
        
        dpc.IdentifyHalo();
        if(strcmp(argv[5], "True") == 0)
            dpc.printClusterAssignment();

        clock_gettime(CLOCK_MONOTONIC, &t_endHalo);
        timeHalo += 1000000*(t_endHalo.tv_sec-t_startHalo.tv_sec)+(t_endHalo.tv_nsec-t_startHalo.tv_nsec)/1000;

    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    
	timeTotal = 1000000*(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_nsec-t_start.tv_nsec)/1000;
    
    printf("it took %lf ms\n", timeTotal / 1000.0);
    printf("Dc took %lf ms\n", timeDc / 1000.0);
    printf("Rho took %lf ms\n", timeRhoDelta / 1000.0);
    printf("Halo took %lf ms\n", timeHalo / 1000.0);

    return 0;
}

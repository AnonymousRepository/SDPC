#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <cstring>
#include "DPC.h"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 6)
    {
        cout<<"Error! Usage: ./DPC [filename] [number of data points] [number of clusters] [number of dc] [print clustering results: True or False]"<<endl;
        return 0;
    }
    
    string file(argv[1]);
    int N = atoi(argv[2]);
    int K = atoi(argv[3]);
    int P = atoi(argv[4]);

    timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    timespec t_startRho, t_endRho, t_startDelta, t_endDelta, t_startHalo, t_endHalo;
    timespec t_startDc, t_endDc;                   //for the measurement of time spent on getting Dc
	long timedifRho = 0, timedifDelta = 0, timedifHalo = 0, timedifDc = 0;

    DPC dpc(N);
    clock_gettime(CLOCK_MONOTONIC, &t_startDc);
    dpc.readDist(file);
    clock_gettime(CLOCK_MONOTONIC, &t_endDc);
    timedifDc = 1000000*(t_endDc.tv_sec-t_startDc.tv_sec)+(t_endDc.tv_nsec-t_startDc.tv_nsec)/1000;
    dpc.setK(K);

    if(strcmp(argv[5], "True") == 0)
        dpc.setPrintOption(true);
    else if(strcmp(argv[5], "False") == 0)
        dpc.setPrintOption(false);
    else
    {
        cout<<"Error! Print option should be \"True\" of \"False\" (no quotation marks)"<<endl;
        return 0;
    }

    vector<double> percents; 
    for(double i = 0; i < P; i++)
        percents.push_back(1 + i / P * 4);
    //vector<double> percents = {0.5};

    for(int i = 0; i < percents.size(); i++)
    {
        dpc.setPercent(percents[i]);

        clock_gettime(CLOCK_MONOTONIC, &t_startRho);
        dpc.calRho(false);
        clock_gettime(CLOCK_MONOTONIC, &t_endRho);
        dpc.sortRho();
        //clock_gettime(CLOCK_MONOTONIC, &t_endRho);
        timedifRho += 1000000*(t_endRho.tv_sec-t_startRho.tv_sec)+(t_endRho.tv_nsec-t_startRho.tv_nsec)/1000;

        clock_gettime(CLOCK_MONOTONIC, &t_startDelta);
        dpc.calDelta();
        clock_gettime(CLOCK_MONOTONIC, &t_endDelta);
        timedifDelta += 1000000*(t_endDelta.tv_sec-t_startDelta.tv_sec)+(t_endDelta.tv_nsec-t_startDelta.tv_nsec)/1000;
    
        dpc.calGamma();
        dpc.pickCenter();
        dpc.assignation();

        clock_gettime(CLOCK_MONOTONIC, &t_startHalo);
        dpc.findHalo();
        clock_gettime(CLOCK_MONOTONIC, &t_endHalo);
        timedifHalo += 1000000*(t_endHalo.tv_sec-t_startHalo.tv_sec)+(t_endHalo.tv_nsec-t_startHalo.tv_nsec)/1000;
    }
    clock_gettime(CLOCK_MONOTONIC, &t_end);
	long timedif = 1000000*(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_nsec-t_start.tv_nsec)/1000;
    printf("it took %lf ms\n", timedif / 1000.0);
    printf("Dc took %lf ms\n", timedifDc / 1000.0);
    printf("Rho took %lf ms\n", timedifRho / 1000.0);
    printf("Delta took %lf ms\n", timedifDelta / 1000.0);
    printf("Halo took %lf ms\n", timedifHalo / 1000.0);
    return 0;
}

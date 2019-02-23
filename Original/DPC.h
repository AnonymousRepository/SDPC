#ifndef _DPC_H
#define _DPC_H
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <limits>

using namespace std;

class DPC
{
    int N; //the number of points
    int K; //the number of clusters;
    vector<vector<double>> dist; //distance btw two points
    vector<double> dist2; //used for finding dc;
    double percent;
    double dc;
    vector<double> rho; //density of each point
    vector<double> delta; //the distance to the nearest neighbour with higher density
    vector<int> nneigh; //the nearest neighbour of each point
    //vector<double> rho_sorted; //sorted rho
    vector<int> ordrho; //index corresponding to sorted rho
    vector<double> gamma; //gammas
    vector<int> ordgamma;
    vector<int> cl; //cluster for each point 
    vector<int> icl; //index of cluster
    vector<int> halo; //-1 means halo, otherwise indicating which cluster it is in
    vector<double> bord_rho; //border density for each cluster
    bool printable; //indicating whether print the clustering results
    
public:
    DPC(int n);
    void readDist(string fileName);
    void setPercent(double Pct);
    double getDC();
    void calRho(bool isGaussKernel); //calculate the density for each point
    void gaussKernel();
    void cutOffKernel();
    void sortRho();
    void calDelta();
    void calGamma();
    void pickCenter();
    void assignation();
    void findHalo();
    //set the number of clusters
    void setK(int K);
    //set whether print clustering results
    void setPrintOption(bool para);
};
#endif

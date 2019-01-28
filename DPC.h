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
#include <queue>

using namespace std;

class DPC
{
    int N; //the number of points
    int K; //the number of clusters;
    int nDist; //the number of distances in the input file
    int nDistLessDc; //the number of distances in the input file less than dc
    int prevNDistLessDc; //nDistLessDc for last percent value
    vector<vector<double>> dist; //distance btw two points
    vector<pair<double, pair<int, int>>> distList; //list of all distances
    vector<pair<double, pair<int, int>>> omitted_dist; //for distances hava same value with the dc which are omitted from current iteration \
    but will be used in the next iteration.
    double percent; //the percentage of the average number of neighbors compared with all data points
    double dc; //cut-off distance
    vector<double> rho; //local density of each point
    vector<double> delta; //the distance to the nearest data point of higher density 
    vector<int> ndp; //the nearest data point of higher density
    vector<int> q; //index corresponding to sorted rho
    vector<double> gamma; //gammas
    vector<int> ordgamma; //index corresponding to sorted gamma
    vector<int> cl; //cluster for each point 
    vector<int> icl; //index of cluster center
    vector<int> halo; //-1 means halo, otherwise indicating which cluster it is in
    vector<double> bord_rho; //border density threshold for each cluster
    vector<vector<int>> nb; //neighbor list for every point with points within one dc
    vector<vector<int>> nbS; //list of neighbor size of each data point under different dc
    vector<double> dcs; //contain all the dc that has been calculated
    int nTrial; //indicating which trail is running

public:
    DPC(int n);
    //read distances from input file to dist and distList
    void readDist(string fileName); 
    //set the percentage to be tested 
    void setPercent(double Pct);
    //get current dc
    double getDC();
    //calculate the density for each data point
    void calRho(bool isGaussKernel); 
    //calculate local density with gaussian kernel
    void gaussKernel();
    //calculate local density with cut-off kernel
    void cutOffKernel();
    //sort data point according to their local density
    void sortRho();
    //calculate delta for each data point
    void calDelta();
    //update delta for each data point
    void updateDelta();
    //calculate gamma for each data point
    void calGamma();
    //pick the cluster centers according to gamma value
    void pickCenter();
    //assign cluster label for each data point
    void assignation();
    //find the corresponding position of dc in the dist list
    void findDcPos(vector<double> percents);
    //identify halo part of each cluster
    void IdentifyHalo();
    //set the number of clusters
    void setK(int K);
    //print cluster assignment for each element without identify cluster halo
    void printClusterAssignment();
    //output the decision graph, so that users can draw the graph use gnuplot
    void outputDecisionGraph();
};
#endif
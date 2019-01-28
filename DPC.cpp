//DPC Express implementation. DPC Express is a fast drop-in replacement of\
Density Peak Clustering via redundancy removal in Tuning\
@Copyright:© 2018 Shuai Yang
#include "DPC.h"

DPC::DPC(int n)
{ 
    N = n;
    dist = vector<vector<double>>(N, vector<double>(N, -1));
    rho = vector<double>(N, 0);
    delta = vector<double>(N, numeric_limits<double>::max());
    ndp = vector<int>(N);
    q = vector<int>(N);
    iota(q.begin(), q.end(), 0);
    gamma = vector<double>(N);
    ordgamma = vector<int>(N);
    nb = vector<vector<int>>(N);
    nbS = vector<vector<int>>(N);
    prevNDistLessDc = -1;
    nTrial = 0;

}

void DPC::IdentifyHalo()
{
    halo = cl;
    bord_rho = vector<double>(K, 0);

    if(K > 1)
    {
        double rho_avg;
        for(int i1 = 0; i1 < N; i1++)
        {
            int i = q[i1];
            if(delta[i] > dc)  //delta[i] > dc means there is no higher density data point in its neighbors, therefore, continue to next data point
                continue;
            
            int n = nb[i].size();
            for(int k= 0; k < n; k++)
            {
                int j = nb[i][k];
                if(cl[i] != cl[j] && rho[i] <= rho[j])
                {
                    rho_avg = (rho[i] + rho[j]) / 2;

                    if(rho_avg > bord_rho[cl[i]])
                        bord_rho[cl[i]] = rho_avg;
                    if(rho_avg > bord_rho[cl[j]])
                        bord_rho[cl[j]] = rho_avg;
                }  
            }
        }
        for(int i = 0; i < N; i++)
        {
            if(rho[i] < bord_rho[cl[i]]) 
                halo[i] = -1;
        }
    }
    /*if(printable)
    {
        cout<<"Clustering results for each data point:"<<endl;
        for(auto i : halo)
            cout<<i<<endl;
        cout<<endl<<endl;
    }*/
}

void DPC::printClusterAssignment()
{
    if (nTrial == 1)
        system("mkdir Result");
    ofstream myfile;
    myfile.open("Result/" + to_string(nTrial) + ".res");

    for (auto i : halo)
        myfile << i << endl;
    myfile.close();
}

void DPC::assignation()
{
    //assign each point to the same cluster as its nearest neighbour of higher density's
    for(int i = 0; i < N; i++)
    {
        if(cl[q[i]] == -1)
            cl[q[i]] = cl[ndp[q[i]]];
    }
}

void DPC::pickCenter()
{
    icl = vector<int>(K);
    cl = vector<int>(N, -1);
    for(int i = 0; i < K; i++)
    {
        icl[i] = ordgamma[i]; //store the index of the cluster center
        cl[ordgamma[i]] = i;  //label cluster center
    }
}

void DPC::setK(int K)
{
    this->K = K;
}

void DPC::outputDecisionGraph()
{
    if(nTrial == 1)
        system("mkdir DecisionGraph");

    ofstream myfile;
    myfile.open("DecisionGraph/DecisionGraph."+to_string(nTrial)+".dat");

    for(int i = 0; i < rho.size(); i++)
        myfile<<rho[i]<<" "<<delta[i]<<endl;
}
void DPC::calGamma()
{
    for(int i = 0; i < N; i++)
        gamma[i] = rho[i] * delta[i];
    iota(ordgamma.begin(), ordgamma.end(), 0);
    sort(ordgamma.begin(), ordgamma.end(), [&](const int &a, const int &b){return gamma[a] > gamma[b];}); //sort the index of data point according to gamma values
}

void DPC::calDelta()
{
    delta = vector<double>(N, numeric_limits<double>::max());
    ndp = vector<int>(N, numeric_limits<int>::max());
    delta[q[0]] = -1;
    ndp[q[0]] = 0; //the nearest data point of higher density of the highest density data point is undefined, just set it to 0.
    double maxDelta = numeric_limits<double>::min();
    
    for(int i = 1; i < N; i++)
    {
        int n = nb[q[i]].size();
        bool foundInDc = false;

        if(n <= q[i]) //if index of this point in the sorted list is smaller than its neighbor size, then use previous method
        {
            for(int k = 0; k < n; k++) //firstly check all its neighbors (distance < dc)
            {
                int j = nb[q[i]][k];

                if((rho[q[i]] < rho[j] || rho[q[i]] == rho[j] && j < q[i]) &&  //go through the neighbour list to check whether there is a new nearest neighbor of higher density
                (dist[q[i]][j] < delta[q[i]] || dist[q[i]][j] == delta[q[i]] && j < ndp[q[i]]))
                {
                    delta[q[i]] = dist[q[i]][j];
                    ndp[q[i]] = j;
                    foundInDc = true; //if found one in its neighbors, no need to explore non-neighbour.
                }
            }
        }

        if(!foundInDc) //if no new nearest neighbor of higher density is found in its neighbors or the index of this point in the sorted list is smaller than its neighbor size, then use previous method.
        {
            for(int j = 0; j < i; j++)
            {
                if(dist[q[i]][q[j]] < delta[q[i]] || dist[q[i]][q[j]] == delta[q[i]] && q[j] < ndp[q[i]])
                {
                    delta[q[i]] = dist[q[i]][q[j]];
                    ndp[q[i]] = q[j];
                }
            }
        }
        maxDelta = max(maxDelta, delta[q[i]]);
    }

    delta[q[0]] = maxDelta;
}


void DPC::updateDelta()
{
    vector<int> parStart(N, dcs.size()); //the start partition for each points
    //double maxDelta = numeric_limits<double>::min();

    for(int i = 1; i < N; i++)
    {
        int i1 = i;
        bool needCheckAll = false;
        while(i1 > 0 && (rho[q[i1]] > rho[q[i1-1]] || rho[q[i1]] == rho[q[i1 - 1]] && q[i1] < q[i1 - 1]))
        {
            q[i1] ^= q[i1 - 1];
            q[i1 - 1] ^= q[i1];
            q[i1] ^= q[i1 - 1];
            //update delta and nearest neighbor
            if(dist[q[i1]][q[i1 - 1]] < delta[q[i1]] || dist[q[i1]][q[i1 - 1]] == delta[q[i1]] && q[i1 - 1] < ndp[q[i1]])
            {
                delta[q[i1]] = dist[q[i1]][q[i1 - 1]];  
                ndp[q[i1]] = q[i1 - 1];
            }
            //if the one switched to the place behind is the ndp, then we need to find the new ndp
            if(ndp[q[i1 - 1]] == q[i1])
            {
                needCheckAll = true;
                
                for(int i2 = 0; i2 < dcs.size(); i2++)
                {
                    if(delta[q[i1 - 1]] < dcs[i2])
                    {
                        parStart[q[i1 - 1]] = i2; //determine which partition current ndp is in, then new ndp will not be less than it, so it searches from the beginning of this partition and ignore all prior partitions
                        break;
                    }
                }
            }
            i1--;
        }
        if(needCheckAll)
        {
            int n = nb[q[i1]].size();
            delta[q[i1]] = numeric_limits<double>::max();
            bool foundInDc = false;
            if(n <= q[i1] && parStart[q[i1]] < dcs.size()) //if adj size > the rank of the point in the sorted list or the ndp is outside all dc partition and locate in the remaining partition, then skip this if statement
            {   
                int prev; //prev record the index of the last point of the previous seg, next time we search from prev + 1
                if(parStart[q[i1]] == 0)
                    prev = -1;
                else
                    prev = nbS[q[i1]][parStart[q[i1]] - 1];

                for(int t = parStart[q[i1]]; t < nbS[q[i1]].size() && !foundInDc; t++)
                {
                    for(int k = prev + 1; k <= nbS[q[i1]][t]; k++)
                    {
                        int j = nb[q[i1]][k];

                        if((rho[q[i1]] < rho[j] || rho[q[i1]] == rho[j] && j < q[i1]) && 
                        (dist[q[i1]][j] < delta[q[i1]] || dist[q[i1]][j] == delta[q[i1]] && j < ndp[q[i1]]))
                        {
                            delta[q[i1]] = dist[q[i1]][j];
                            ndp[q[i1]] = j;
                            foundInDc = true;  //only if we find one in this partition, we know it's no need to iterate the next partition, the ndp must be in this partition
                        }
                    }
                    prev = nbS[q[i1]][t];
                }
            }
            if(!foundInDc) //get ndp using the standard way
            {
                for(int j = i1 - 1; j >= 0; j--)
                {
                    if(dist[q[j]][q[i1]] < delta[q[i1]] || dist[q[j]][q[i1]] == delta[q[i1]] && q[j] < ndp[q[i1]])
                    {
                        delta[q[i1]] = dist[q[j]][q[i1]];
                        ndp[q[i1]] = q[j];
                    }
                }
            }
        }
    }
    delta[q[0]] = numeric_limits<double>::min();
    ndp[q[0]] = 0;
    
    for(int i = 1; i < N; i++)
        delta[q[0]] = max(delta[q[0]], delta[q[i]]);

}



void DPC::sortRho()
{
    sort(q.begin(), q.end(), [&](const int &a, const int &b){
        if(rho[a] == rho[b])
            return a < b;
        return rho[a] > rho[b];
        });
}

void DPC::gaussKernel()
{
    for(int i = 0; i < N - 1; i++)
    {
        for(int j = i + 1; j < N; j++)
        {
            double delta_rho = exp(-(dist[i][j] / dc) * (dist[i][j] / dc));
            rho[i] += delta_rho;
            rho[j] += delta_rho;
        }
    }
}


void DPC::cutOffKernel()
{
    int position = round(percent / 100.0 * distList.size());

    dc = (distList.begin() + position)->first;
    //cout<<"dc:"<<dc<<endl;
    dcs.push_back(dc); //store new dc into wllocate list

    nDistLessDc = position;
    
    if(prevNDistLessDc == -1) //first dc
        prevNDistLessDc = 0;
    else if(distList[prevNDistLessDc].first == dc) //for special case, dc does not change
        return;
    //The omitted_dist is for the special case that some distances equal to dc and ranked in front of the corresponding position of dc.
    //They are not considered in the previous iteration but should be take into account for this iteration (new dc)
    for(int i = 0; i < omitted_dist.size(); i++)
    {
        int x = omitted_dist[i].second.first;
        int y = omitted_dist[i].second.second;
        rho[x]++;
        rho[y]++;
        nb[x].push_back(y);
        nb[y].push_back(x);
    }
            
    omitted_dist.clear();
    for(int i = prevNDistLessDc; i < nDistLessDc; i++)  //prevNDistLessDc indicates the first distances that have not been checked in previous iterations 
    {
        if(distList[i].first == dc)  //If the distance equal to dc, store it to ommited_dist, they will be used in the next iteration of dc
        {
            omitted_dist.push_back(distList[i]);
            continue;
        }
        int x = distList[i].second.first;
        int y = distList[i].second.second;
        
        rho[x]++;
        rho[y]++;
        nb[x].push_back(y);
        nb[y].push_back(x);
    }
    prevNDistLessDc = nDistLessDc;

    for(int i = 0; i < N; i++)
        nbS[i].push_back(nb[i].size() - 1);

}

void DPC::calRho(bool isGaussKernel)
{
    nTrial++;
    if(isGaussKernel)
        gaussKernel();
    else
        cutOffKernel();
}

//read distances from input file to dist and distList
//dist is the distance matrix and distList is the list of all distances
void DPC::readDist(string fileName)
{
    ifstream in(fileName);
    int i1,i2; //index of two data points

    while(!in.eof())
    {
        in>>i1>>i2;
        in>>dist[i1 - 1][i2 - 1];
        dist[i2 - 1][i1 - 1] = dist[i1 - 1][i2 - 1];
        distList.push_back(make_pair(dist[i2 - 1][i1 - 1], make_pair(i2 - 1, i1 - 1)));
    }
}

void DPC::findDcPos(vector<double> percents)
{
    auto endIter = distList.end();
    //find all dc positions from the largest dc to the smallest
    for(int i = percents.size() - 1; i >= 0; i--)
    {
        int position = round(percents[i] / 100.0 * distList.size());
        nth_element(distList.begin(), distList.begin() + position, endIter);
        endIter = distList.begin() + position;
    }
}

void DPC::setPercent(double Pct)
{
    percent = Pct;
}

double DPC::getDC()
{
    return dc;
}
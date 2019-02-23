#include "DPC.h"

DPC::DPC(int n)
{ 
    N = n;
    dist = vector<vector<double>>(N, vector<double>(N, -1));

    gamma = vector<double>(N);
    ordgamma = vector<int>(N);

}

void DPC::findHalo()
{
    halo = cl;
    bord_rho = vector<double>(K, 0);
    int count = 0;
    if(K > 1)
    {
        double rho_aver;
        for(int i = 0; i < N - 1; i++)
        {
            for(int j = i + 1; j < N; j++)
            {
                count++;
                if(cl[i] != cl[j] && dist[i][j] < dc)
                {
                    rho_aver = (rho[i] + rho[j]) / 2;

                    if(rho_aver > bord_rho[cl[i]])
                        bord_rho[cl[i]] = rho_aver;
                    if(rho_aver > bord_rho[cl[j]])
                        bord_rho[cl[j]] = rho_aver;
                }
            }
        } 
        for(int i = 0; i < N; i++)
        {
            if(rho[i] < bord_rho[cl[i]])
                halo[i] = -1;
        }
    }
    if(printable)
    {
        cout<<"Clustering results for each data point:"<<endl;
        for(auto i : halo)
            cout<<i<<endl;
        cout<<endl<<endl;
    }
}

void DPC::assignation()
{
    for(int i = 0; i < N; i++)
    {
        if(cl[ordrho[i]] == -1)
            cl[ordrho[i]] = cl[nneigh[ordrho[i]]];
    }
}

void DPC::setK(int K)
{
    this->K = K;
}

void DPC::setPrintOption(bool para)
{
    printable = para;
}

void DPC::pickCenter()
{
    this->K = K;
    icl = vector<int>(K);
    cl = vector<int>(N, -1);
    for(int i = 0; i < K; i++)
    {
        icl[i] = ordgamma[i];
        cl[ordgamma[i]] = i;
    }
}

void DPC::calGamma()
{
    for(int i = 0; i < N; i++)
        gamma[i] = rho[i] * delta[i];
    iota(ordgamma.begin(), ordgamma.end(), 0);
    sort(ordgamma.begin(), ordgamma.end(), [&](const int &a, const int &b){return gamma[a] > gamma[b];});
    /*for(auto i:ordgamma)
        cout<<i<<endl;
    cout<<endl<<"ttttt"<<endl;*/
}

void DPC::calDelta()
{
    delta = vector<double>(N, numeric_limits<double>::max());
    nneigh = vector<int>(N, numeric_limits<int>::max());

    delta[ordrho[0]] = -1;
    nneigh[ordrho[0]] = 0;
    double maxDelta = numeric_limits<double>::min();

    for(int i = 1; i < N; i++)
    {
        for(int j = 0; j < i; j++)
        {
            if(dist[ordrho[i]][ordrho[j]] < delta[ordrho[i]] || dist[ordrho[i]][ordrho[j]] == delta[ordrho[i]] && ordrho[j] < nneigh[ordrho[i]])
            {
                delta[ordrho[i]] = dist[ordrho[i]][ordrho[j]];
                nneigh[ordrho[i]] = ordrho[j];
            }
        }
        maxDelta = max(maxDelta, delta[ordrho[i]]);
    }

    delta[ordrho[0]] = maxDelta;
    
    /*for(auto i : delta)
        cout<<i<<endl;
    cout<<endl<<endl<<endl<<"xxxxx";*/
}

void DPC::sortRho()
{
    ordrho = vector<int>(N);
    iota(ordrho.begin(), ordrho.end(), 0);
    sort(ordrho.begin(), ordrho.end(), [&](const int &a, const int &b){
        if(rho[a] == rho[b])
            return a < b;
        return rho[a] > rho[b];
        });
    /*for(auto i:ordrho)
        cout<<i<<endl;
    cout<<"yyyyyyyyyyyy"<<endl<<endl;*/
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

bool flag3  = true;
void DPC::cutOffKernel()
{
    if(flag3)
    {
        timespec t_start, t_end;
	    clock_gettime(CLOCK_MONOTONIC, &t_start);
        sort(dist2.begin(), dist2.end());
        flag3 = false;
        clock_gettime(CLOCK_MONOTONIC, &t_end);
	    long timedif = 1000000*(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_nsec-t_start.tv_nsec)/1000;
        printf("Dc took %lf ms\n", timedif / 1000.0);
    }
    int position = round(percent / 100.0 * dist2.size());
    
    dc = *(dist2.begin() + position);

    rho = vector<double>(N, 0);

    for(int i = 0; i < N - 1; i++)
    {
        for(int j = i + 1; j < N; j++)
        {
            if(dist[i][j] < dc)
            {
                rho[i]++;
                rho[j]++;
            }
        }
    }
    /*for(auto i:rho)
        cout<<i<<endl;
    cout<<endl<<"rrrr"<<endl;*/
}

void DPC::calRho(bool isGaussKernel)
{
    if(isGaussKernel)
        gaussKernel();
    else
        cutOffKernel();
}

void DPC::readDist(string fileName)
{
    ifstream in(fileName);
    int i1,i2; //index of two data points
    
    while(!in.eof())
    {
        if(!(in>>i1>>i2))
	    break;
        in>>dist[i1 - 1][i2 - 1];
        dist[i2 - 1][i1 - 1] = dist[i1 - 1][i2 - 1];
        dist2.push_back(dist[i2 - 1][i1 - 1]);
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

/* process benchmark data to get distances between pairs
*/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>


using namespace std;

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        cout<<"Error! Usage: ./processBenchmark [input file name] [number of dimensions of data]"<<endl;
        return 0;
    }
    string inFileName(argv[1]), outFileName = inFileName + ".data";

    ifstream in(inFileName);
    vector<vector<double>> points;

    int n = atoi(argv[2]);

    while(!in.eof())
    {
        points.push_back(vector<double>());
        double tmpNum;
        for(int i = 0 ; i < n; i++)
        {
            if (in >> tmpNum)
                points[points.size() - 1].push_back(tmpNum);
            else
            {
                points.pop_back();
                break;
            }
                
        }
    }
    in.close();

    ofstream out(outFileName);
    for(int i = 0; i < points.size() - 1; i++)
    {
        for(int j = i + 1; j < points.size(); j++)
        {
            double dist = 0;
            for(int k = 0; k < n; k++)
                dist += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);

            out<<i + 1<<" "<<j + 1<<" "<<
            sqrt(dist)<<endl;
        }
    }
    out.close();
    return 0;
}

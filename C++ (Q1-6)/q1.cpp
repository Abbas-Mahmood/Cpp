#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;


double rofx(double x){
    double r0=0.8*x;        //initial r value
    double r=r0;
    int count =0;       //variable to count number of iterations
    double e=1;          //initialise e equal to arbitrary value greater than tolerance (10^-12).The variable e                      represents difference between successive r values.
    while(e>pow(10,-12)){       //execute loop while e is greater than the tolerance.
        r = x-log(r0-1);       //use r0 to calculate next r value.
        e = abs(r-r0);
        r0 = r;
        count++;
        }
    cout<<"Number of iterations: "<<count<<endl;        //outputs number of iterations.
    return r;       //returns final value of r.
}

int main(){
    cout<<setprecision(13);
    double x=4;
    double R=rofx(x);   //inputting x=4 into rofx
    cout<< "r = "<<R <<endl;
    double check=(R+log(R-1));
    cout<<"This value of r satisfies the original equation with an error of "<<setprecision(2)<<check-4<<endl;
}



#include <iostream>
using namespace std;
#include <cmath>
#include <vector>
#include <iomanip>


double qp(double p)
{
    return p;
}
double pp(double q)
{
    return -q;
}

int main(){
    vector<double>t(1001);      //vector t to store time values; size 1001(from 0 to 100 in 0.1 intervals)
    t[0]=0;                     //initial time
    vector<double>q(1001);      //vector q to store position values; size 1001
    q[0]=0;                     //initial position
    vector<double>p(1001);      //vector p to store momentum values; size 1001
    p[0]=sqrt(2);               //initial momentum.
    vector<double>E(1001);      //vector E to store energy values; size 1001
    E[0]=1;                     //inital energy value.
    
    vector<double>e(1001);      //vector e to store error in energy values.
    e[0]=E[0]-E[0];
    double dt=0.1;              //time interval

    for(int i=0;i<1000;i++){        //loop to solve system of equations using the RK2 method.
        double u1=qp(p[i]);           double v1=pp(q[i]);
        double u2=qp(p[i]+dt*v1*0.5); double v2=pp(q[i]+0.5*dt*u1);
        q[i+1]=q[i]+dt*u2;
        p[i+1]=p[i]+dt*v2;
        E[i+1]=0.5*(pow(p[i+1],2)+pow(q[i+1],2));
        t[i+1]=t[i]+dt;
        e[i+1]=E[i+1]-E[0];
    }
    cout<<"t \t\t q(t) \t\t p(t) \t\t E(t) \t\t e(t)"<<endl;
    for(int i=0;i<=1000;i++){     //prints q(t),p(t),E(t) and e(t) for t=0,1,10,100.
        if(i>0&i<10)
            continue;
        else if(i>10&i<100)
            continue;
        else if(i>100&i<1000)
            continue;
        cout<<setprecision(6)<<setw(8)<<left<<t[i]<<setw(12)<<q[i]<<setw(12)<<p[i]<<setw(12)<<E[i]<<e[i]<<endl;
    }
}

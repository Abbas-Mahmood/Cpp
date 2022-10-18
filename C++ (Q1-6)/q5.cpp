#include <iostream>
#include <cmath>
#include<vector>
#include<iomanip>
#include<random>
using namespace std;
double IP(vector<double>u,vector<double>v){
    double P=0;
    for(int i=0; i<v.size(); i++)
    P += u[i]*v[i];
    return P;
}
double g(double x)
{
    return (1/(1+25*pow(x,2)));
}

int main(){
    double Exact=0.4*atan(5);       //exact integral value.
    cout<<"Exact integral = "<<setprecision(10)<<Exact<<"\n"<<endl;
    
    double N=63; double a=-1; double b=1;   //N represents number of intervals. a and b are the limits of intergation
    double dx=(b-a)/N;      //width of each interval
    double x=a;             //inital x value; takes of the lower limit.
    vector<double>f;          //vector to store f(x) values.
    for(int i=0;i<=N;i++)    // loop to assign f(x) values.
    {
        f.push_back(g(x));
        x+=dx;
    }
//part(a)
    vector<double> w(N+1,(b-a)/N);      //vector w; represents weights for the trapezium rule; size is N+1.
    w[0]=(b-a)/(2*N);       //initial weight
    w[N]=(b-a)/(2*N);       //final weight
    double TrapeziumInt=IP(f,w);        //Trapezium integration calculation
    cout<<"Trapezium integration = "<<TrapeziumInt<<endl;
    cout<<"Trapezium error = "<<TrapeziumInt-Exact<<"\n"<<endl;
    
//part(b)
    double fpa=-50*a*pow((1+25*pow(a,2)),-2);       //f'(a) value
    double fpb=-50*b*pow((1+25*pow(b,2)),-2);       //f'(b) value
    double HermiteInt=TrapeziumInt+(pow(dx,2)/12)*(fpa-fpb);        //Hermite integration calculation.
    cout<<"Hermite integration = "<<HermiteInt<<endl;
    cout<<"Hermite error = "<<HermiteInt-Exact<<"\n"<<endl;
    
//part(c)
    vector<double>u;        // vector u to store f(x) values for C-C integration.
    for(int i=0;i<=N;i++)
    {
        u.push_back(g(-cos(M_PI*i/N)));
    }
    vector<double> Wcc;     //vector wcc; represents weights for the Clenshaw-Curtis rule.
    for(int i=0;i<=N;i++)       //loop to assign values to Wcc
    {
        double s=0;
        for(int k=1;k<=(N-1)/2;k++)
        {
            s+=2*cos(2*k*M_PI*i/N)/(4*pow(k,2)-1);
        }
        Wcc.push_back((2/N)*(1-s));
    }
    Wcc[0]=1/pow(N,2);      //initial weight
    Wcc[N]=1/pow(N,2);      //final weight
    double ClenshawCurtis=IP(u,Wcc);         //C-C integration calculation.
    cout<<"C-C integration = "<<ClenshawCurtis<<endl;
    cout<<"C-C error:"<<ClenshawCurtis-Exact<<"\n"<<endl;
//part(d)
    int s = 31;
    mt19937_64 mtrand(s);
    uniform_real_distribution<double> unif(0.0, 1.0);       //generate random values between 0 and 1

    int M = 10000;      //M represents number of samples
    int j = 0;          // variable to count number of random points below graph.
    for (int i = 0; i<M; i++) {
        double x = -1.0 + 2.0*unif(mtrand);     //random x values between -1 and 1
        double y = unif(mtrand);             // random y values between 0 and 1
        if (y<g(x))
            j++;
    }
    double MonteCarlo=2.0*double(j)/double(M);      //Monte-Carlo calculation.
    cout << "Monte Carlo = "<<MonteCarlo << endl;
    cout << "Monte Carlo error = "<<MonteCarlo-Exact << endl;
}
            

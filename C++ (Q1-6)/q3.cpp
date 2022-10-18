#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double IP(vector<double>u,vector<double>v){
    double P=0;
    for(int i=0; i<v.size(); i++)
    P += u[i]*v[i];
    return P;
}

double lm(vector<double>u,vector<double>v,int m){
    double P=0;
    for(int i=0; i<v.size(); i++)
        P += pow(abs(u[i]*v[i]), float(m)/2);
    return pow(P, 1/float(m));
}

vector<double> e(double j){
    vector <double> x;       //initialise vector to store x values
    vector <double> f;       //initialise vector to store f(x) values.
    vector <double> dfa;     //initialise vector to store analytically evaluated f'(x) values.
    vector <double> dfn;     //initialise vector to store numerically evaluated f'(x) values.
    vector <double> error;   //initialise vector to store error values.
    double N=pow(2,j)-1;     //Number of strips: 2^4-1=15 for part(a), N=2^j-1 is convenient for part(b).
    double dx=2/N;           //width of each strip
    for(int i=0;i<=N; i++){     //loop to assign x,f and dfa values.
        x.push_back((2*i-N)/N);
        f.push_back(pow(exp(1),-pow(x[i],2)));
        dfa.push_back(-2*x[i]*f[i]);
    }
    for(int j=0;j<=N; j++){     //loop to assign dfn values.
        if(j==0)
            dfn.push_back((-3*f[j]+4*f[j+1]-f[j+2])/(2*dx));
        else if(j<N)
            dfn.push_back((f[j+1]-f[j-1])/(2*dx));
        else if(j==N)
            dfn.push_back((f[j-2]-4*f[j-1]+3*f[j])/(2*dx));
    }
    for(int k=0;k<=N; k++){     //loop to assign error values. error=dfn-dfa.
    error.push_back(dfn[k]-dfa[k]);
    }
    return error;
}

int main(){
//part(a)
    vector<double>b;
    cout<<"e[i]"<<endl;
    b=e(4);
    for(int i=0;i<b.size();i++){
        cout<<b[i]<<endl;
    }
//part(b)
    vector<double>C;    //initalise vector C to examine 2nd order convergence.
    vector<double>a;    //initalise vector to store error for different N=2^j-1 values.
    cout<<"\n C[i] "<<endl;
    for(int j=4;j<=14;j++){     //loop to assign C values
        a=e(j);
        int N = pow(2,j)-1;
        C.push_back(pow(N,1.5)*lm(a,a,2));
    }
    for(int i=0;i<C.size();i++){        //loop to output C values.
        cout<<C[i]<<endl;
    }
}




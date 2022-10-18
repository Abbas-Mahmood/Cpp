#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

double IP(vector<double>u,vector<double>v){      //function takes two vectors and returns their inner product.
    double P=0;     //initialise P to store the value of the inner product.
    for(int i=0; i<v.size(); i++)   //calculate the inner product by summing the product of u[i] and v[i] terms.
    P += u[i]*v[i];
    return P;       //returns inner product.
}

double lm(vector<double>u,vector<double>v,int m){    //function takes two vectors and returns the m-th weighted norm
    double P=0; //initialise P to store the value of the (weighted norm)^m.
    for(int i=0; i<v.size(); i++)
        P += pow(abs(u[i]*v[i]), float(m)/2);
    return pow(P, 1/float(m));      //returns weighted norm.
}
int main(){
    vector<double>u{2,7,2};
    vector<double>v{3,1,4};
    double w=lm(u,v,1);
    double z=lm(u,v,2);
    cout<<"IP(u,v) = "<<IP(u,v) <<endl;  //prints inner product of u and v
    cout<<"lm(u,v,1) = "<< w <<endl;    //prints weighted norm(1) of u and v
    cout<<"lm(u,v,2) = "<< z <<endl;    //prints weighted norm(2) of u and v
    if(z*z==IP(u,v))
    cout<< "lm(u,v,2)^2 = u*v" <<endl;  //prints weighted norm(2)^2 of u and v
}


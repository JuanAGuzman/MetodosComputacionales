#include<iostream>
#include <math.h>  /* sin */

using namespace std;

class Derivada{
    private:
        float h;
        float d = 0;
    public:
        void def_h(float h_);
        float central(float (*funcion)(float), float x);
};

float Derivada::central(float (*funcion)(float), float x){
    d = 0;
    if(h!=0){
        d = (funcion(x+h) - funcion(x-h)) / 2*h;
    };
    return d;
}

void Derivada::def_h(float h_){
    h = h_;
}


float f(float x){
    return x - sin(x);
}

int main(){
    Derivada derivada;
    derivada.def_h(0.05);
    cout<<'La derivada en 0 es '<< derivada.central(f, 0);
    return 0;
}
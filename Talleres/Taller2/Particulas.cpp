#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <random>

using namespace std;

class Particle{
	
	public:
		Particle();
		Particle(double x_,double y_, double vx_, double vy_, double m_, double r_, double ID_);
		~Particle();

		void SetWallLimits(double Wxmin_, double Wxmax_, double Wymin_, double Wymax_);
		void CheckWallLimits();
		void Print();

		void Move(double t_, double deltat, int it);

		void ChangeWallLimits(double deltat, double v, double lim_);

	private:

		double x,y;
		double t;

		double vx,vy;
		double ax,ay;
		double Fx,Fy;

		double m,r;

		int ID;

		std::ofstream *File;

		double Wxmin,Wxmax;
		double Wymin,Wymax;

		double K = 100.;

	protected:


};


Particle::Particle(){
}

Particle::Particle(double x_,double y_, double vx_, double vy_, double m_, double r_, double ID_): x(x_), y(y_), vx(vx_), vy(vy_), m(m_), r(r_), ID(ID_){
 
}

Particle::~Particle(){
}

// Methods

void Particle::SetWallLimits(double Wxmin_, double Wxmax_, double Wymin_, double Wymax_){
Wxmin = Wxmin_;
Wxmax = Wxmax_;
Wymin = Wymin_;
Wymax = Wymax_;
}

void Particle::ChangeWallLimits(double deltat, double v, double lim_){

	if(Wxmax >= lim_)
		Wxmax += v*deltat;
}

void Particle::CheckWallLimits(){

	if( (x+r) >= Wxmax || (x-r) <= Wxmin ) vx = -vx;
	if( (y+r) >= Wymax || (y-r) <= Wymin ) vy = -vy;


}

void Particle::Print(){
 std::cout <<" , "<<x<<"+"<<r<<"*cos(t),"<<y<<"+"<<r<<"*sin(t)";
}

void Particle::Move(double t_, double deltat, int it){


	t = t_;
	
	x += vx*deltat;
	y += vy*deltat;

}


void StartAnim(){

  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange [-50:150]"<<std::endl;
  std::cout<<"set yrange [-50:50]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}

void StartLine(){
        std::cout <<"plot 0,0";
}

void EndLine(){
        std::cout<<std::endl;
}


int main(int argc, char *argv[]){

	Particle *p1 = new Particle(2,3,2,3.3,5,6,0);
	p1->SetWallLimits(-50,50,-50,50);

	// Evolucion 
	double deltat = 0.0001;

	double tmax = std::stof(argv[1]);
	int it = 0;

        int films = 5000;

        
	int NParticles = 50;

	double time = 0.;

	Particle *AllParticles[NParticles];

	for(int i = 0; i< NParticles; i++){
		Particle *p = new Particle(50*rand(),50*rand(), 50.*rand(), 50*rand(),5,6,i);
		p->SetWallLimits(-50,50,-50,50);
		AllParticles[i] = p;
	}


	// Evolucion
	while (time < tmax){

        if(it%films == 0){
	StartAnim();
	StartLine();
	}

	for(int i = 0; i< NParticles; i++){
	 AllParticles[i]->ChangeWallLimits(deltat,-1.0,0.);
	}

	for(int i = 0; i< NParticles; i++){
		AllParticles[i]->Move(time,deltat,it);
		AllParticles[i]->CheckWallLimits();
	

	//p1->Move(time,deltat,it);
	//p1->CheckWallLimits();

        if(it%films == 0)
	{
	//p1->Print();
	
	AllParticles[i]->Print();
	}

	}

	time += deltat;
	it ++;
        if(it%films == 0)	
	EndLine();
	}
	
	return 0;

}
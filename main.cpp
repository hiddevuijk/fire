
#include "fire2.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class E {
public:
	E(double a, double b, double c,double d):a(a),b(b),c(c),d(d) {}
	double a,b,c,d;

	void dE( vector<double>& F, const vector<double>& x)
	{
		F[0] = -2*c*(x[0] - a);
		F[1] = -2*d*(x[1] - b);
	}

    double energy(const vector<double>& x)
    {
        double e = 0;
        e += c*(x[0]-a)*(x[0]-a);
        e += d*(x[1]-b)*(x[1]-b);
        return e;
    }

	void EdE(double& U, vector<double>& F, const vector<double>& x)
	{

		F[0] = -2*c*(x[0] - a);
		F[1] = -2*d*(x[1] - b);
		U = 0;	
        U += c*(x[0]-a)*(x[0]-a);
        U += d*(x[1]-b)*(x[1]-b);
	}
};

int main()
{
	E energy(1,2,1,1);

	Fire<E> fire(2, &energy);
	fire.error = 1e-12;
	fire.dt0 = 0.01;

	fire.error = 1e-8;
	fire.error_max = 1e-9;

	fire.dt0 = 1.e-2;
	fire.dtmax = 100;
	fire.dtmin = 1e-5;
	fire.finc = 1.1;
	fire.fdec = 0.5;
	fire.Nmin = 5;
	fire.alpha0 = 0.1;
	fire.falpha = 0.99;
	fire.m = 1.;

	vector<double> x(2);
	x[0] = 100;
	x[1] = 0;
	//fire.minimizeEE(x);
	fire.minimizeSIE(x);
	//fire.minimizeVV(x);

	cout << fire.x[0] << '\t' << fire.x[1] << "\t" << fire.Fnorm << endl;
    cout << energy.energy(fire.x) << endl;
	cout << fire.Fnorm << endl;
	cout << fire.NF <<endl;


	return 0;
}


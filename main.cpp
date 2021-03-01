
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
};

int main()
{
	E energy(1,2,1,1);

	Fire<E> fire(2, &energy);
	fire.error = 1e-12;

	vector<double> x(2);
	x[0] = 1;
	x[1] = 1;
	fire.minimizeSIE(x);

	cout << fire.x[0] << '\t' << fire.x[1] << endl;
    cout << energy.energy(fire.x) << endl;
	cout << fire.Fnorm << endl;


	return 0;
}


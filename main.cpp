
#include "fire.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class Parameters {
public:
	Parameters( double a, double b, double c)
		: a(a), b(b), c(c) {}
	double a,b,c;
};

void dE( vector<double>& F,
		const vector<double>& x, void *params)
{

	
	//F[0] = -p->a - 2*(p->c)*x[0];
	//F[1] = -p->b - 2*(p->c)*x[1];
	double a = 1;
	double b = 2;
	double c = 3;
	F[0] = -a - 2*c*x[0];
	F[1] = -b - 2*c*x[1];


}


class E {
public:
	E(double a, double b, double c,double d):a(a),b(b),c(c),d(d) {}
	double a,b,c,d;

	void dE( vector<double>& F, const vector<double>& x)
	{
		F[0] = -2*c*(x[0] - a);
		F[1] = -2*d*(x[1] - b);
	}
};

int main()
{
	E energy(1,2,10,20);

	Fire<E> fire(2, &energy);
	fire.error = 1e-6;

	vector<double> x(2);
	x[0] = 1;
	x[1] = 1;
	fire.minimize(x);

	cout << fire.x[0] << endl << fire.x[1] << endl;
	cout << fire.Fnorm << endl;
	cout << fire.dt << endl;


	return 0;
}


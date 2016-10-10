#include <cstdlib>
#include <iostream>
#include <cmath>
using namespace std;

int main(){
	int n, nstop=30;
	double db, p;

	for (n=5; n<=nstop; n++){
		db=2*M_PI/n;
		p=1.0-sin(db)/db;
		cout << n << " " << p;
		db=db/2;
		p=1.0-sin(db)/db;
		cout << " " << p << endl;
	}

	return 0;
}

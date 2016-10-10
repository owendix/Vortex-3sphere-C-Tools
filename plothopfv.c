#include "linalg.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
using namespace std;
/*************************************************
Calculate trajectory due to hopf vector field: v=(-y,x,-w,z)
Projects that R^4 trajectory to either hopf or hyperspherical 
coords: method=0,1 respectively
**************************************************
Input: g++ -g plothopfv.c linalg.c -o hopfv.out
		./hopfv.out > tmp.dat
Output: (xx1):
		out0 out1 out2	[x(steps+1)]
		...
		(xx2):
		out0 out1 out2 	[x(steps+1)]
		...

		repeated n times

	-> splits into xx1 xx2 ... xxn
**************************************************/

void hopffib(point, double*);
void hspher(point, double*);


int main(){
	double d, r[4], dt, out[3];
	int n=80, i, j, k, s, limhigh=1, limlow=0, steps=1000, fac=1;
	point o, q, onew, qnew, pnew, p, pold, v;
	ofstream fout;
	string fbase="xx", fname, strnum;
	stringstream strhandle;

	int method=0;	//=0, hopffib; =1, hspher

	srand((unsigned)time(0));

	dt=2*M_PI/steps;	//make it just close if its a great circle

	//for all trajectories on 3sphere
	for (i=1;i<=n;i++){
		strhandle << i;
		strnum = strhandle.str();
		fname=fbase+strnum;
		strhandle.str("");
		fout.open(fname.c_str());
		//Generate initial value: r of trajectory
		d=2.0;
		while(d>1.0){
			s = (rand()%(limhigh+1))+limlow;//rand ints fm limlow->high
			s = 2*s-1;	//s = 1 or -1
			for (j=0;j<3;j++){
				//can't seed it right before: produces bad numbers 
    			d=((double)rand()+(double)rand()/(double)RAND_MAX)/
					((double)RAND_MAX);
    			//random nums 0<=d<=1: added factor is to get finer spacing
    			//otherwise won't produced numbers between 0 and 1/RAND_MAX
			 	r[j]=2.0*d-1.0;
			}
			d=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
		}
		r[3]=s*sqrt(1-d);
		//first pt
		q.x=r[0]; q.y=r[1]; q.z=r[2]; q.w=r[3];
		o=q;	//Other methods of doing numerical integration
		p=q;
		fout.precision(15);
		
		/*
		//1st method
		if (!method)	//project: hopf fibration, hyperspherical
			hopffib(q, out);
		else
			hspher(q, out);
		for (k=0;k<3;k++){
			fout << out[k];
			if (k!=2)
				fout << " ";
			else
				fout << endl;
		}
		for (j=0;j<fac*steps;j++){
			//calculate velocity, v
			v.x=-q.y; v.y=q.x; v.z=-q.w; v.w=q.z;
			//q method: uses geom-series only (fm expansion): exact?
			qnew=(q+v*dt)/(1+dt*dt);
			q=qnew;
			if (!method)	//project: hopf fibration, hyperspherical
				hopffib(q, out);
			else
				hspher(q, out);
			for (k=0;k<3;k++){
				fout << out[k];
				if (k!=2)
					fout << " ";
				else
					fout << endl;
			}
		}
		*/
		//2nd method: more accurate than method 1
		if (!method)	//project: hopf fibration, hyperspherical
			hopffib(o, out);
		else
			hspher(o, out);
		for (k=0;k<3;k++){
			fout << out[k];
			if (k!=2)
				fout << " ";
			else
				fout << endl;
		}
		for (j=0;j<fac*steps;j++){
			//calculate velocity, v
			v.x=-o.y; v.y=o.x; v.z=-o.w; v.w=o.z;
			//o method: uses sqrt only (fm geometry): exact?
			onew=o*sqrt((1-dt)*(1+dt))+v*dt;
			o=onew;
			if (!method)	//project: hopf fibration, hyperspherical
				hopffib(o, out);
			else
				hspher(o, out);
			for (k=0;k<3;k++){
				fout << out[k];
				if (k!=2)
					fout << " ";
				else
					fout << endl;
			}
		}
		/*
		//3rd method (approx): about as accurate as 2nd method
		if (!method)	//project: hopf fibration, hyperspherical
			hopffib(p, out);
		else
			hspher(p, out);
		for (k=0;k<3;k++){
			fout << out[k];
			if (k!=2)
				fout << " ";
			else
				fout << endl;
		}
		//calculate velocity, v
		v.x=-p.y; v.y=p.x; v.z=-p.w; v.w=p.z;
		//p method: uses sqrt to get 2nd pt, then approx
		pold=p;
		//need 2nd pt
		p = pold*sqrt((1-dt)*(1+dt)) + v*dt;
		if (!method)	//project: hopf fibration, hyperspherical
			hopffib(p, out);
		else
			hspher(p, out);
		for (k=0;k<3;k++){
			fout << out[k];
			if (k!=2)
				fout << " ";
			else
				fout << endl;
		}
		for (j=0;j<(fac*steps-1);j++){
			//calculate v at new pt
			v.x=-p.y; v.y=p.x; v.z=-p.w; v.w=p.z;
			pnew = pold + v*(2*dt);
			pold=p;
			p=pnew;
			if (!method)	//project: hopf fibration, hyperspherical
				hopffib(p, out);
			else
				hspher(p, out);
			for (k=0;k<3;k++){
				fout << out[k];
				if (k!=2)
					fout << " ";
				else
					fout << endl;
			}
		}
		*/
		fout.close();
	}
	return 0;
}

void hopffib(point r, double* h){
	//projects great circles onto points on 2sphere
	//valid for all r0
	h[0]=2*(r.x*r.z+r.y*r.w);
	h[1]=2*(r.y*r.z-r.x*r.w);
	h[2]=(r.x-r.z)*(r.x+r.z)+(r.y-r.w)*(r.y+r.w);

	//projected onto 2sphere of radius: r0*r0
	return;
}

void hspher(point r, double* h){
	//hyperspherical coordinates
	h[0]=atan2(sqrt(r.w*r.w+r.z*r.z+r.y*r.y),r.x);
	h[1]=atan2(sqrt(r.w*r.w+r.z*r.z),r.y);
	h[2]=atan2(r.w,r.z);

	return;
}

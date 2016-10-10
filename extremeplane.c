#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include "functextrem.h"
#ifndef LIN_ALG
	#define LIN_ALG
	#include "linalg.h"
#endif
#include <vector>
#include <list>
using namespace std;

//#define HESSIANCHECK
/*******************************************************
 * Finds minimum and maximum angles between planes
 * *****************************************************
 * Compile with linalg.c
 * functextrem.h and linalg.c include linalg.h
 * *****************************************************
 * Input:
 * 		./a.out plnes/plane.filename.#.dat
 * Output:
 * 		minang    maxang    xxyy
 *
 * 	xx= plane 1 number
 * 	yy= plane 2 number
 *******************************************************
 * */
struct angplane{
	double min;
	double max;
	int pnum;
};

list<angplane>::iterator sortinsert(const double, list<angplane> *);

int main(int argc, char * argv[]){
	ifstream fplanein;
	plane tmpplane;
	point tmppt;
	double y1, y2, y3, y4, jnk=1.0, TINY=1e-14;
	double term1, term2, tmpd;
	double x, y, x0, y0, fmax, fmin;
	int ycnt=0, inserted=0;
	angplane tmpangplane;
	#ifdef HESSIANCHECK
	double det, a, b, eigval1, eigval2, fextrem;
	#endif

	if (argc==2)
		fplanein.open(argv[1]);
	else
		cout << "Include: filename";
	vector<plane> vecplanes;
	vector<plane>::iterator vj=vecplanes.begin();
	//read in planes data (spanning vectors)
	while(1){
		fplanein >> tmpplane.pnum;
		if (fplanein.eof()) break;
		fplanein >> tmpplane.p1;
		fplanein >> tmpplane.p2;
		fplanein >> tmppt;	//Orthogonal vectors
		fplanein >> tmppt;
		vecplanes.push_back(tmpplane);
		++vj;
	}
	fplanein.close();
	
	//for all nonidentical permutations n(n-1)/2
	vector<plane>::iterator vi=vecplanes.begin();
	vj=vecplanes.begin();
	++vj;
	//for storing/sorting
	list<angplane> anglelist;
	list<angplane>::iterator ilist;
	int count=0;
	while(1){
		//solve for extremal vals: 2 dim (x=theta,y=phi)
		y1=dot((*vi).p2,(*vj).p1)+dot((*vi).p1,(*vj).p2);
		y2=dot((*vi).p1,(*vj).p1)-dot((*vi).p2,(*vj).p2);
		y3=dot((*vi).p2,(*vj).p1)-dot((*vi).p1,(*vj).p2);
		y4=dot((*vi).p1,(*vj).p1)+dot((*vi).p2,(*vj).p2);
		if (abs(y2)>TINY && abs(y4)>TINY){
			term1=atan(y1/y2);
			term2=atan(y3/y4);	//Base values
		}else{ 
			if (abs(y2) < TINY && abs(y4) < TINY){
				jnk=1.0;
				if (y1*y2<0) jnk=-1.0;
				term1=jnk*M_PI-atan(y2/y1);
				jnk=1.0;
				if (y3*y4<0) jnk=-1.0;
				term2=jnk*M_PI-atan(y4/y3);
			}else if (abs(y2) < TINY){
				jnk=1.0;
				if (y1*y2<0) jnk=-1.0;
				term1=jnk*M_PI-atan(y2/y1);
				term2=atan(y3/y4);
			}else{
				term1=atan(y1/y2);
				jnk=1.0;
				if (y3*y4<0) jnk=-1.0;
				term2=jnk*M_PI-atan(y4/y3);
			}
		}
		//determine type of extremum: f-value 
		x=x0=0.5*(term1+term2);
		y=y0=0.5*(term1-term2);
		fmax=-1.0;
		fmin=1.0;
		ycnt=0;
		//establish designation code ppqq
		tmpangplane.pnum=100*(*vi).pnum+(*vj).pnum;
	
		Functgh f=Functgh((*vi).p1, (*vi).p2, (*vj).p1, (*vj).p2);
		//Go through the 9 permutations of (x,y) combos
		for (int i=0;i<9;i++){
			//test function values for max and min
			if ((tmpd=f(x,y))>fmax) fmax=tmpd;
			if (tmpd<fmin) fmin=tmpd;
		#ifdef HESSIANCHECK
			//Check determinant of hessian to look for extrema
			a=-tmpd;	//Equals f(x,y)
			b=f.dxdy(x,y);
			eigval1=a-b;
			eigval2=a+b;	
			det=eigval1*eigval2;
			if (det>TINY){
				fextrem=tmpd;
			}
			cout.setf(ios_base::showpoint);
			cout << det << " " << acos(tmpd) << " ";
			cout.width(4);
			cout.fill('0');
			cout << tmpangplane.pnum << endl;
			cout.fill(' ');
			//cout << det << " " << eigval1 << " " << eigval2 << " ";
			//cout << fmax << " " << fmin << " " << tmpangplane.pnum << endl;
		#endif
			x+=M_PI/2.0;
			if (i%3==2){ 
				x=x0;
				ycnt++;
				if (ycnt==1) y=y0-M_PI/2.0;
				if (ycnt==2) y=y0+M_PI/2.0;
			}
		}
		fmax=acos(fmax);
		fmin=acos(fmin);
		//doubly linked list, sort by fmax, while inserting
		tmpangplane.min=fmax;	//max dotprod is min angle
		tmpangplane.max=fmin;
	#ifdef HESSIANCHECK
		//prints min angle max angle and designation ####
		cout.setf(ios_base::showpoint);
		cout << fmax << " " << fmin << " ";
		cout.width(4);
		cout.fill('0');
		cout << tmpangplane.pnum << endl;
		cout.fill(' ');
	#else
		if (!count){ 
			anglelist.push_back(tmpangplane);
		}else{
			//Build list of angplane structures: anglelist
			inserted=0;
			//Check beginning
			ilist=anglelist.begin();
			if (tmpangplane.min<(*ilist).min){
				anglelist.push_front(tmpangplane);
				inserted=1;
			}
			//Check end
			ilist=anglelist.end();
			ilist--;
			if (!inserted && tmpangplane.min>(*ilist).min){
				anglelist.push_back(tmpangplane);
				inserted=1;
			}
			//Insert somewhere in between begin() and end()
			if (!inserted)
				anglelist.insert(sortinsert(tmpangplane.min, &anglelist), 
					tmpangplane);
		}
	#endif
		vj++;
		if (vj==vecplanes.end()){
			++vi;
			vj=vi+1;
		}
		if (vj==vecplanes.end())
			break;
		count++;
	}
#ifndef HESSIANCHECK
	cout.setf(ios_base::showpoint);
	cout.precision(15);
	//Need to print
	for (ilist=anglelist.begin();
			ilist!=anglelist.end(); ilist++){
		cout << (*ilist).min << " " << (*ilist).max << " ";
		cout.width(4);
		cout.fill('0');
		cout << (*ilist).pnum << endl;
		cout.fill(' ');
	}
#endif
}
//sorts as it inserts: returns pointer to location for insert in list
//divide and conquer
list<angplane>::iterator sortinsert(const double testmin, 
										list<angplane> * angs){
	//Guaranteed upon entering: begin < testmin < end
	//Hence, also guaranteed upon entering: n >=2
	
	//Initialize limits
	list<angplane>::iterator low=(*angs).begin();
	list<angplane>::iterator high=(*angs).end();
	--high; 		//end() points to location SUCCEEDING end of list
	list<angplane>::iterator mid;
	int n=(*angs).size(), i;	//how many between low and high (inclusively)
	//Limits: guarantees low < testmin < high
	while(n>2){
		//Reassign mid
		//Works for n even or odd
		for (i=0, mid=low; i<(n/2); i++, mid++);	
		//Check spot in middle of range
		if (testmin<(*mid).min){
			high=mid;
			n=n/2+1;	//Works for n even or odd
		}else{
			low=mid;
			if (n%2==0) //even
				n=n/2;	
			else		//odd 
				n=n/2+1;
		}
	}
	//Impossible to get n<2
	return high;	//Return pointer to one greater
}

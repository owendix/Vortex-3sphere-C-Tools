#include "functors.h"
#include "quadrature.h"
#include <vector>
#include <map>
#include <fstream>

//#define CHECKMINNUM		//uncomment to check minimum by brute force testing
//#define CHECKMINVIS	//uncomment to check minimum by printing values
//#define CHECK_ORTHOGONAL	//uncomment prints orthogonal.out if planes orthog
//#define CHECKDIRVIS
//#define CHECKDIRNUM

//Finds area between two vortices
//------------------------------------------------------
//Compile with linalg.c
//functors.h and linalg.c include linalg.h
//------------------------------------------------------
//Input:
//		./a.out plns/plane.filename.#.dat
//
//Output:
//		area_val	xxyy
//
//	xx = vortex 1 number
//	yy = vortex 2 number
//-------------------------------------------------------


int main(int argc, char *argv[]){
	ifstream fplanein;
	int pm, jnk=1;
	point tmppt;
	double x, xmin, area, r0=.05, TINY=1e-3;
	plane tempplane;
	if (argc==2)
		fplanein.open(argv[1]);
	else
		cout << "Include: filename\n";
#ifdef CHECK_ORTHOGONAL
	ofstream forthog;
#endif	
	//Map class
	map<double, int> mapplanes;	//Key value will be area: map sorts by key
	//Vector class
	vector<plane> vecplanes;
	//Read in data
	vector<plane>::iterator vj=vecplanes.begin();	//Vector iterator
	while(1){
		fplanein >> tempplane.pnum;
		if (fplanein.eof())
			break;
		fplanein >> tempplane.p1;
		fplanein >> tempplane.p2;
		fplanein >> tmppt;
		fplanein >> tmppt;
		vecplanes.push_back(tempplane);
		++vj;
	}
	fplanein.close();	
	//Go through all nonidentical permutations n(n-1)/2
	vector<plane>::iterator vi=vecplanes.begin();
	vj=vecplanes.begin();
	++vj;
	double y1, y2, y3, y4, x0=0.0;
	double xmin1, xmin2;	
	while (1){
		//Analytically minimize function
		y1=dot((*vj).p1,(*vi).p1);
		y2=dot((*vj).p2,(*vi).p1);
		if (abs(y1)<TINY && abs(y2)<TINY){
#ifdef CHECK_ORTHOGONAL
				if(!forthog.is_open())
					forthog.open("orthogonal.out");
				forthog << "Orthogonal vectors:" << endl;
				forthog << (*vi).pnum << " " << (*vj).pnum;
				forthog << ": {p1, q1}, {p1, q2}" << endl;
#endif
			//plane j orthogonal to plane i's vector: p1
			x0=M_PI/2.0;
			y3=dot((*vj).p1,(*vi).p2);
			y4=dot((*vj).p2,(*vi).p2);
			//eval f(M_PI/2.0), for both xmin1, xmin2 AND pm=1,-1
			if (abs(y3)<TINY && abs(y4)<TINY){
				//plane j orthogonal to plane i
				//distance unaffected by xmin (and theta)
#ifdef CHECK_ORTHOGONAL
				if(!forthog.is_open())
					forthog.open("orthogonal.out");
				forthog << "Orthogonal planes:" << endl;
				forthog << (*vi).pnum << " " << (*vj).pnum << endl;
#endif
				xmin=0.0;
			}else if (abs(y3)<TINY){
#ifdef CHECK_ORTHOGONAL
				if(!forthog.is_open())
					forthog.open("orthogonal.out");
				forthog << "Orthogonal vectors:" << endl;
				forthog << (*vi).pnum << " " << (*vj).pnum;
				forthog << ": {p1, q1}, {p1, q2}, {p2, q1}" << endl;
#endif
				jnk=1;
				if (y3*y4<0) jnk=-1;
				xmin1=jnk*M_PI/2.0-atan(y3/y4);
				xmin2=xmin1+M_PI;
				Funcd g=Funcd((*vi).p1,(*vi).p2,
						(*vj).p1,(*vj).p2, x0);
				xmin=(g(xmin1)<g(xmin2)) ? xmin1 : xmin2;
			}else{
#ifdef CHECK_ORTHOGONAL
				if(!forthog.is_open())
					forthog.open("orthogonal.out");
				forthog << "Orthogonal vectors:" << endl;
				forthog << (*vi).pnum << " " << (*vj).pnum;
				forthog << ": {p1, q1}, {p1, q2}" << endl;
#endif
				xmin1=atan(y4/y3);
				xmin2=xmin1+M_PI;
				Funcd g=Funcd((*vi).p1,(*vi).p2,
						(*vj).p1,(*vj).p2, x0);
				xmin=(g.d2x(xmin1)>0) ? xmin1 : xmin2; //check curvature
			}
		}else if(abs(y1)<TINY){
#ifdef CHECK_ORTHOGONAL
				if(!forthog.is_open())
					forthog.open("orthogonal.out");
				forthog << "Orthogonal vectors:" << endl;
				forthog << (*vi).pnum << " " << (*vj).pnum;
				forthog << ": {p1, q1}" << endl;
#endif
			x0=0.0;
			jnk=1;
			if (y1*y2<0) jnk=-1;
			xmin1=jnk*M_PI/2.0-atan(y1/y2);
			xmin2=xmin1+M_PI;
			Funcd g=Funcd((*vi).p1,(*vi).p2,(*vj).p1,(*vj).p2,x0);
			xmin=(g.d2x(xmin1)>0) ? xmin1 : xmin2; 	//check curvature
		}else{
			xmin1=atan(y2/y1);
			xmin2=xmin1+M_PI;
			Funcd g=Funcd((*vi).p1,(*vi).p2,(*vj).p1,(*vj).p2,x0);
			xmin=(g.d2x(xmin1)>0) ? xmin1 : xmin2;	//check curvature
		}

#ifdef CHECKMINNUM
		//Check minimum numerically (brute force)
		while (xmin>2.0*M_PI) xmin-=2.0*M_PI;
		while (xmin<0.0) xmin+=2.0*M_PI;
		double dx=M_PI/50.0, gmin;
		Funcd g1((*vi).p1, (*vi).p2, (*vj).p1, (*vj).p2, x0);
		gmin=g1(xmin);
		for (x=0.0;x<2.0*M_PI;x+=dx)
			if (g1(x)<gmin){
				cout << "Minimum not properly attained in vortices: ";
				cout << (*vi).pnum << ", " << (*vj).pnum << endl;
			}
		
#elif defined(CHECKMINVIS)
		//Check minimum by printing and visually inspecting
		while (xmin>2.0*M_PI) xmin-=2.0*M_PI;
		while (xmin<0.0) xmin+=2.0*M_PI;
		double dx=M_PI/50.0;
		Funcd g1((*vi).p1, (*vi).p2, (*vj).p1, (*vj).p2, x0);
		cout << xmin << " " << g1(xmin) << endl;
		for (x=0.0;x<2.0*M_PI;x+=dx)
			cout << x << " " << g1(x) << endl;
 		cout <<  "blah" << endl;
		
#else
		x=M_PI/10.0;
		pm=1;
		Funcdist f1 = Funcdist((*vi).p1, (*vi).p2, (*vj).p1, 
					(*vj).p2, x0, xmin, pm);
		pm=-1;
		Funcdist f2 = Funcdist((*vi).p1, (*vi).p2, (*vj).p1, 
					(*vj).p2, x0, xmin, pm);
		//Pick whichever direction minimizes the function
	#ifdef CHECKDIRVIS
		double xtmp, dxtmp=M_PI/50.0;
		int dir = f1(x)<f2(x) ? 1 : 2 ;
		for (xtmp=0.0;xtmp<2*M_PI;xtmp+=dxtmp){
			cout << xtmp << " " << f1(xtmp) << " " << f2(xtmp) <<  " ";
			cout << dir << endl;	
		}
		cout << "blah" << endl;
	#elif defined(CHECKDIRNUM)
		double area1, area2;
		area1=qsimp(f1, 0.0, 2.0*M_PI, 1.0e-8);
		area2=qsimp(f2, 0.0, 2.0*M_PI, 1.0e-8);
		if (area1<area2 && f1(x)>f2(x)){
			cout << "Wrong direction picked for vorts: ";
			cout << (*vi).pnum << ", " << (*vj).pnum << " " << endl;
		}
	#else
		Funcdist f = (f1(x)<f2(x) ? f1 : f2);

		try{//Integrate distance function from 0-2M_PI
			area=r0*qsimp(f, 0.0, 2.0*M_PI, 1.0e-8);
		}
		catch(const char * s){
			cout << s << ": " << (*vi).pnum << " " << (*vj).pnum << endl;
		}
		//Map class 
		mapplanes[area]=100*(*vi).pnum+(*vj).pnum;	//map sorts by area
	#endif
#endif
		vj++;
		if(vj==vecplanes.end()){
			++vi;
			vj=vi+1;
		}
		if (vj==vecplanes.end())
			break;	
	}
	cout.setf(ios_base::showpoint);		//allow trailing zeros in floats
	for (map<double,int>::iterator mj=mapplanes.begin();
				mj!=mapplanes.end();mj++){
		cout << (*mj).first << " ";
		cout.width(4);		//Affects only next item displayed
		cout.fill('0');
		cout << (*mj).second << endl;
		cout.fill(' ');		//Revert to default (space)
	}

#ifdef CHECK_ORTHOGONAL
		forthog.close();
#endif
	return 0;
}

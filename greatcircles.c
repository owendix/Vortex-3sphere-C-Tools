#include "linalg.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>

//#define OUTPUTSPLIT	//define if you want to plot the hopf projection
	//with several files: or put 2 spaces between them then gnuplot
	//will separate
					
using namespace std;
/*************************************************
Calculate trajectory due to hopf vector field: v=(-y,x,-w,z)
If project=1, projects that R^4 trajectory to either hopf or hyperspherical 
coords: projectmethod=0,1 respectively
Can choose to have input from file (selected during runtime), 
or randomly chosen
**************************************************
Input: g++ -g greatcircles.c linalg.c -o gc.out
		./gc.out
Output: #ifndef OUTPUTSPLIT
		(To: gcircles.r4d.dat, gcircles.hopf.dat, or gcircles.hsph.dat)
		out0 out1 out2	[x(steps+1)]
		out0 out1 out2 	[x(steps+1)]
		... (n times)

		#ifdef OUTPUTSPLIT
		(xxr4d1 or xxhopf1 or xxhsph1):
		out0 out1 out2	[x(steps+1)]
		...
		(xxr4d2 or xxhopf2 or xxhsph2):
		out0 out1 out2 	[x(steps+1)]
		...

		repeated n times

	-> splits into xx1 xx2 ... xxn
**************************************************/
double r0;

point randS3point();
void hopffib(point, double*);
void hspher(point, double*);
void stereo(point, double*);
bool hitwall(double*, double *);

int main(){
	double d, r[4], dt, proj[3], prevproj[3];
	int n=15, i, j, k, s, steps=80, fac=1, tmpint;
	char buf[22];
	double tmpdb;
	point s3point, v;
	ofstream fout, fplot, finproj;
	ifstream fin;
	string infname;
	
	int project;	//=0, don't project, =1:  hopffib; =2, hspher,=3 stereo
	cout << "Set project = 0 (R4), 1 (hopf), 2 (hyperspher), 3 (stereo) \n";
	cout << "project: ";
	cin >> project;
	if (!cin.good()){	//returns 0 if problem encountered
		perror("Invalid parameter format\n");
		exit(1);
	}else if (project!=0 && project!=1 && project!=2 && project!=3){
		perror("Invalid parameter value\n");
		exit(1);
	}
	int infromfile;	//=1, use restart file from 3sphere output for start
						//point, =0, use random points (set n)
	cout << "Choose random start points (input: 0), or input from restart ";
	cout << "file (input: 1)\n";
	cout << "Input source (0,1): ";
	cin >> infromfile;
	if (!cin.good()){
		perror("Invalid parameter format\n");
		exit(1);
	}else if (infromfile!=0 && infromfile!=1){
		perror("Invalid parameter value\n");
		exit(1);
	}

#ifdef OUTPUTSPLIT
	string fbase="xx", fname, strnum;
	if (!project)
		fbase+="r4d";
	else if (project==1)
		fbase+="hopf";
	else if (project==2)
		fbase+="hsph";
	else
		fbase+="ster"
	stringstream strhandle;
#else
	if (!project)
		sprintf(buf,"gcircles.r4d.dat");
	else if (project==1)
		sprintf(buf,"gcircles.hopf.dat");
	else if (project==2)
		sprintf(buf,"gcircles.hsph.dat");
	else
		sprintf(buf,"gcircles.ster.dat");
	fout.open(buf);
#endif
	
	if (!infromfile){
		srand(time(0));
		r0=1;	//randomly generated point on S^3, r0=1
	}else if (infromfile==1){//will only be 0 or 1
		//read in from a restart file
		cout << "Which restart file?: ";
		cin >> infname;
		if (!cin.good()){
			perror("Invalid in-filename assignment\n");
			exit(1);
		}
		fin.open(infname.c_str(),ios_base::in);
		if (!fin.good()){
			perror("Failed attempt to open file\n");
			exit(1);
		}
		//deal with header information
		//set n
		fin >> n;	//totalpts
		fin >> tmpint;	//nvort
		fin >> tmpdb;	//time
		fin >> tmpdb;	//dt
		fin >> tmpint; fin >> tmpint; fin >> tmpint;	//start,end,vtype
		finproj.open("infile.projected.dat");
	}
	//srand(1);
	dt=2*M_PI/(steps-1);	//make it just close if its a great circle

#ifdef OUTPUTSPLIT
	cout << "OUTPUTSPLIT on: output to " << fbase << "1-";
	cout << fbase << n << endl;
#else
	cout << "OUTPUTSPLIT off: output to " << buf << endl;
#endif
	if (project){//hopffib or hspher
		if (project==1){//hopffib
			cout << "Run: $gnuplot plot.hopf.gcircles\n";
			fplot.open("plot.hopf.gcircles");
		}else if (project==2){
			cout << "Run: $gnuplot plot.hsph.gcircles\n";
			fplot.open("plot.hsph.gcircles");
		}else{
			cout << "Run: $gnuplot plot.ster.gcircles\n";
			fplot.open("plot.ster.gcircles");
		}
		if (project==1)//hopffib
			fplot << "set size 0.7,1.0\n";
		else if (project==2)
			fplot << "set size 0.5,1.0\n";
		fplot << "set xlabel \'x\'\n";
		fplot << "set ylabel \'y\'\n";
		fplot << "set zlabel \'z\'\n";
		fplot << "set ticslevel 0\n";
		fplot << "set pointsize 0.7 \n";
		fplot << "set title \"n=" << n << "\"\n";
		if (project==3){//stereographic, limit size
			tmpdb=2.5;
			fplot << "set xrange [-" << tmpdb*r0 << ":";
			fplot << tmpdb*r0 << "]\n";
			fplot << "set yrange [-" << tmpdb*r0 << ":";
			fplot << tmpdb*r0 << "]\n";
			fplot << "set zrange [-" << tmpdb*r0 << ":";
			fplot << tmpdb*r0 << "]\n";
		}
#ifndef OUTPUTSPLIT//gnuplot's cool feature
		fplot << "splot for [i=1:" << n << "] \"" << buf;
		fplot << "\" index i w linespoints pt 7 notitle" << endl;
#endif
		/*pt:0=dot,1= +,2= x,3= *,4=hollow square w/ dot, 
          5=solid square,6=hollow circle w/ dot,7=solid circle,
          8=hollow triangle w/ dot, 9=solid triangle,
          10=hollow inverted triangle w/dot, 11=solid triangle
          12=hollow rhombus w/dot, 13=solid rhombus*/
	}
	//for all initial points on 3sphere
	for (i=1;i<=n;i++){	//if infromfile==1, n=totalpts
		if (!infromfile){
			//Generate random value for point on S^3 (radius 1)
			/*Requires seeding srand(time(0))
			Requires MACROS defined: RAND_INT(MIN,MAX), RAND_DB*/
			s3point=randS3point();
		}else{//uses restart file as initial points
			//read in last four entries from line, the point structure
			//for each point on vortex
			fin >> tmpint; fin >> tmpint;	//index,recon#
			fin >> s3point;
			if (i==1){
				r0=sqrt((s3point.x*s3point.x+s3point.y*s3point.y)
					+(s3point.z*s3point.z+s3point.w*s3point.w));
				if (project==1){
					fplot << "set xrange [-" << r0*r0 << ":";
					fplot << r0*r0 << "]\n";
					fplot << "set yrange [-" << r0*r0 << ":";
					fplot << r0*r0 << "]\n";
					fplot << "set zrange [-" << r0*r0 << ":";
					fplot << r0*r0 << "]\n";
				}else if (project==3){
					tmpdb=1.5;
					fplot << "set xrange [-" << tmpdb*r0 << ":";
					fplot << tmpdb*r0 << "]\n";
					fplot << "set yrange [-" << tmpdb*r0 << ":";
					fplot << tmpdb*r0 << "]\n";
					fplot << "set zrange [-" << tmpdb*r0 << ":";
					fplot << tmpdb*r0 << "]\n";
				}
			}
		}
		fout.precision(17);
	
		//project s3point
		if (project){//hopf or hspher
			if (project==1)	//hopf
				hopffib(s3point, proj);
			else if (project==2)	//hspher
				hspher(s3point, proj);
			else //stereo
				stereo(s3point, proj);
		}
		//plot file and outfile stuff (if split)
#ifdef OUTPUTSPLIT
		strhandle << i;
		strnum = strhandle.str();
		fname=fbase+strnum;
		strhandle.str("");
		fout.open(fname.c_str());
		if (project){//hopffib, hspher,stereo: can't print unless projected
			if (i==1){
				fplot << "splot";
			}
			fplot << " \"" << fname.c_str() << "\" w points";
			if (i!=1)
				fplot << " notitle";
			else
				fplot << " t \'xx1-xx" << n << "\'";
			if (i!=n)
				fplot << ",";
			else if (infromfile)
				fplot << ",";
			else
				fplot << endl;
		}
		//write projection to original infile
		if (infromfile==1){
			//write original infile
			if (project){
				//hopffib or hspher: can't plot unless projected
				if (i==n){
					fplot << " 'infile.projected.dat'";
					fplot << " w linespoints lw 3 lc 1 t '";
					fplot << infname.c_str() << "'" << endl;
				}
				for (k=0;k<3;k++){
					finproj << proj[k];
					if (k!=2)
						finproj << " ";
					else
						finproj << endl;
				}
			}
		}
#endif
		//plot projected
		if (project){
			for (k=0;k<3;k++){
				if (project==2){//hsphere: to close
					prevproj[k]=proj[k];
				}
				fout << proj[k];
				if (k!=2)
					fout << " ";
				else
					fout << endl;
			}
		}else{
			fout << s3point << endl;
		}
		//print all the points iterated from 1st initial point
		for (j=0;j<fac*steps;j++){
			//calculate velocity, v
			v.x=-s3point.y; v.y=s3point.x; 
			v.z=-s3point.w; v.w=s3point.z;
			/*Integration method: assumes vector to new point is v*dt.
			Scales initial s3point, such that the perpendicular 
			displacement vector v*dt can be added to yield the new 
			point on s3. This method only works because the path will 
			not be leaving a great circle.*/ 

			//Works even if s3point is not unit radius, because
			//v is based on s3point, so would also not be of unit radius
			s3point=s3point*sqrt((1-dt)*(1+dt))+v*dt;
			if (project){//hopf or hspher or stereo
				if (project==1){	//hopf
					hopffib(s3point, proj);
				}else if (project==2){	//hspher, make so I can use linespoints
					hspher(s3point, proj);
					if (hitwall(proj,prevproj)){
						fout << "\n";	//hit wall, add ONE blank line
					}
				}else{//stereographic projection
					stereo(s3point, proj);
				}
				for (k=0;k<3;k++){
					prevproj[k]=proj[k];//reassign prevproj[k]
					fout << proj[k];
					if (k!=2)
						fout << " ";
					else
						fout << endl;
				}
			}else{
				fout << s3point << endl;
			}
			
		}//end for all points from one initial point
		
#ifdef OUTPUTSPLIT
		fout.close();
#else//want 2 spaces between vortices so gnuplot will print different colors
		if (i!=n)//not for last vortex
			fout << "\n\n";
#endif
	}//end for all initial points
#ifndef OUTPUTSPLIT
	fout.close();
	if (infromfile)
		finproj.close();
#endif
	if (infromfile){
		fin.close();
	}
	if (project){
		fplot.close();
	}
	return 0;
}

/*Creates a random point on S^3, (R^4 with radius 1)
Requires seeding srand((unsigned)time(0))
Requires MACROS defined: RAND_INT(MIN,MAX), RAND_DB*/
point randS3point(){
	const static int max=1, min=0;
	static int s, j;
	static double tmpdb1,tmpdb2, r[4];
	static point q;
	
	//Generate random value for point on S^3 (radius 1)
	do {
		//random integer between min and max
		s=rand()%(max-min+1)+min;	//s = 0, or 1 (max=1,min=0)
		s=2*s-1;	//s = +1 or -1
		tmpdb2=0.0;
		for (j=0;j<3;j++){
			//can't seed it right before: produces bad numbers 
			tmpdb1=((double)rand()+(double)rand()/(double)RAND_MAX)
				/((double)RAND_MAX); //rand double from 0->1 inclusive
			//random nums 0<=d<=1: added factor is to get finer spacing
    		//otherwise won't produced numbers between 0 and 1/RAND_MAX
		 	r[j]=2.0*tmpdb1-1.0;	//stretch random doubles from -1 -> +1
			tmpdb2+=r[j]*r[j];
		}
	} while (tmpdb2>1.0);
	//complete valid point on 3sphere
	r[3]=s*sqrt(1-tmpdb2);
	//make it a point
	q.x=r[0]; q.y=r[1]; q.z=r[2]; q.w=r[3];
	
	return q;	
	
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
	h[0]=M_PI-atan2(sqrt(r.w*r.w+r.z*r.z+r.y*r.y),r.x);
	h[1]=M_PI-atan2(sqrt(r.w*r.w+r.z*r.z),r.y);
	h[2]=M_PI-atan2(r.w,r.z);

	return;
}

void stereo(point r, double* h){
	//stereographic coordinates
	static double denom;
	denom=1+r.w/r0; //assumes projection point is -w
	//denom=1-r.w/r0; //assumes projection point is +w
	if (denom<1e-12)
		denom=1e-12;
	else if (denom)
	h[0]=r.x/denom;
	h[1]=r.y/denom;
	h[2]=r.z/denom;
	
	return;
}

//takes hyperspherical coordinates (doubles array of length 3)
//returns 1 if the points intersect a wall, 0 if they don't
bool hitwall(double* p1, double * p2){
	//if distance between walls is more than half the distance between 
	//limits of that variable, then it hit the wall
	return (abs(p1[0]-p2[0])>M_PI/2.0 || abs(p1[1]-p2[1])>M_PI/2.0 || 
		abs(p1[2]-p2[2])>M_PI);
	//3rd hyperspherical angle has different limits
}

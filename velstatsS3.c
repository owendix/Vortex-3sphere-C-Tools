#include "ubiq.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <cerrno>
using namespace std;
/* 3-sphere analysis of velocity stats. Operates on restart files. 
Outputs to stdout
*************************************************************** 
If  params files were pre aug0412, need to include
absmax_fac: absmax2=absmax_fac*absmax_fac*r0*r0;
absmin_fac: absmin2=absmin_fac*absmin_fac*emax*emax*r0*r0;
limbiot_fac: limbiot=limbiot_fac*r0*r0;	//different for 3sphere than torus
After aug0412, these are part of the params files
Compile:
	Requires velstats.c linalg.c ubiq.h, linalg.h
	g++ -g -O3 velstats.c linalg.c -o vels.out

Input: 	
	Pre-aug0412 trials:
	./vels.out ../../code/pfiles/parjul3112a.dat absmax_fac absmin_fac limbiot_fac ../../code/jul3112a nskip

	//For reference:
	//absmax2=absmax_fac*absmax_fac*r0*r0;
	//absmin2=absmin_fac*absmin_fac*emax*emax*r0*r0;
	//limbiot=limbiot_fac*r0*r0;//different from torus for accuracy

	Post-aug0412 trials:
	./vels.out ../../code/pfiles/paraug0512a.dat ../../code/aug0512a nskip

Output:
	//cout: 
	//standard deviation of preceding quantity is printed on even columns
    Column    Quantity (Averages of)
        [1]     [time]
        2       dot(tangent,u1)
        4       radius      
        6       beta
        8       |vlocal|
		10		vlocal-u1
		12		vlocal-u2
		14		vlocal-u3
        16      |v_nonlocal|
        18      v_nl-u1
        20      v_nl-u2
        22      v_nl-u3
        24      |v_fric|
        26      v_fric-u1
        28      v_fric-u2
        30      v_fric-u3
        32      |v_total|
        34      v_tot-u1
        36      v_tot-u2
        38      v_tot-u3
		40		Fsn.u1
		42		Fsn.u2
		44		Fsn.u3
		46		spacing (no std dev)
***********************************************************************
\rho(total)=145kg/m^3 (0<T<T\lambda)
\rho=\rho_n+\rho_s;
Barenghi,Donnelly,Vinen. J. Low Temp.Phys., 52,189 (1983):
@T=1.6K:\alpha=0.098,\alpha'=0.012,\rho=0.1452g/cm^3, \rho_s=0.1213 g/cm^3
***********************************************************************
Reads if NONLOCAL was on from params file, and sets global vble accordingly 
Assumes VSHOPF was off from params file
Assumes SPACCUTOFF was on: variable pt spacing and reconnect dist
*************************************************************************/

//IO function prototypes**************************************************
void readparams(string, int);
void readcores(point[], ptadmin[],ifstream&);
//rkf stuff
void nonlocalvel(point*, point*, double);
void derivs(point[],point[]);
point biot(point, point, point, point);
point localrz(point, point, point, point&, double&, double&);
//************************************************************************

/* Global variables for the module, entered in inout.c from params.dat */
int abscomps=0;//=1, use abs(components) in derivs, or =0
vortex startcore[NVORTMAX];
int outwrite;
double r0, absmax2, absmin2;

double tinit;	//really t: for a particular file
double vn, vs, emax, skew;
double a0, alpha;
double absmin_fac, absmax_fac, limbiot_fac;
int ilocal, ifric, inotan, ncycles;
int nonlocalon;

/* Other global variables */
double const1 = KAPPA / 4. / M_PI, const2;
double rhos=0.1213;//g/cm^3: if \alpha\approx0.1 (T=1.6K)
/* Global variables */
int nvort, totalpts;
double radius[NMAX];
point core[NMAX];
ptadmin periph[NMAX];

#ifdef VELSTRUCT
char velstrfile[39];
#endif

int main(int argc, char* argv[]) 
{
	/*argv order
	IF params file is later than aug 04 2012 (argc=3)
	./velstats ../../code/pfiles/parfilename.dat ../../code/filename nskip
	IF params file is earlier than aug 04 2012 (argc=5)
	./velstats ../../code/pfiles/parfilename.dat absmax_fac absmin_fac limbiot_fac ../../code/filename nskip */
	string infilepre, infilemid, infilebase, infile, paramsfile, strnum;
	stringstream strhandle;
	ifstream infhandle;
	double limbiot;
	int loc, nskip, done, num, numlow;
	point  nonloc[NMAX];

	if (argc==4 || argc==7){	//argv[argc] is null pointer
		nskip=atoi(argv[argc-1]);
		infilepre=(string)argv[argc-2];//no . or final /
		paramsfile=(string)argv[1];	//full params file path in pfiles/
	}else if (argc==5 || argc==8){
		nskip=atoi(argv[argc-1]);
		abscomps=atoi(argv[argc-2]);//=1, use abs(components) in derivs, or =0
		infilepre=(string)argv[argc-3];//no . or final /
		paramsfile=(string)argv[1];	//full params file path in pfiles/
	}else{
		cerr << "If params file is pre-aug0412, Input:\n";
		cerr << "  paramsfile absmax_fac absmin_fac limbiot_fac restartfilepre";
		cerr << " nskip\n";
		cerr << "If params file is post-aug0412, Input:\n";
		cerr << "  paramsfile restartfilepre nskip\n";
		cerr << "argc=" << argc << endl;
		cerr << "argv[argc-1]=" << argv[argc-1] << endl;
		exit(1);
	}
	if (nskip<1)
		nskip=1;
	if (abscomps!=1)
		abscomps=0;
	cerr << "infilepre = " << infilepre.c_str() << endl;
	cerr << "paramsfile = " << paramsfile.c_str() << endl;
	cerr << "nskip = " << nskip << endl;
	cerr << "abscomps = " << abscomps << endl;

	//Concatenate filename to path
    if ((loc=infilepre.find_last_of("/"))==string::npos)
        perror("Error in find_last_of()\n");
    infilemid=infilepre;
	//loc is position (last / before filename), 9 is length of "/mon3112a"
	//assign replaces old content with this "/mon3112a"
    infilemid.assign(infilepre,loc+1,8);//filename with leading /
    infilebase=(infilepre+"/")+infilemid; 
		//concat with num to make infile, inside while loop

	/*I want readparams to read in r0, NONLOCAL, absmax_fac, absmin_fac,
	limbiot_fac. Some I may have to input myself with argv. If 
	argc certain size, override
	*/
	//Read params from a specific file, included in argv
	//This params file is common to all in the list of filenames (restarts)
	readparams(paramsfile, argc);
	const2=exp(0.25)*a0;//for derivs()
	
	//read absmax_fac, absmin_fac, limbiot_fac from argv, or from params 
	//if params file is later than aug 04 2012
	if (argc==7 || argc==8){//strtod more robust than atof (sets errno)
		//sets to 0.0 upon failure
		absmax_fac=strtod(argv[2],NULL);
		absmin_fac=strtod(argv[3],NULL);
		limbiot_fac=strtod(argv[4],NULL);
		if (absmax_fac==0 || absmin_fac==0 || limbiot_fac<0 || errno){
			cerr << "Problem assigning vals from argv[], ";
			cerr << strerror(errno) << endl;
			exit(1);
		}
		cerr << "absmax=" << absmax_fac << "*r0" << endl;
		cerr << "absmin=" << absmin_fac << "*emax*r0" << endl;
		cerr << "limbiot=" << limbiot_fac << "*r0*r0" << endl;
		cerr << "r0=" << r0 << ", vn=" << vn << ", vs=";
		cerr << vs << ", emax=" << emax << endl;
		if (nonlocalon)
			cerr << "NONLOCAL on" << endl;
		else
			cerr << "NONLOCAL off" << endl;
	}
	absmax2=absmax_fac*absmax_fac*r0*r0;
	absmin2=absmin_fac*absmin_fac*emax*emax*r0*r0;
	limbiot=limbiot_fac*r0*r0;//different from torus for accuracy

	#ifdef VELSTRUCT
    sprintf(velstrfile,"%s.%s.%s","velstruct",infilemid.c_str(),"dat");
    #endif
	
	//while opening the filename file was successful
	num=15839;
	numlow=num;
	cerr << "Column    Quantity (Averages of)\n";
	cerr << "[1]    [time]\n";
	cerr << "2      dot(tangent,u1)\n";
	cerr << "4      radius\n";
	cerr << "6      beta\n";
	cerr << "8      |vlocal|\n";
	cerr << "10		vlocal-u1\n";
	cerr << "12		vlocal-u2\n";
	cerr << "14		vlocal-u3\n";
	cerr << "16     |v_nonlocal|\n";
	cerr << "18     v_nl-u1\n";
	cerr << "20     v_nl-u2\n";
	cerr << "22     v_nl-u3\n";
	cerr << "24     |v_fric|\n";
	cerr << "26     v_fric-u1\n";
	cerr << "28     v_fric-u2\n";
	cerr << "30     v_fric-u3\n";
	cerr << "32     |v_total|\n";
	cerr << "34     v_tot-u1\n";
	cerr << "36     v_tot-u2\n";
	cerr << "38     v_tot-u3\n";
	cerr << "40		Fsn.u1\n";
	cerr << "42		Fsn.u2\n";
	cerr << "44		Fsn.u3\n";
	cerr << "46		spacing (no std dev)\n";
	while (1){
		strhandle << num;
		strnum=strhandle.str();
		infile=infilebase+"."+strnum;
		strhandle.str("");//clears it for next use
		infhandle.open(infile.c_str());
		if (!infhandle.good()){//would exit upon eofbit being set
			//Not necessarily an error, but I'll be > stdout to file
			cerr << "Unable to open file " << infile.c_str() << endl;
			exit(0);
		}
	
		//Read file from "restart" file: one of a list included in argv
		//may not be able to pass ifstream handle
		readcores(core, periph, infhandle);
	
		//Calculate the nonlocal velocity field
		nonlocalvel(core, nonloc, limbiot);
		//Call derivs which prints either velstats (stdout) 
		//or velstats (stdout) AND velstruct (velstruct.filename.dat)
		derivs(core,nonloc);

		infhandle.close();
		num+=nskip;
	}	
	
	return (0);
}

void nonlocalvel(point* y, point* nonloc, const double limbiot){
	static int i, j, ivort;
	static double dist2, skiplength, jnkdb;
	//double limbiot=2*r0*r0; Set in argv
	static point errsum, tmppt, ypt, basisv1;
	static const double maxlfactor=0.2, minlfactor=0.0833333333333;
	static double arclength, arclimbiot;
	//the values are taken from reconnect.c
	//y is the core value
	
	if (limbiot==2.0*r0*r0)
        arclimbiot=r0*M_PI/2.0; //quarter of a great circle
    else
        arclimbiot=r0*acos(1.0-limbiot/(2.0*r0*r0));
	for (i=0; i<totalpts; i++){//i is point calculating nonlocal ON, by pts j
		nonloc[i].x = 0; nonloc[i].y = 0; nonloc[i].z = 0; nonloc[i].w=0;
		//Error for compensated summation: Higham (2002) pg 84
		errsum.x=0.0; errsum.y=0.0; errsum.z=0.0; errsum.w=0.0;
		basisv1 = y[i]/sqrt(dot(y[i],y[i]));
		for (ivort=0; ivort<nvort; ivort++){
			j=startcore[ivort].start;
			do{
				if (j!=i && periph[j].ip !=i){//j is source of nonlocal on i
					if ((dist2=dot(y[j]-y[i],y[j]-y[i]))<limbiot){
						tmppt = nonloc[i];
						ypt = errsum + biot(y[j],y[periph[j].ip],
								y[i], basisv1);
						nonloc[i] = tmppt + ypt;
						errsum = (tmppt - nonloc[i]) + ypt;
						//w/ compensated sum method
					}else{
						arclength=r0*acos(1.0-dist2/(2.0*r0*r0));
						skiplength=arclength-arclimbiot;
						while (j!=startcore[ivort].end && skiplength > 0){
							j=periph[j].ip;
							jnkdb=maxlfactor*radius[j];
							if (jnkdb*jnkdb>absmax2)
								jnkdb=sqrt(absmax2);
							else if (jnkdb*jnkdb<absmin2)
								jnkdb=sqrt(absmin2);
							skiplength -= jnkdb;
						}
					}
				}
				j=periph[j].ip;
			}while (j!=startcore[ivort].start);
		}//end for (ivort)
		nonloc[i] = nonloc[i] + errsum;	//Final correction: comp sum
	}
}

point biot(point stseg, point endseg, point testpt, point v1)
{
    static point l1, l2, l3, deldot, v[4], tmppt, dir;
    static double dot11, dot12, dot22, dot33, dot23, factor1, factor2;
    static double tmpd[4], tmpv[4], angle, c, s, Rij, psum, va[4][4];
    static int i, j;

    /*testpt is the point we are calculating the velocity field at, 
    due to the vortex segment l2=endseg-stseg.
    *deldot is the integrand for this one vortex segment, solved 
    as a piece of the integral assuming a straight vortex segment, l2.
    we perform this integral, with the cross product, at the segment 
    l2 location, either at stseg or endseg. then parallel transport the 
    result to testpt, so we can add all the integral pieces together*/

    l1 = testpt - stseg;
    l2 = endseg - stseg;
    l3 = testpt - endseg;
    dot11=r0*acos(1.0-dot(l1,l1)/(2.0*r0*r0));  //arclength
    dot33=r0*acos(1.0-dot(l3,l3)/(2.0*r0*r0));
    l1 = l1-dot(l1,stseg)*stseg/(r0*r0);    //Find tangent vec at stseg
    l1 = l1*dot11/sqrt(dot(l1,l1));     //Normalize to arclength
    l3 = l3-dot(l3,endseg)*endseg/(r0*r0);  //Find tangent vec at endseg
    l3 = l3*dot33/sqrt(dot(l3,l3));     //Normalize to arclength

    //Neglects curvature effects over distance mag(l2)  

    dot11 = dot11*dot11;    //need arclength squared
    dot33 = dot33*dot33;    //need arclength squared
    dot12 = dot(l1, l2);    //angle of vecs at approx stseg
    dot22 = dot(l2, l2);
    dot23 = dot12 - dot22;

    factor1 = dot11*dot22-dot12*dot12;
    //deldot calculated at particular point on 3sphere
    if (fabs(factor1/dot11/dot22) > 1.0e-14){
        factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
        deldot = (factor2 / factor1) * cross(l1, l2, stseg)/r0;
        angle=-sqrt(dot11)/r0;  //Angle between stseg and testpt
        v[1]=stseg - dot(v1,stseg)*v1;
    }else{
        factor1 = dot22*dot33-dot23*dot23;
        if  (fabs(factor1/dot22/dot33) > 1.0e-14){
            factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
            deldot = (factor2 / factor1) * cross(l3, l2, endseg)/r0;
			 angle=-sqrt(dot33)/r0;  //Angle between endseg and testpt
            v[1]=endseg - dot(v1,endseg)*v1;
        }else{
            if (dot11*dot33 > 1.0e-12 && dot12*dot23 > 0.0) {
                factor2 = fabs(dot11-dot33)/dot11/dot33/sqrt(dot22);
                deldot = (factor2 / 2.0) * cross(l1, l2, stseg)/r0;
                angle=-sqrt(dot11)/r0;  //Angle between stseg and testpt
                v[1]=stseg - dot(v1,stseg)*v1;
            }else if (dot12*dot23 <= 0.0){
                factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
                dir = cross(l1, l2, stseg)/r0;
                factor1 = dot(dir, dir);
                deldot = (factor2 / factor1) * dir;
                return deldot;  //ignore angle diff over dist l2
            }else{
                deldot.x=0.0; deldot.y=0.0; deldot.z=0.0; deldot.w=0.0;
                return deldot;
            }
        }
    }
    //The rotation matrix rotates deldot back to the testpt to perform 
    //a wedge product between it and some other vectors located there.
    //Now, need to rotate deldot back to testpt
    //Establish basis vectors of plane of rotation
    v[0]=v1;
    v[1] = v[1]/sqrt(dot(v[1],v[1]));
    //Go through euclid vecs to find 3rd rotation basis vec
    v[2].x=1.0; v[2].y=0.0; v[2].z=0.0; v[2].w=0.0;
    for (i=0;i<3;i++){
        tmppt = v[2]-dot(v[0],v[2])*v[0]-dot(v[1],v[2])*v[1];
        tmppt=tmppt/sqrt(dot(tmppt,tmppt));
        if(dot(tmppt,tmppt)<1.0e-8){
            //If orthogonal, go on to next euclid vec
            if (i==0){
                v[2].x=0.0;
                v[2].y=1.0;
            }else if(i==1){
                v[2].y=0.0;
                v[2].z=1.0;
            }
        }else{
            v[2]=tmppt;
            break;
		}
    }
    //4th rotation vec
    v[3]=cross(v[0],v[1],v[2]);
    //To iterate Rij, need to make point structs into arrays:
    //v[i].x,y,z,w = va[i][0,1,2,3]
    for (i=0;i<4;i++){
        va[i][0]=v[i].x;    //Notation va = v_array
        va[i][1]=v[i].y;
        va[i][2]=v[i].z;
        va[i][3]=v[i].w;
    }

    tmpv[0]=deldot.x;   //Makes matrix mult easier
    tmpv[1]=deldot.y;
    tmpv[2]=deldot.z;
    tmpv[3]=deldot.w;

    c=cos(angle);
    s=sin(angle);
    //Mult rotation matrix w pairwise summation over j
    for (i=0;i<4;i++){
        psum=0.0;
        for(j=0;j<4;j++){
            //Confirmed, accurate to 1e-15, 02-09-11: cpptools/rotate4d.cpp
            Rij = ((c*va[0][i]+s*va[1][i])*va[0][j] + va[2][i]*va[2][j])
                + ((-s*va[0][i]+c*va[1][i])*va[1][j] + va[3][i]*va[3][j]);
            //pairwise summation
            psum += Rij*tmpv[j];
            if (j==1){
                tmpd[i] = psum;
                psum=0.0;
            }
        }
        tmpd[i] +=psum;
    }
    deldot.x = tmpd[0];
    deldot.y = tmpd[1];
    deldot.z = tmpd[2];
    deldot.w = tmpd[3];

    return deldot;
}

/* Evaluates the derivatives (RHS) of the system of ODE's x' = f(t,x) */

void derivs(point corein[], point nonloc[]) {
    //corein contains input values, coretemp is dummy array: they are often
    //BUT NOT ALWAYS the same -> saves space and time: no need to 
    //create a temporary variable (like coretemp[]) in derivs
    //to store scaled data, no need to copy
    //ybfr5[i]=y[i], before first call to derivs, and no need to
    //scale back coretemp at end, because its a dummy array
	static point f[NMAX], coretemp[NMAX];
    static int i, j;
    static point sp; /* unit vector approximating tangent to vortex line */
    static point tmp;
	static double scale[NMAX];
    static double beta, space, spaceprev, weight, sumweight, tempwt, Rwt;
	static point fD;
    static const double cunity=1.1;
    //const2 = exp(.25)*a0;//now initialized once in main()
    static point u1, u2, u3, temppt, norm, bnorm;
    //Print, for each vortex, the avg and std dev in several
    //terms of the velocity, as component WRT orthonormal velocity vectors
    static int k, l, ivort;
    static double num[22], dif;
	static struct statsdb{
		double m;   //Like avg
		double q;   //Like variance
	} stats[22];
	//initialize stats
	for (k=0;k<22;k++){
		stats[k].m=0.0;
		stats[k].q=0.0;
	}
	sumweight=0.0;//line length
#ifdef VELSTRUCT
    //prints, for each point on a vortex, components of the 
    //tangent vector (WRT orthonormal velocity vectors), and 
    //different velocity terms (local, nonlocal, friction, total).
    // No averaging, lots of printing: expensive
    ofstream fvstruct(velstrfile, ios_base::out | ios_base::app);
#endif
    // Project all core points onto the three-sphere of radius r0.
	for (ivort=0;ivort<nvort;ivort++){
		i=0;
		j=startcore[ivort].start;
		do{
			scale[i] = r0 / sqrt(dot(corein[j], corein[j]));
			coretemp[j] = corein[j]*scale[i];
			i++;
			j=periph[j].ip;
		}while (j!=startcore[ivort].start);
		i=0;
		//j=startcore[ivort].start;		//should be true anyway
    	do{
			//Initialize the derivative vector 
			f[i].x = 0.;
			f[i].y = 0.;
			f[i].z = 0.;
			f[i].w = 0.;

			/* Find the local self-induced field */
			if (ilocal+ifric+inotan != 0){
				f[i] = localrz(coretemp[periph[j].im], coretemp[j],
						coretemp[periph[j].ip], sp, radius[j], beta);
				if (ilocal !=0){
					if (!nonlocalon){
						//beta = const1*log(r0/a0);
						beta = cunity*const1*log(radius[j]/a0);
					}
					f[i] = beta*f[i];
				}
			}
			//establish basis vecs: velocity dir
			u1.x=-coretemp[j].y/r0; u1.y=coretemp[j].x/r0;
			u1.z=-coretemp[j].w/r0; u1.w=coretemp[j].z/r0;
			u2.x=u1.w; u2.y=u1.z; u2.z=-u1.y; u2.w=-u1.x;
			u3.x=u1.z; u3.y=-u1.w; u3.z=-u1.x; u3.w=u1.y;
			
			temppt=coretemp[periph[j].ip]-coretemp[j];
			space=sqrt(dot(temppt,temppt));//spacing for weighted avg
			//establish weighting factor for averages/stdevs
			if (i==0){
                temppt=coretemp[j]-coretemp[periph[j].im];
                spaceprev=sqrt(dot(temppt,temppt));
            }
            weight=0.5*(space+spaceprev);
			num[0]=dot(sp,u1);//alignment of tangent vector: also from Ipar
			num[1]=radius[j];
			num[2]=beta;
            num[3]=sqrt(dot(f[i],f[i]));//local velocity
			num[4]=dot(f[i],u1);//vlocal-u1
			num[5]=dot(f[i],u2);//vlocal-u2
			num[6]=dot(f[i],u3);//vlocal-u3
			if (nonlocalon){
				//Add nonlocal field
				//Remove radial component
				tmp=nonloc[j]-(dot(nonloc[j],coretemp[j])/(r0*r0))*coretemp[j];
				f[i] = f[i] + const1*tmp;
				num[7]=const1*sqrt(dot(tmp,tmp));	//v_nl
				num[8]=const1*dot(tmp,u1);// v_nl||u1
				num[9]=const1*dot(tmp,u2);// v_nl||u2
				num[10]=const1*dot(tmp,u3);// v_nl||u3
			}else{
				num[7]=0.0;
				num[8]=0.0;
				num[9]=0.0;
				num[10]=0.0;
			}
			
			tempwt=weight+sumweight;//do this once, now 
            for (k=0;k<11;k++){
				if (k>3 && k!=7 && abscomps)//only components
                	dif = abs(num[k]) - stats[k].m;
				else
                	dif = num[k] - stats[k].m;
                Rwt=dif*weight/tempwt;
                stats[k].m += Rwt;
                stats[k].q += sumweight*dif*Rwt;
            }
            //sumweight=tempwt;do this at the end
			/* Add imposed superfluid flow velocity */
		/*#ifdef VSHOPF: code assumes vshopf is off
			tmp.x = -coretemp[j].y;
			tmp.y = coretemp[j].x;
			tmp.z = -coretemp[j].w;
			tmp.w = coretemp[j].z;
			f[i] = f[i] + vs * tmp / r0;
		#else   //VSHOPF*/
			tmp.x = 0; tmp.y = 0; tmp.z = 0; tmp.w = 0;
			f[i] = f[i] + vs * tmp;
		//#endif  //VSHOPF

			// Add normal component for friction term
			temppt=f[i];	//Save prior to adding friction
		/*#ifdef VSHOPF
			tmp.x = 0;
			tmp.y = 0;
			tmp.z = 0;
			tmp.w = 0;
			// Find the friction-induced field
			if (ifric != 0)
				f[i] = f[i] + (alpha / r0) * cross(sp, (vn * tmp - f[i]),
						coretemp[j]);
		#else   //VSHOPF*/
			tmp.x = -coretemp[j].y; tmp.y = coretemp[j].x;
			tmp.z = -coretemp[j].w; tmp.w = coretemp[j].z;

			//tmp has magnitude r0, divide it out: want |vn*(dir)| = vn
			//Need |vn|*tmp2/r0
			/* Find the friction-induced field */
			if (ifric != 0)
			f[i] = f[i] + (alpha/r0)*cross(sp,(vn*tmp/r0-f[i]),coretemp[j]);
		//#endif  //VSHOPF
			temppt=f[i]-temppt;//get friction term;
			//get drag force on each vortex segment
			fD=-rhos*KAPPA*cross(sp,temppt,coretemp[j])/r0;
            num[11]=sqrt(dot(temppt,temppt));//friction magnitude
            num[12]=dot(temppt,u1);
            num[13]=dot(temppt,u2);
            num[14]=dot(temppt,u3);

			/* Reject tangential velocity */
			if (inotan != 0)
				f[i] = f[i] - dot(sp, f[i]) / dot(sp, sp) * sp;
			
			//rescale the total: should be small, for angular velocity constant
			f[i]=f[i]/scale[i];

			//total velocity and its fractional components
            num[15] = sqrt(dot(f[i],f[i]));
            num[16] = dot(f[i],u1);
            num[17] = dot(f[i],u2);
            num[18] = dot(f[i],u3);
			//Fsn mutual friction force density
			num[19] = dot(fD,u1);
			num[20] = dot(fD,u2);
			num[21] = dot(fD,u3);
            for (k=11;k<22;k++){
				if (k>11 && k!=15 && k<19 && abscomps)//only components
                	dif = abs(num[k]) - stats[k].m;
				else
                	dif = num[k] - stats[k].m;
                Rwt=dif*weight/tempwt;
                stats[k].m += Rwt;
                stats[k].q += sumweight*dif*Rwt;
            }
            sumweight=tempwt;//wait until now, the end, to do this
        #ifdef VELSTRUCT//Doesn't use abs(vel components)
            //Printing to VELSTRUCT starts here
            fvstruct << "(ivort,i,j)=(" << ivort << ",";
            fvstruct << i << "," << j << ") ";
            fvstruct << "(x,y,z)=(" << coretemp[j].x << ",";
            fvstruct << coretemp[j].y << "," << coretemp[j].z << ") ";
            for (k=0;k<22;k++){
                fvstruct << num[k];
				if (k!=21)
					fvstruct << " ";
				else
					fvstruct << endl;
			}
        #endif

			i++;
			j=periph[j].ip;
			spaceprev=space;
    	} while (j!=startcore[ivort].start);
		//end: for all points within each vortex
	}//end: for all vortices
	
	//Checking velocity statistics
    //avg = stats[i].m
    //var = stats[i].q/sumweight;
    //cout: (fvstruct prints only values, one line each quantity)   
    //standard deviation of preceding quantity is printed on even columns
    /*Column    Quantity (Averages of)
        [1]     [time]
        2       dot(tangent,u1)
        4       radius      
        6       beta
        8       |vlocal|
		10		vlocal-u1
		12		vlocal-u2
		14		vlocal-u3
        16      |v_nonlocal|
        18      v_nl-u1
        20      v_nl-u2
        22      v_nl-u3
        24      |v_fric|
        26      v_fric-u1
        28      v_fric-u2
        30      v_fric-u3
        32      |v_total|
        34      v_tot-u1
        36      v_tot-u2
        38      v_tot-u3
		40		Fsn.u1
		42		Fsn.u2
		44		Fsn.u3:no need to parallel transport? add comp's individually
		46		spacing (no std dev)
    */
    //Output time, from ifname (tinit)
	
    cout.precision(8);
    cout << tinit << " ";
    for (k=0;k<22;k++){
		if (k<19){
        	cout << stats[k].m << " " << sqrt(stats[k].q/sumweight) << " ";
		}else{//var for Fsn has factor of sumweight
			cout << (stats[k].m*sumweight/(2*M_PI*M_PI*r0*r0*r0)) << " ";
            cout << sqrt(stats[k].q/(2*M_PI*M_PI*r0*r0*r0)) << " ";
		}
	}
	cout << sumweight/(totalpts-1) << endl;
	
#ifdef VELSTRUCT
    fvstruct.close();
#endif
    return;
}

point localrz(point start, point mid, point end, point& tangent, double& rad,
        double& prefac)
{
    static point lp, lm, vel, curv, tmppt;
    static double lpdotp, lmdotm, lpdotm, norm;
    static double tmpdb, tmpdb1, tmpdb2;

    static point sumpt, combopt;
    static double denom, dp, dm;
    const static int curvfix=1;
                    //=1 remove component of curvature along r0, redo rad
                    //=0 don't adjust radius from 3sphere curvature
    lp = end - mid;
    lm = mid - start;

    lpdotm = dot(lp, lm);
    lpdotp = dot(lp, lp);
    lmdotm = dot(lm, lm);
    //New definition of tangent vector with error 
    //2nd order in neighboring arclength
    //Previous definition was 1st order
    tangent = lmdotm*lp + lpdotp*lm;
    //Remove component along radius
    tangent = tangent - (dot(tangent,mid)/(r0*r0))*mid;
    tangent = tangent/sqrt(dot(tangent,tangent));
    //same as newpt() and newfakept()
    if ((tmpdb=acos(lpdotm/sqrt(lpdotp*lmdotm)))>0.2){
        denom = lpdotp*lmdotm - lpdotm*lpdotm;
        sumpt=lp+lm;
        dp = lmdotm*dot(lp,sumpt);
        dm = lpdotp*dot(lm,sumpt);
        combopt=dp*lp - dm*lm;
        curv = 2*denom/dot(combopt,combopt)*combopt;
    }else{
        if (tmpdb<1e-15){
            curv.x=0.0; curv.y=0.0; curv.z=0.0; curv.w=0.0;
        }else{
            /*Linear in s_{+-}, derived by Aarts 3.9 & 3.10: 
            YIELDED BAD RESULTS AT HIGH ANGLES (truncation error), 
            EXCELLENT AT LOW ANGLES (low rounding error)*/
            curv=2*(lp/sqrt(lpdotp)-lm/sqrt(lmdotm))
                /sqrt(lpdotp+lmdotm+2*lpdotm);
        }
    }
    norm = dot(curv, curv);
	if (norm > 1e-20) { //below this approx value, tmpdb<1e-15
        vel = cross(tangent, curv, mid)/r0; //Divide out r0, from mid
        if (curvfix){   //correct radius for 3sphere curvature
            curv = curv-(dot(curv, mid)/(r0*r0))*mid;
            norm = dot(curv,curv);
            if (norm>1e-20)
                rad = 1/sqrt(norm);
            else
                rad = 1e10;
        }else{  //don't correct radius for 3sphere curvature
            rad=1/sqrt(norm);
        }
    }else{  //small curvature, small velocity
        vel.x = 0; vel.y = 0; vel.z = 0; vel.w = 0;
        if (curvfix)    //correct radius for 3sphere curvature
            rad = 1e10;
        else    //don't correct radius
            rad = r0; //indicates infinite radius, great circle
    }

    prefac = const1*log(2.*sqrt(sqrt(lpdotp*lmdotm))/const2);

    if (isnan(vel.x) || isnan(tangent.x) || isnan(rad)){
        cerr.precision(16);
        cerr << "localrz contains nan" << endl;
        cerr << "totalpts=" << totalpts << " nvort=" << nvort << endl;
        cerr << "start=" << start << endl;
        cerr << "mid=" << mid << endl;
        cerr << "end=" << end << endl;
        cerr << "lpdotp=" << lpdotp << " lpdotm=" << lpdotm;
        cerr << " lmdotm=" << lmdotm << endl;
        cerr << "rad=" << rad << " norm=" << norm << endl;
        cerr << "tangent=" << tangent << endl;
        cerr << "curv=" << curv << endl;
        exit(1);
    }

    return vel;
}

/* The point of dummy is so the compiler correctly identifies the class T. 
 This routine is first because I haven't figured out how to prototype template
 classes properly. */
template<class T> T readvar(T dummy, ifstream& ffrom, string varname)
{	//passes ffrom by reference	
	/*params must be in same order as readvar(), but some in 
	params.dat can be ignored in readvar() list.
	When using old params.dat or not, must be able to skip useless
	quantities.
	When using old params.dat need to skip quantities to get to NONLOCAL
	When using new params.dat just need to skip intermediate quantities*/
	string phrase; T thenum;
	//Not used for setting nonlocalon when "NONLOCAL on" detected

	//ifstream operator>> is the extraction operator
	while(ffrom.good()){//allows for skipping of varnames I don't want
					//set if badbit, failbit or eofbit
		ffrom >> phrase;	//ends at valid whitespace
		while (phrase[0]=='#') {
			if (phrase[1]=='#'){//reached the macros added to params.dat
				//go back (-) what was just put down (gcount) from "cur"(rent) 
				//position
				ffrom.seekg((streamoff)(-(int)ffrom.gcount()),ios_base::cur);
				return -1;//no parameters have negative val, use as error
			}else if (!ffrom.good()){
				if (ffrom.eof()){
					cerr << "Reached EOF in params file, in readvar\n";
				}else{
					cerr << "Either badbit or failbit set, in readvar\n";
				}
				exit(1);
			}
			getline(ffrom, phrase);//gets rest of line until \n
			ffrom >> phrase;//get beginning of new line until whitespace
		}
		if (phrase==varname) {//already took varname (string)
			ffrom >> phrase;	//then takes '=' sign
			ffrom >> thenum;	//then takes thenum
			break;
		}
	}
	if (!ffrom.good()){
		cerr << "Problem in params readvar finding "<< varname << endl;
		exit(1);
	}

	return thenum;
}

void readparams(string params, int argc) {
	int i;
	size_t found;
	string line;
	ifstream from(params.c_str());
	double cutoff_fac;//not needed as global variable in this program

	//in this particular version, variables must be in order 
	//but there can be variables in the params file that aren't 
	//present here
	vn = readvar(vn, from, (string) "vn");
	vs = readvar(vs, from, (string) "vs");
	ilocal = readvar(ilocal, from, (string) "ilocal");
	ifric = readvar(ifric, from, (string) "ifric");
	inotan = readvar(inotan, from, (string) "inotan");
	a0 = readvar(a0, from, (string) "a0");
	alpha = readvar(alpha, from, (string) "alpha");
	emax = readvar(emax, from, (string) "emax");
	r0 = readvar(r0, from, (string) "r0");
	//If others =-1, r0==-1 too:
	if (r0==-1){//never set/read r0 from params file
		cerr << "params file in readvar bad, r0=-1" << endl;
		exit(1);
	}
	if (argc==4 || argc==5){//Didn't include these in argv[] in main()
		absmax_fac = readvar(absmax_fac, from, (string)"absmax_fac");
		absmin_fac = readvar(absmin_fac, from, (string)"absmin_fac");
		cutoff_fac = readvar(cutoff_fac, from, (string)"cutoff_fac");
		limbiot_fac = readvar(limbiot_fac, from, (string)"limbiot_fac");
	}
	//look for NONLOCAL on, if so, set nonlocalon=1
	//ifstream operator>> is the extraction operator
	i=0;//if last one should end long before ##, if 
	while(from.good()){//found r0 properly, then exited
		getline(from,line);//retrieves entire line, unlike 
		if (line[1]=='#')//operator>> which ends at valid whitespace
			if(i==0)//2nd occurrence of ## at beginning of line
				i++;
			else
				break;
	}//#defined macros at end of params file
	nonlocalon=0;
	if (line[1]=='#'){//second occurrence of ## at beginning of line
		found=line.find("NONLOCAL on");//finds first occurrence of string
		if (found!=string::npos)//returned if string not found
			nonlocalon=1;
	}

	return;
}

void readcores(point core[], ptadmin periph[], ifstream& from) {
	int i, j, jmin, ivort;
	double hinit;	//normally global but useless here

	from >> totalpts; //totalpts; must not exceed NMAX in ubiq.h
	from >> nvort; //nvort; it must not exceed NVORTMAX in ubiq.h

	vortex * vlims = new vortex [nvort];
	// restart file overrides tinit from params.dat
	from >> tinit;
	from >> hinit;

	for (i = 0; i < nvort; i++) {
		from >> vlims[i].start;
		from >> vlims[i].end;
		from >> startcore[i].term;
	}

	for (ivort = 0; ivort < nvort; ivort++) {
		for (i = vlims[ivort].start; i != vlims[ivort].end; i++) {
			from >> j;
			from >> periph[j].recon;
			from >> core[j];
			if (i == vlims[ivort].start) {
				startcore[ivort].start = j;
			} else {
				periph[j].im = jmin;
				periph[jmin].ip = j;
			}
			jmin = j;
		}
		from >> j;
		from >> periph[j].recon;
		from >> core[j];
		startcore[ivort].end = j;
		periph[j].im = jmin;
		periph[jmin].ip = j;
		periph[j].ip = startcore[ivort].start;
		periph[startcore[ivort].start].im = j;
	}
	delete [] vlims;
	return;
}


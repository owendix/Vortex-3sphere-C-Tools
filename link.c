#include "ubiq.h"
/************************************
 * Computes circulation of vortex 1 on vortex 2 (note r0 value)
 * **********************************
 * Compile with linalg.c
 * Includes "ubiq.h" -> "linalg.h"
 * **********************************
 * Input: 
 * 		./a.out (restart)filename vortex 1 vortex 2
 * Output:
 * 		circulation(1->2)	circulation(2->1)
 * **********************************
 * */
void readcores(point *, char *, ptadmin *);
point biot(point &, point &, point &);

vortex startcore[NVORTMAX];
point core[NMAX];
ptadmin periph[NMAX];
int totalpts, nvort;
double r0=.05;		//3sphere radius in cm

int main(int argc, char *argv[]){
	int v1, v2;
	//Read in restart file
	if (argc==4){
		readcores(core, argv[1], periph);
		v1=atoi(argv[2]); v2=atoi(argv[3]);
		if (v1>=nvort || v2>=nvort || v1<0 || v2<0){
			cout << "vort1, vort2 must be in range: [0, ";
			cout << nvort-1 << endl;
			exit(0);
		} else if (v1==v2){
			cout << "vort1, vort2 must be distinct\n";
			exit(0);
		}
	} else {
		cout << "Include: ./a.out filename vort1 vort2\n";
		exit(0);
	}
	
	//Calculate circulation from vort2 (all) on vort1 (segment)
	//Repeat for all vort1 segments => compute circulation around vort1
	point vel, vellast, velstart, velavg;
	double circulate1=0.0;
	int i=startcore[v1].start, j;
	bool ok=1;
	while (ok){
		vel.x=0.0; vel.y=0.0;
		vel.z=0.0; vel.w=0.0;
		j=startcore[v2].start;
		//Calculate velocity of all vortex j on one point on vortex i
		vel = biot(core[j], core[periph[j].ip], core[i]);
		for (j=periph[j].ip; j!=startcore[v2].start; j=periph[j].ip)
			vel = vel + biot(core[j], core[periph[j].ip], core[i]);
		
		if (i==startcore[v1].start){
			//Save first point's velocity for end of while(ok)
			velstart=vel;
		}else{
			//vel is avg of begin and end of position vector: i->i+1
			velavg=0.5*(vel+vellast);
			//Compute increment of circulation: dot(vel, dl)
			circulate1=circulate1+dot(velavg,core[i]-core[periph[i].im]);
		}
		//Save for later
		vellast=vel;
		i=periph[i].ip;
		if (i==startcore[v1].start){
			//Last circulation increment: closes loop
			velavg=0.5*(vellast+velstart);
			circulate1=circulate1+dot(velavg,core[i]-core[periph[i].im]);
			ok=0;
		}
	}
	
	//Check:
	//Calculate circulation from vort1 (all) on vort2 (segment)
	//Repeat for all vort2 segments => compute circulation around vort2
	double circulate2=0.0;
	i=startcore[v2].start;
	ok=1;
	while (ok){
		vel.x=0.0; vel.y=0.0;
		vel.z=0.0; vel.w=0.0;
		j=startcore[v1].start;
		//Calculate velocity of all vortex j on one point on vortex i
		vel = biot(core[j], core[periph[j].ip], core[i]);
		for (j=periph[j].ip; j!=startcore[v1].start; j=periph[j].ip)
			vel = vel + biot(core[j], core[periph[j].ip], core[i]);
		
		if (i==startcore[v2].start){
			//Save first point's velocity for end of while(ok)
			velstart=vel;
		}else{
			//vel is avg of begin and end of position vector: i->i+1
			velavg=0.5*(vel+vellast);
			//Compute increment of circulation: dot(vel, dl)
			circulate2=circulate2+dot(velavg,core[i]-core[periph[i].im]);
		}
		//Save for later
		vellast=vel;
		i=periph[i].ip;
		if (i==startcore[v2].start){
			//Last circulation increment: closes loop
			velavg=0.5*(vellast+velstart);
			circulate2=circulate2+dot(velavg,core[i]-core[periph[i].im]);
			ok=0;
		}
	}
	
	//Two calculations should match?
	cout.setf(ios_base::showpoint);		//Allow trailing zeros
	cout << circulate1 << " " << circulate2 << endl;
	return 0;
}	

point biot(point & stseg, point & endseg, point & testpt)
{
	point l1, l2, l3, deldot, unitpt={1.0,0.0,0.0,0.0};
	double dot11, dot12, dot22, dot33, dot23, factor1, factor2;
	double tmp, root2=sqrt(2);
	//Calculates the 
	l1 = testpt - stseg;
	l2 = endseg - stseg;
	l3 = testpt - endseg;      
	
	//dot11=r0*acos(1-dot(l1,l1)/2/r0/r0);
	//dot33=r0*acos(1-dot(l3,l3)/2/r0/r0);
	//New way reduces catestrophic cancellation
	dot11 = r0*acos(dot(unitpt-l1/r0/root2, unitpt+l1/r0/root2));	//arclength
	dot33 = r0*acos(dot(unitpt-l3/r0/root2, unitpt+l3/r0/root2));	//arclength
	l1 = l1-dot(l1,stseg)*stseg/r0/r0;	//Find tangent vec at stseg
	l1 = l1*dot11/sqrt(dot(l1,l1));		//Normalize to arclength
	l3 = l3-dot(l3,endseg)*endseg/r0/r0;	//Find tangent vec at endseg
	l3 = l3*dot33/sqrt(dot(l3,l3));		//Normalize to arclength

	//Neglects curvature effects over distance mag(l2)	

	dot11 = dot11*dot11;	//need arclength squared
	dot33 = dot33*dot33;	//need arclength squared
	dot12 = dot(l1, l2);	//angle of vecs at approx stseg
	dot22 = dot(l2, l2);
	dot23 = dot12 - dot22;

	factor1 = dot11*dot22-dot12*dot12;

      if (fabs(factor1/dot11/dot22) > 1e-24)
      {
        factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
        deldot = (factor2 / factor1) * cross(l1, l2, stseg)/r0;
      }
      else
      {
        factor1 = dot22*dot33-dot23*dot23;
        if  (fabs(factor1/dot22/dot33) > 1e-24)
        {
            factor2 = dot23/sqrt(dot33) - dot12/sqrt(dot11);
            deldot = (factor2 / factor1) * cross(l3, l2, endseg)/r0;
        }
        else
        {
            if (dot11*dot33 > 1e-12) {
                factor2 = fabs(dot11-dot33)/dot11/dot33/sqrt(dot22);
                deldot = (factor2 / 2) * cross(l1, l2, stseg)/r0;
            }
            else {
                deldot.x = 0; deldot.y = 0; deldot.z = 0; deldot.w=0;
            }
        }
      }
      return deldot;
}

void readcores(point core[], char *file, ptadmin periph[]) {
	int i, j, jmin, ivort;
	double tinit, hinit;
	ifstream from(file);

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
		//numpts[i] = vlims[i].end-vlims[i].start+1;
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

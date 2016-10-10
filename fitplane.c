#include "ubiq.h"
/*************************************
 * Fits a plane to a vortex lying nearly on a great circle (note r0)
 * ***********************************
 * Compile with functors.c, linalg.c
 * ***********************************
 * Input:
 *		a.out filename r0(optional)
 *		filename is restart file from 3sphere containing
 *		vortices to be fit
 * ***********************************
 * Output:
 * 		vortex number
 * 		tangent plane vector 1 (t1)
 * 		tangent plane vector 2 (t2)
 * 		normal plane vector 1 (n1)
 * 		normal plane vector 2 (n2)
 * 		(point azimuthal ang)	n1-comp		n2-comp		|dev fm plane|
 * ***********************************
 * */ 

void readcores(point *, char *, ptadmin *);
point tangent(point, point, point);

vortex startcore[NVORTMAX];
point core[NMAX];
ptadmin periph[NMAX];
int totalpts, nvort, numpts[NVORTMAX];
double r0 = .05; //3sphere radius in cm

int main(int argc, char *argv[]){
	int i, j, ivort, maxnumpts=0;
	point vecs[4], tmppt;	//Four orthogonal vectors
	double t1comp, t2comp, n1comp, n2comp;
	double azim, orthog;
	//Read in vortex data from file option
	
	if (argc==2){
		readcores(core, argv[1], periph);
	}
	else if (argc==3){
		readcores(core, argv[1], periph);
		r0 = atof(argv[2]);		//Can specify r0 in cm
	} else {
		cout << "Include exactly 1 filename option\n";
		exit(0);
	}
	for (ivort=0;ivort<nvort; ivort++)
		if (numpts[ivort]>maxnumpts)
			maxnumpts=numpts[ivort];

	point *coreinv = new point[maxnumpts];
	point *taninv = new point[maxnumpts];
	point *rem = new point[maxnumpts];
	
	cout.precision(15);

	for (ivort=0; ivort<nvort; ivort++){
		//Find basis vecs for tangentplane and normalplane
		vecs[0].x=0; vecs[0].y=0;
		vecs[0].z=0; vecs[0].w=0;
		vecs[1].x=0; vecs[1].y=0;
		vecs[1].z=0; vecs[1].w=0;
		//find tangent vects
		for (i=0, j=startcore[ivort].start; i<numpts[ivort];
				i++, j=periph[j].ip){
			taninv[i] = tangent(core[periph[j].im], 
					core[j], core[periph[j].ip]);
			coreinv[i] = core[j]/r0;
			//Take the inversion of half the position
			//and tangent vectors
			if (i >= (numpts[ivort]-(numpts[ivort]%2))/2){
				taninv[i] = -1*taninv[i];
				coreinv[i] = -1*coreinv[i];
			}
			//Find vector sum of pos/tan vecs
			vecs[0] = vecs[0] + coreinv[i];
			vecs[1] = vecs[1] + taninv[i];
		}
		//Don't need to divide by numpts: will normalize = 1
		vecs[0] = vecs[0]/sqrt(dot(vecs[0],vecs[0]));
		//Check for orthogonality
		vecs[1] = vecs[1]-dot(vecs[1],vecs[0])*vecs[0];
		vecs[1] = vecs[1]/sqrt(dot(vecs[1],vecs[1]));

		//Go through euclidean basis vecs to find 3rd vector
		vecs[2].x=1.0; vecs[2].y=0.0;		
		vecs[2].z=0.0; vecs[2].w=0.0;		
		for (i=0;i<3;i++){
			tmppt=vecs[2]-dot(vecs[0],vecs[2])*vecs[0]
					-dot(vecs[1],vecs[2])*vecs[1];
			tmppt=tmppt/sqrt(dot(tmppt,tmppt));
			if (dot(tmppt,tmppt)<1.0e-8){
				//If orthogonal, go to next euclid vec
				if (i==0){
					vecs[2].x=0.0;
					vecs[2].y=1.0;
				}else if (i==1){
					vecs[2].y=0.0;
					vecs[2].z=1.0;
				}
			}else{
				vecs[2]=tmppt;
				break;
			}
		}

		vecs[3] = cross(vecs[0], vecs[1], vecs[2]);
		//Find components along vecs
		cout << ivort << endl;
		for (i=0; i<4; i++)
			cout << i << " " << vecs[i] << endl;
		for (i=0, j=startcore[ivort].start;i<numpts[ivort];
				i++, j=periph[j].ip){
			t1comp = dot(core[j],vecs[0])/r0;
			t2comp = dot(core[j],vecs[1])/r0;
			rem[i] = core[j]/r0 - t1comp*vecs[0] - t2comp*vecs[1];
		
			azim = atan(t2comp/t1comp);
			if (t1comp<=0)
				azim += M_PI;		
			else if(t1comp>=0 && azim < 0)
				azim += 2*M_PI;

			n1comp = dot(rem[i],vecs[2]);
			n2comp = dot(rem[i],vecs[3]);
			orthog = sqrt(n1comp*n1comp+n2comp*n2comp);
			cout << azim << " " << n1comp << " " << n2comp;
			cout << " " << orthog << endl;
		}

	}

	delete [] coreinv;
	delete [] taninv;
	delete [] rem;

	//Find vectors spanning plane of each vortex
	return 0;
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
		numpts[i] = vlims[i].end-vlims[i].start+1;
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
point tangent(point start, point mid, point end)
{
	point lp, lm, tan;
	double lpdotp, lmdotm;

	lp = end - mid;
	lm = mid - start;

	lpdotp = dot(lp, lp);
	lmdotm = dot(lm, lm);
	//New definition of tangent vector with error 2nd 
	//order in neighboring arclength
	//Previous definition was 1st order
	tan = lmdotm*lp + lpdotp*lm;
	tan = tan/sqrt(dot(tan,tan));

	return tan;
}

#include "arrangevorts.h"
/* ***************************************************
 This program arranges desired vortices from a 3sphere
 output file and prints them with point ordering
******************************************************
	Input (e.g.):
			./a.out in_filename out_filename v1 v2 v3 ...
			in_filename file is a restart file from program: 3sphere
	Output:
			out_filename file is formatted as a restart file
			which can be used to run 3sphere. It contains only 
			the vortices included upon program execution and the 
			index values are ordered
			
*/
void readcores(point *, char *, ptadmin *);
void writesort(point *, char *, ptadmin *);

vortex startcore[NVORTMAX];
point core[NMAX];
ptadmin periph[NMAX];
int totalpts, nvort, numpts[NVORTMAX], printvorts[NVORTMAX];
double tinit, hinit, r0 = .05; //3sphere radius in cm
int nvortout;

int main(int argc, char *argv[]){
	char *outfile;
	int i;
	//Read in vortex data from file option
	if (argc>3){
		readcores(core, argv[1], periph);
		outfile= argv[2];
		for (i=3,nvortout=0;i<argc;i++,nvortout++)
			printvorts[nvortout]=atoi(argv[i]);
	} else {
		cout << "Include: in_filename out_filename"; 
		cout << " ordered-nonduplicate-vortex-list\n";
		exit(0);
	}
	//Safety catches
	for (i=0;i<nvortout;i++){
		if (printvorts[i]<0 || printvorts[i]>(nvort-1)){
			cout << "Vortex range: (0," << (nvort-1)<< ")\n";
			exit(0);
		}
		if (i>0)
			if (printvorts[i]<=printvorts[i-1]){
				cout << "Vortices must be ordered, and unique\n";
				exit(0);
			}
	}
	//Write desired vortices in arranged order
	writesort(core, outfile, periph);
	
	return 0;
}

void readcores(point core[], char *file, ptadmin periph[]) {
	int i, j, jmin, ivort;
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

void writesort(point core[], char writefile[], ptadmin periph[]) {
	int i, j, k, count;
	ofstream tofile( writefile);

	count=0;
	for (i=0;i<nvortout;i++)
		count+=numpts[printvorts[i]];

	tofile.precision(15);
	/* file header should include number of points, nvort, time */
	tofile << count << "\n";
	tofile << nvortout << "\n";
	tofile << tinit << "\n";
	tofile << hinit << "\n";
	
	count=0;
	for (i = 0; i < nvortout; i++) { //now write these numbers of points
			tofile << count << " ";
			count=count+numpts[printvorts[i]]-1;
			tofile << count << " ";
			tofile << startcore[printvorts[i]].term;
			tofile << "\n";
			count++;
	}
	k=0;
	for (i = 0; i < nvortout; i++) { //now write the points themselves
		j = startcore[printvorts[i]].start;
		tofile << k << " " << periph[j].recon << " ";
		tofile << core[j] << "\n";
		while (j != startcore[printvorts[i]].end) {
			j = periph[j].ip;
			k++;
			tofile << k << " " << periph[j].recon << " ";
			tofile << core[j] << "\n";
		}
		k++;
	}
	return;
}

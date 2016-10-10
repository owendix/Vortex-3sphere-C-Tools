#include "ubiq.h"
#include <sstream>
//#define TESTROTATE
/**********************************************
Measures anisotropy of vortex tangles in trial. See 
bottom for what is printed. If isotropic, 
ipar = iperp1 = iperp2 = 2/3 and ibinpar=0. If lying 
within plane normal to velocity field, ipar=1, 
iperp1 =? iperp2 =? 1/2 and ibinpar is more complicated. 
The following relation should hold: ipar/2 + iperp1,2 = 1, 
I think.

Rotates binormal vector to r2nd's location to add together
for ibinpar.
***********************************************
Compile with linalg.c
Input: g++ -g anisotropy.c linalg.c -o isot.out
		./isot.out ../../code/mar1011a 10
-> use trial mar1011a only analyzing every 10th file
Output: (cout) 

time dist(optional) ipar iperp1 iperp2 |ibinpar| ipar/2+iperp1 ipar/2+iperp2
************************************************/
using namespace std;

void localrz(point, point, point, point&, point&);
void gramschmidt4d(point []);
point rotate4d(point, point*, double);

double r0;

int main(int argc, char* argv[]){

	int printdist=1;	//=1, print dist (c.f. ldens), =0, don't print
	
	int num, numlow, nskip, done, totalpts, nvort;
	int line, i, j, ivort, tmpint1, tmpint2, vlims[NVORTMAX];
	double a0=1.3e-8, volinv, angle;
	double dist, ipar, iperp1, iperp2;
	double time, tmpdb1, tmpdb2, tmpdist, fracdb;
	point t, b, r1st, r2nd, r, rm, rmm;
	point ibin, ibinpar, v, vp1, vp2, tmppt, tmppt1;
	point bvecs[4], ibinfrac, tmppt2;
	string ifnamemid, ifnamepre, ifnamebase, ifname;
	int loc;
	string strnum;
	stringstream strhandle;
	ifstream ifhandle;
#ifdef TESTROTATE
	point testvec;
	ofstream ftest("test.anisotropy.dat");
	ftest.precision(12);
#endif
	//directory where files are located	
	//ifnamepre="../../code/";

	if (argc==3){	//argv[argc] is null pointer
		nskip=atoi(argv[argc-1]);
		if (nskip<1)
			nskip=1;
		ifnamepre=(string)argv[1];
	}else if (argc==2){
		nskip=1;
		ifnamepre=(string)argv[1];
	}else{
		cerr << "Input requires 2 arguments: ./a.out filename nskip\n";
		exit(1);
	}
	//Concatenate filename to path
    if ((loc=ifnamepre.find_last_of("/"))==string::npos)
        perror("Error in find_last_of()\n");
    ifnamemid=ifnamepre;
    ifnamemid.assign(ifnamepre,loc,9);//filename with leading /

    ifnamebase=ifnamepre+ifnamemid; //concat with num to make ifname

	num=1;
	numlow=num;
	done=0;
	while (!done){
		//fracdb=0.0;
		dist=0.0;
		ipar=0.0;
		iperp1=0.0;
		iperp2=0.0;
		ibin.x=0.0;ibin.y=0.0;ibin.z=0.0;ibin.w=0.0;
		//ibinpar.x=0.0;ibinpar.y=0.0;ibinpar.z=0.0;ibinpar.w=0.0;
		strhandle << num;
		strnum=strhandle.str();
		ifname=ifnamebase+"."+strnum;
		strhandle.str("");	//clears it for next use
		ifhandle.open(ifname.c_str());
		if (!ifhandle.good()){
			done=1;
			cerr << "Last file checked: " << ifname << endl;
		}else{
			ifhandle >> totalpts;
			ifhandle >> nvort;
			ifhandle >> time;
			ifhandle >> tmpdb1;
			//get vortex limits
			for (ivort=0;ivort<nvort;ivort++){
				ifhandle >> tmpint1;
				ifhandle >> tmpint2;
				vlims[ivort]=tmpint2-tmpint1+1;	//num of pts per vortex
				ifhandle >> tmpint1;
			}
			//start going through pts
			line=0;
			for (ivort=0; ivort<nvort; ivort++){
				//first point
				ifhandle >> tmpint1; ifhandle >> tmpint1;
				ifhandle >> r1st;	//1st pt within vortex
				rm=r1st;
				//determine size of 3sphere
				if (num==numlow && ivort==0){
					r0=sqrt(dot(r1st, r1st));
					volinv=1/(2*M_PI*M_PI*r0*r0*r0);
				}
				line++;
				//r1st is actually i=0 pt, r2nd is i=1 pt
				for (i=1;i<vlims[ivort];i++){	//i: pt number within vortex
					ifhandle >> tmpint1; ifhandle >> tmpint1;
					ifhandle >> r;
					if (i==1){
						r2nd=r;	//2nd pt within vortex
						/*will be rotating vector b (binormal) to r2nd, to 
						compare at same point on manifold*/
						bvecs[0]=r2nd;
					}
					tmppt=r-rm;	//diff btw pt and prev pt
					tmpdb1=dot(tmppt,tmppt);
					tmpdist=r0*acos(1-tmpdb1/(2*r0*r0));
					dist += tmpdist;
					
					//compute anisotropy measures
					if (i>1){ //r, rm, rmm all stored by this iteration
						//associate rm with distance rm->r, calculate
						//tangent and binormal vectors at rm
						localrz(rmm, rm, r, t, b);	//sets tangent, t, and 
													//binormal, b
						//calculate applied velocity field at rm (unit vec)
						v.x=-rm.y/r0; v.y=rm.x/r0; v.z=-rm.w/r0; v.w=rm.z/r0;
						//calculate perpendicular directions to applied field
						vp1.x=rm.z/r0;vp1.y=-rm.w/r0;
						vp1.z=-rm.x/r0;vp1.w=rm.y/r0;
						vp2.x=-rm.w/r0;vp2.y=-rm.z/r0;
						vp2.z=rm.y/r0;vp2.w=rm.x/r0;
						//ipar
						tmpdb1=dot(t,v);
						ipar += (1-tmpdb1)*(1+tmpdb1)*tmpdist;
						//iperp1
						tmpdb1=dot(t,vp1);
						iperp1 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
						//iperp2
						tmpdb1=dot(t,vp2);
						iperp2 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
						//ibin
						//need to rotate b to r2nd
						
						//fracdb+=dot(b,v)/sqrt(dot(b,b));//use comp. paral. to v
						//somehow find fraction while integrating:
						//perhaps of magnitude		
						if (i!=2){//i==2: rm=r2nd: no need to rotate
							/*2nd basis vector: rotate b from rm to r2nd*/
							bvecs[1]=rm;
							//negative angle rotates fm bvecs[1] to bvecs[0]
							angle=-acos(dot(rm,r2nd)/(r0*r0));
							gramschmidt4d(bvecs);
							b=rotate4d(b,bvecs,angle);
							//tmppt1=rotate4d(tmppt1,bvecs,angle);
#ifdef TESTROTATE
							testvec=rm;
							ftest << "testvec= " << testvec;
							ftest << " r2nd= " << r2nd;
							testvec=rotate4d(testvec,bvecs,angle);
							ftest << " rot_testvec= " << testvec << endl;
#endif
						}
						ibin = ibin + b*tmpdist;
						//ibinpar, |v|=1
						//ibinpar = ibinpar + tmppt1*tmpdist;
					}
		
					//circulate old pts
					rmm=rm;
					rm=r;
					
					line++;	//cumulative pt number within tangle
				}
				//Complete ring: get last distance rlast->r1st
				//and calculate anisotropy at rm=rlast, and rm=r1st
				for (j=0;j<2;j++){
					if (j==0)
						r=r1st;
					else
						r=r2nd;
					tmppt=r-rm;	//diff btw pt and prev pt
					tmpdb1=dot(tmppt,tmppt);
					tmpdist=r0*acos(1-tmpdb1/(2*r0*r0));
					if (j==0)
						dist += tmpdist;	//close ring: with distance
					
					//compute anisotropy measures
					//tangent and binormal vectors at rm
					localrz(rmm, rm, r, t, b);	//sets tangent, t, and 
												//binormal, b
					//calculate applied velocity field at rm (unit vec)
					v.x=-rm.y/r0; v.y=rm.x/r0; v.z=-rm.w/r0; v.w=rm.z/r0;
					//calculate perpendicular directions to applied field
					vp1.x=rm.z/r0;vp1.y=-rm.w/r0;vp1.z=-rm.x/r0;vp1.w=rm.y/r0;
					vp2.x=-rm.w/r0;vp2.y=-rm.z/r0;vp2.z=rm.y/r0;vp2.w=rm.x/r0;
					//ipar
					tmpdb1=dot(t,v);
					ipar += (1-tmpdb1)*(1+tmpdb1)*tmpdist;
					//iperp1
					tmpdb1=dot(t,vp1);
					iperp1 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
					//iperp2
					tmpdb1=dot(t,vp2);
					iperp2 += (1-tmpdb1)*(1+tmpdb1)*tmpdist;	
				
					//fracdb+=dot(b,v)/sqrt(dot(b,b));//use comp. paral. to v
					//somehow find fraction while integrating:
					//perhaps of magnitude		
					//tmppt1=dot(b,v)*v;//use component parallel to v
					/*2nd basis vector: rotate b from rm to r2nd*/
					bvecs[1]=rm;
					//negative angle rotates fm bvecs[1] to bvecs[0]
					angle=-acos(dot(rm,r2nd)/(r0*r0));
					gramschmidt4d(bvecs);
					b=rotate4d(b,bvecs,angle);
					//tmppt1=rotate4d(tmppt1,bvecs,angle);
#ifdef TESTROTATE
					testvec=rm;
					ftest << "testvec= " << testvec;
					ftest << " r2nd= " << r2nd;
					testvec=rotate4d(testvec,bvecs,angle);
					ftest << " rot_testvec= " << testvec << endl;
#endif
					ibin = ibin + b*tmpdist;
					//ibinpar, |v|=1
					//ibinpar = ibinpar + tmppt1*tmpdist;
					//circulate old pts
					rmm=rm;
					rm=r;
				}
			}
			//|v|=1
			//ibin uses b
			tmpdb1=sqrt(dot(ibin,ibin));
			//ibinpar uses dot(b,v)*v: DOESN'T MAKE SENSE
			//what we ultimately want is for all the binormal vectors 
			//to be along the velocity field, on average. I'd need to 
			//think about the implications of the velocity field being
			//nonuniform. We can't find the component at the end.
			//tmpdb2=sqrt(dot(ibinpar,ibinpar));
			//fracdb/=totalpts;	//average fractional magnitude:
						//how much does b align with v

			ipar /= dist;
			iperp1 /= dist;
			iperp2 /= dist;
			dist *= volinv;
			tmpdb1 *= volinv/(dist*sqrt(dist)); //schwarz has /dist^(3/2)
			//tmpdb2 *= volinv/(dist*sqrt(dist));
	
			cout.precision(17);
	
			cout << time << " ";
			if (printdist)
				cout << dist << " ";
			cout.precision(7);
			cout << ipar << " " << iperp1;
			//print just component parallel to v,
			cout << " " << iperp2 << " " << tmpdb1 << " ";
			//then print what fraction this is of total
			//cout << fracdb << " ";			

			tmpdb1 = ipar/2 + iperp1;
			cout << tmpdb1 << " ";
			tmpdb1 = ipar/2 + iperp2;
			cout << tmpdb1 << endl;
			
			ifhandle.close();
			num += nskip;
		}
	}

#ifdef TESTROTATE
	ftest.close();
#endif
	return 0;
}

void localrz(point start, point mid, point end, point& tangent, point& vel){
	
	static point lp, lm, curv, tmppt, sumpt, combopt;
    static double lpdotp, lmdotm, lpdotm, norm;
    static double tmpdb1, tmpdb2, tmpdb;
	static double denom, dp, dm;
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
    }else{  //small curvature, small velocity
        vel.x = 0; vel.y = 0; vel.z = 0; vel.w = 0;
    }

    return;
}

void gramschmidt4d(point v[]){
    static int i;
    static point tmppt;
    /*Performs gram-schmidt process in 4D: requires first 2 basis vecs to be 
    set to some value, not equal to each other*/

    v[0] = v[0]/sqrt(dot(v[0],v[0]));   //unit vector to testpt
    v[1] = v[1]-dot(v[0],v[1])*v[0];            //remove component along v[0]
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
}

point rotate4d(point v, point bvecs[], double angle){
    static double bva[4][4], va[4], c, s, psum, Rij, tmpd[4];
    static point fv;
    static int i, j;
    /*Rotates vector v through angle in plane bvecs[0], bvecs[1]
    positive angle rotates from bvecs[0] to bvecs[1]*/

    //To iterate Rij, need to make point structs into arrays:
    //v[i].x,y,z,w = va[i][0,1,2,3]
    for (i=0;i<4;i++){
        bva[i][0]=bvecs[i].x;    //Notation va = v_array
        bva[i][1]=bvecs[i].y;
        bva[i][2]=bvecs[i].z;
        bva[i][3]=bvecs[i].w;
    }
    va[0]=v.x;   //Makes matrix mult easier
    va[1]=v.y;
    va[2]=v.z;
    va[3]=v.w;

    c=cos(angle);
    s=sin(angle);
    //Mult rotation matrix w pairwise summation over j
    for (i=0;i<4;i++){
        psum=0.0;
        for(j=0;j<4;j++){
            //produces the correct vectors: confirmed!
            Rij=((c*bva[0][i]+s*bva[1][i])*bva[0][j] + bva[2][i]*bva[2][j])
            + ((-s*bva[0][i]+c*bva[1][i])*bva[1][j] + bva[3][i]*bva[3][j]);
            /*  
            //from the current copy of nonlocal.pdf: has old, bad formula
            Rij = ((c*va[j][0]-s*va[j][1])*va[i][0] + va[i][2]*va[j][2])
                + ((s*va[j][0]+c*va[j][1])*va[i][1] + va[i][3]*va[j][3]);
            */
                    //pairwise summation
            psum += Rij*va[j];
            if (j==1){
                tmpd[i] = psum;
                psum=0.0;
            }
        }
        tmpd[i] +=psum;
    }
    fv.x = tmpd[0];
    fv.y = tmpd[1];
    fv.z = tmpd[2];
	fv.w = tmpd[3];

    return fv;
}

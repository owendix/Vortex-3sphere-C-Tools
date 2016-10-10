#include <cstdlib>
#include <ctime>
#include "linalg.h"
using namespace std;
/****************************************************************
*Test different ways to compute quantities within newpt:
*denoma compares current (05-08-11) method with several methods
*for various angle of separation and length difference between vectors
**Note, tried a few long doubles, minimal help unless I convert everything
**to long double
******************************************************************
*Compile: g++ -g newpttest.c linalg.c
*Input (run): ./a.out > outfile.dat 
*Parameters to adjust: 	eps  (|lp1|-|lm1|=eps)
						len=1e-5 (|lm1|=len)
*
*Output: (standard out) ./a.out
*printed columns:
*1: eps (|lp1-lm1|)
*2: ang (input)
*3: |p1curv| (current) w current combopt1
*4: cp1 method1: old way
*5: cm1
*6: denom1
*7: |p1curv|
*8: denom1 method2: (method1 w tweak1 to denom1)
*9: |p1curv|
*10: |p1curv| method3: old aarts way: SAID IT YIELDED AWFUL RESULTS, BUT I 
	GOT THAT THIS IS STABLE, AND BY FAR THE BEST!! NO NAN VALUES ANYWHERE,
	AND BEHAVES VERY WELL UNTIL ang=1e-15: AAH, IT gives a bit too big 
	a curvature vector for ang>~0.1 -> too small a radius, more nan's in the 
	equation to get the new point location on the curve!!! => Use old way, 
	above a certain angle (~0.2), and old aarts way (method3) below
11: |p1curv| method4: cleaner version of old way, check higher ang vals: 
	same high ang values as 7, excellent
12: ang(measured after rotation): CONFIRMED, EXCELLENT AGREEMENT
*********************************************************************/

double heron(double, double, double);

int main(){
	point tmppt;
	int newcombo=0;
	int i, j, k, n, flag, neps=101;
	double len=1e-5, eps=0, deps, ang=0, dang=1e-20, angmax=1;
	double scale=RAND_MAX+1.0, tmpdb;
	//rotation quantities
	point v[4];
	double Rij, psum, tmpd[4], tmpv[4], va[4][4], c, s;
	//newpt quantities
	point lp1, lm1, combopt1, p1curv, sumpt1;
	double lmdotm1, lpdotp1, lpdotm1, dp1, dm1, denoma;
	
	//Build difference vectors
	//Pick random lm1 direction
	srand((unsigned)time(0));	//Comment to keep same seed
	rand();		//Can't seed right before (?) produces bad numbers
	while(1){
		lm1.x=((double)rand()+(double)rand()/scale)/scale;	
		lm1.y=((double)rand()+(double)rand()/scale)/scale;	
		lm1.z=((double)rand()+(double)rand()/scale)/scale;	
		lm1.w=((double)rand()+(double)rand()/scale)/scale;
		if ((tmpdb=dot(lm1,lm1))>1.0e-12){
			//fix lm1 and lp1 to set length
			lm1=lm1*len/sqrt(tmpdb);
			break;
		}
	}
	lmdotm1=dot(lm1, lm1);
	//to rotate lp1, need plane in which to rotate: must contain lm1
	v[0]=lm1/sqrt(lmdotm1);
	tmppt.x=1.0; tmppt.y=0.0; tmppt.z=0.0; tmppt.w=0.0;
	flag=0;
	for (i=0;i<3;i++){
		tmppt=tmppt-dot(tmppt,v[0])*v[0];
		if (dot(tmppt,tmppt)<1.0e-8){
			if (flag==0){
				tmppt.x=0.0; tmppt.y=1.0;//tmppt, can only be parallel v[0] once
			}else{
				tmppt.y=0.0; tmppt.z=1.0;
			}
		}else{
			if (flag==0){
				v[1]=tmppt/sqrt(dot(tmppt,tmppt));
				flag=1;
			}else{
				v[2]=tmppt/sqrt(dot(tmppt,tmppt));
				break;
			}
		}
	}
	v[3]=cross(v[0],v[1],v[2]);
	//vary length difference
	deps=len/(neps-1);
	for (k=0;k<neps;k++){
		eps=deps*(k-neps/2);
		//rotate lp1, wrt lm1, by small angles
		n=0;
		ang=0.0;
		while (ang<angmax){
			//reinitialize lp1 so error doesn't pile up
			lp1=lm1*(len+eps)/sqrt(dot(lm1,lm1));	//lp1 unrotated initially
			if (n>0){
				//ang=dang on first entry here
		    	//v[i].x,y,z,w = va[i][0,1,2,3]
		    	for (i=0;i<4;i++){
		    	    va[i][0]=v[i].x;    //Notation va = v_array
	    		    va[i][1]=v[i].y;
	    		    va[i][2]=v[i].z;
	    		    va[i][3]=v[i].w;
	    		}
			
	    		tmpv[0]=lp1.x;   //Makes matrix mult easier
	    		tmpv[1]=lp1.y;
	    		tmpv[2]=lp1.z;
	    		tmpv[3]=lp1.w;
		
		    	c=cos(ang);
		    	s=sin(ang);
				//Mult rotation matrix w pairwise summation over j
		    	for (i=0;i<4;i++){
		    	    psum=0.0;
		    	    for(j=0;j<4;j++){
		    	    //Confirmed, accurate to 1e-15, 02-09-11: 
					//cpptools/rotate4d.cpp
		        	    Rij=((c*va[0][i]+s*va[1][i])*va[0][j]+va[2][i]*va[2][j])
		        	    +((-s*va[0][i]+c*va[1][i])*va[1][j]+va[3][i]*va[3][j]);
		        	    //pairwise summation
		        	    psum += Rij*tmpv[j];
		        	    if (j==1){
		        	        tmpd[i] = psum;
		        	        psum=0.0;
		        	    }
		        	}
		        	tmpd[i] +=psum;
		    	}
		    	lp1.x = tmpd[0];
		    	lp1.y = tmpd[1];
		    	lp1.z = tmpd[2];
		    	lp1.w = tmpd[3];
	
			}//end if(cnt>0)
	
			//calculate useful quantities common to both methods
		 	lpdotp1 = dot(lp1, lp1);
			//lmdotm1 already calc'd, doesn't change
		    lpdotm1 = dot(lp1, lm1);
			//current (05-08-11) combopt1
			sumpt1=lp1+lm1;
			dp1 = lmdotm1*dot(lp1,sumpt1);
	    	dm1 = lpdotp1*dot(lm1,sumpt1);
			combopt1=dp1*lp1 - dm1*lm1;
	
			cout.precision(16);
			cout << eps << " " << ang << " ";
	
			//now test newpt quantities
			//current (05-08-11) denoma, original combopt1
			denoma = lpdotp1 * lmdotm1 - lpdotm1 * lpdotm1;	
			p1curv = (2*denoma/dot(combopt1,combopt1))*combopt1;
			tmpdb=sqrt(dot(p1curv,p1curv));

			cout << tmpdb << " ";	

			//method1: old way (feb2010)
			double cp1, cm1, denom1;
			//denoma = lpdotp1*lmdotm1 - lpdotm1*lpdotm1;
			cp1 = .5 * lmdotm1 * (lpdotp1 + lpdotm1) / denoma;
    		cm1 = .5 * lpdotp1 * (lmdotm1 + lpdotm1) / denoma;	
			denom1=cp1*cp1*lpdotp1+cm1*cm1*lmdotm1-2*cp1*cm1*lpdotm1;
			p1curv=(cp1/denom1)*lp1 - (cm1/denom1)*lm1;
			tmpdb=sqrt(dot(p1curv,p1curv));
		
			cout << cp1 << " " << cm1 << " " << denom1 << " " << tmpdb << " ";

			//method2: old way w/ tweak
			tmppt=cp1*lp1-cm1*lm1;
			denom1=dot(tmppt,tmppt);
			p1curv=(cp1/denom1)*lp1 - (cm1/denom1)*lm1;
			tmpdb=sqrt(dot(p1curv,p1curv));
		
			cout << denom1 << " " << tmpdb << " ";
			
			//method3: aarts derived way "YIELDED AWFUL RESULTS" means 
			//above ang=0.1:
			//WRONG, YIELDED GREAT RESULTS below ang=0.1
			p1curv=2*(lp1/sqrt(lpdotp1)-lm1/sqrt(lmdotm1))
				/sqrt(lpdotp1+lmdotm1+2*lpdotm1);
			tmpdb=sqrt(dot(p1curv,p1curv));

			cout << tmpdb << " ";

			//method4: my cleaner, less stable, version of old way
			denoma = lpdotp1 * lmdotm1 - lpdotm1 * lpdotm1;
			sumpt1=lp1+lm1;
			dp1 = lmdotm1*dot(lp1,sumpt1);
			dm1 = lpdotp1*dot(lm1,sumpt1); 
			combopt1=dp1*lp1 - dm1*lm1;
			p1curv = 2*denoma/dot(combopt1,combopt1)*combopt1;
			tmpdb=sqrt(dot(p1curv,p1curv));

			cout << tmpdb << " ";

			cout << acos(lpdotm1/sqrt(lpdotp1*lmdotm1));
			//use as measure of when invalid (ang<1e-15)

			cout << endl;

			//iterate for next time
			if (n==0)
				ang=dang;		//first pass, ang=dang;
			else
				ang*=2;
		
			n++;
		}//end while(ang<angmax)
	}//end for (k=0;k<neps;k++)

	return 0;
}

double heron(double tmp1, double tmp2, double tmp3){
	double a, b, c, a2;

	if (tmp1>tmp2){
		if (tmp1>tmp3){
			a=tmp1;
			if (tmp2>tmp3){
				b=tmp2; c=tmp3;
			}else{
				b=tmp3; c=tmp2;	
			}
		}else{
			a=tmp3;	b=tmp1; c=tmp2;
		}
	}else{
		if (tmp2>tmp3){
			a=tmp2;
			if (tmp1>tmp3){
				b=tmp1; c=tmp3;
			}else{
				b=tmp3; a=tmp1;
			}
		}else{
			a=tmp3; b=tmp2; c=tmp1;
		}
	}

	//parentheses essential
	a2=(a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))/16;

	return a2;
}

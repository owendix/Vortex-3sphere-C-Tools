#include "linalg.h"

ostream& operator<< (ostream& outwhere, const point& p)
{
	outwhere << p.x <<" "<<p.y<<" "<<p.z<<" "<<p.w;
	return outwhere;
}

istream& operator>> (istream& infrom, point& p)
{
	infrom >> p.x;
	infrom >> p.y;
	infrom >> p.z;
	infrom >> p.w;
	return infrom;
}

int operator== (point a, point b)
{
	int test=1;
	if (a.x != b.x) test=0;
	if (a.y != b.y) test=0;
	if (a.z != b.z) test=0;
	if (a.w != b.w) test=0;
	return test;
}

ostream& operator>> (ostream& outwhere, const plane& p){
	outwhere << p.pnum << " " << p.p1 << " " << p.p2;
	return outwhere;
}

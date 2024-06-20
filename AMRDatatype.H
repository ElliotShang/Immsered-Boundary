#ifndef _AMRDatatype_H_
#define _AMRDatatype_H_

struct S_newtag
{
	int remotebkp;
	int newtag;

	S_newtag()
	{}

	S_newtag(const int & r0, const int & n0):
	remotebkp(r0),newtag(n0)
	{}
};

struct S_tag
{
	int aimcell;
	int localbkp;
#ifdef DEBUG
	Point bxyz;
#endif		

	S_tag()
	{}

	S_tag(const int & a0, const int & l0):
	aimcell(a0),localbkp(l0)
	{}

};

struct S_RemoteBKPInfo
{
	int bkpindex;
	int localcell;
#ifdef DEBUG
	Point bxyz;
#endif		

	S_RemoteBKPInfo()
	{}

	S_RemoteBKPInfo(const int & li, const int & rc)
	{
		bkpindex = li;
		localcell = rc;
	}
};

struct S_newfinebp0
{
	int aimbkp;
	Boxson<int> son0;
#ifdef DEBUG
	Boxson<Point> sonpt;
#endif		

	S_newfinebp0()
	{}

	S_newfinebp0(const int & a0, const Boxson<int> & fi0):
	aimbkp(a0),son0(fi0)
	{}

};

struct S_newcoarsebp
{	
	int remotebkp;
	int localcell;

	S_newcoarsebp()
	{}

	S_newcoarsebp(const int & r0, const int & l0):
	remotebkp(r0),localcell(l0)
	{}
};

struct Block_dis
{
	int remotebkp;
	int b0;
	int p0;
	double distance;
	double signdis;

	Block_dis()
	{}

	Block_dis(const int & rmbkp, const int & bd, const int & pat0, const double & dis, const double & sdis):
	remotebkp(rmbkp),b0(bd),p0(pat0),distance(dis),signdis(sdis)
	{}
};

struct BoxIndexBC
{
	Point pt;
	Pointxyz bc;
	BoxIndexBC()
	{}

	BoxIndexBC(Point & mypt, Pointxyz & mybc):
	pt(mypt), bc(mybc)
	{}


};

#endif
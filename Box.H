#ifndef _BOX_H_
#define _BOX_H_

#include "Point.H"
//#include "Boxloc.H"
#include "Point_iterator.H"
#include "Constants.H"
#include "AMRSpaceTime.H"
#include "FlowVariables.H"

struct rftag
{
  int tag;
  int detag;

  rftag()
  {
    tag = -1;
    detag = -1;
  }

  inline void givetag()
  {
    tag = 0;
  }

  inline void givedetag()
  {
    detag = 0;
  }

  inline void Init()
  {
    tag = -1;
    detag = -1;
  }

  bool istag()
  {
    if (tag > -1)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  bool isdetag()
  {
    if (detag > -1)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
};

enum Boxtype
{
  Notype, Dmghost, Blockghost, Normalcell
};

struct BoxtoWall
{
  int body;
  int patch;
  double distance;
  double signdis;
  double distance_to_dm;
  double distance_to_body;

  BoxtoWall()
  {
      body = -1;
      patch = -1;
      distance = 9999.0;
      signdis = 9999.0;
      distance_to_dm = 9999.0;
      distance_to_body = 9999.0;
  }

  void init()
  {
      body = -1;
      patch = -1;
      distance = 9999.0;
      signdis = 9999.0;
      distance_to_dm = 9999.0;
      distance_to_body = 9999.0;
  }

  inline void switchtwo(BoxtoWall & abw)
  {
    if (abw.distance < distance)
    {
      BoxtoWall temp_abw = *this;
      *this = abw;
      abw = temp_abw;
    }
  }
};

enum Pairtype
{
  none, f_pro, f_res1, f_res2, ghost_pro, ghost_res1, ghost_res2
};

struct Box
{

  Point             lowpt;
  int        neib[3][3][3];
  int pts[2][2][2];
  int faces[3][2];
  Boxtype type;
  Pairtype ptype;
  int bkpid;
  bool solid;
  BoxtoWall pair;

	Box()
	{
		lowpt = Point(-1,-1,-1);
    for (Point_iterator p(0,2); p.end() ; ++p)
    {
      pts[p.i][p.j][p.k] = -1;
    }
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        faces[i][j] = -1;
      }
    }
    for (Point_iterator p(0,3); p.end(); ++p)
    {
      neib[p.i][p.j][p.k] = -1;
    }
    type = Notype;
    solid = false;
    ptype = none;
    bkpid = -1;
	}

	Box (Point & lp):
	lowpt(lp)
	{	
    for (Point_iterator p(0,2); p.end() ; ++p)
    {
      pts[p.i][p.j][p.k] = -1;
    }
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        faces[i][j] = -1;
      }
    }
    for (Point_iterator p(0,3); p.end(); ++p)
    {
      neib[p.i][p.j][p.k] = -1;
    }
    type = Notype;
    ptype = none;
    bkpid = -1;
	}
  Box (const int &i, const int &j, const int &k):
  lowpt(Point(i,j,k))
  {
    for (Point_iterator p(0,2); p.end() ; ++p)
    {
      pts[p.i][p.j][p.k] = -1;
    } 
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        faces[i][j] = -1;
      }
    }
    for (Point_iterator p(0,3); p.end(); ++p)
    {
      neib[p.i][p.j][p.k] = -1;
    }
    type = Notype;
    solid = false;
    ptype = none;
    bkpid = -1;
  }

  inline Point & boxpt()
  {
    return lowpt;
  }

  inline void init_pairtype()
  {
    ptype = none;
  }

  void setneib(const int & i0, const int & j0, const int & k0, const int & neib0)
  {
    neib[i0][j0][k0] = neib0;
  }
  void setneib(const int & i0, const int & j0, const int & k0, const int * neib0)
  {
    neib[i0][j0][k0] = *neib0;
  }
  // void setneib(const int & i0, const int & j0, const int & k0, const int & level0, const int & index0)
  // {
  //   neib[i0][j0][k0].level = level0;
  //   neib[i0][j0][k0].index = index0;
  // }

	inline const int & ix() const
	{
		return lowpt.xy[0];
	}
	inline const int & iy() const
	{
		return lowpt.xy[1];
	}
	inline const int & iz() const
	{
		return lowpt.xy[2];
	}

  // inline double p(const int & i0) const
  // {
  //   return double(lowpt.xy[0]+i0)*dh[neib[1][1][1].level][0];
  // }
  // inline double pointy(const int & j0) const
  // {
  //   return double(lowpt.xy[1]+j0)*dh[neib[1][1][1].level][1];
  // }
  // inline double pointz(const int & z0) const
  // {
  //   return double(lowpt.xy[2]*bcnz+z0)*dh[neib[1][1][1].level][2];
  // }

  inline int & point(const int & i0, const int & j0, const int & k0)
  {
    return pts[i0][j0][k0];
  }

  // bool hasfineneib()
  // {
  //   bool fineflag = false;
  //   for (Point_iterator p(0,3); p.end(); ++p)
  //   {
  //     if (neib[p.i][p.j][p.k].level == neib[1][1][1].level+1)
  //     {
  //       fineflag = true;
  //       break;
  //     }
  //   }
  //   return fineflag;
  // }

  // bool hascoarseneib()
  // {
  //   bool coarseflag = false;
  //   if (neib[1][1][1].level > 0)
  //   {
  //     for (Point_iterator p(0,3); p.end(); ++p)
  //     {
  //       if (neib[p.i][p.j][p.k].level == neib[1][1][1].level-1)
  //       {
  //         coarseflag = true;
  //         break;
  //       }
  //     }
  //   }
  //   return coarseflag;
  // }

  void initpoint()
  {
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int k = 0; k < 2; ++k)
        {
          pts[i][j][k] = -1;
        }
      }
    }
  }

  void init_boxtowall()
  {
    pair.init();
  }
};

struct BoxlocChange
{
  int oldloc;
  int newloc;
  
  BoxlocChange()
  {}

  BoxlocChange(int opt, int newpt)
  {
    oldloc = opt;
    newloc = newpt;
  }
};

// struct BoxtoWall_Pair
// {
//   BoxtoWall pair[2];

//   BoxtoWall_Pair()
//   {}

//   void init()
//   {
//     pair[0].init();
//     pair[1].init();
//   }
// };

struct HGCell
{
  Pointxyz pt;
  Pointxyz attach_pt;
  double hgdis;
  int closecell;
  double closedis;
#if DIM == 4  
  int intpcell[4];
  double distointp[4];
  double intpcoef[4];
#endif 
#if DIM < 4
  int intpcell[2][2][2];
  double distointp[2][2][2];
  double intpcoef[2][2][2];
#endif  
  FlowVariables fv;
#ifdef PASSAGE_ANGLE
  double rotangle;
#endif    
#if DIM == 4
  void Comptintpcoef()
  {
    double revs[4];
    double dissum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      revs[i] = 1.0/(distointp[i]*distointp[i]);
      dissum += revs[i];
    }
    for (int i = 0; i < 4; ++i)
    {
      intpcoef[i] = revs[i]/dissum;
      if (isnan(intpcoef[i])) intpcoef[i] = 1.0;
    }
  }
#endif
#if DIM < 4
  void Comptintpcoef()
  {
    double revs[2][2][2];
    double dissum = 0.0;
    for (Point_iterator p(0,2); p.end(); ++p)
    {
        revs[p.i][p.j][p.k] = 1.0/(distointp[p.i][p.j][p.k]*distointp[p.i][p.j][p.k]);
        dissum += revs[p.i][p.j][p.k];
    }
    for (Point_iterator p(0,2); p.end(); ++p)
    {
        intpcoef[p.i][p.j][p.k] = revs[p.i][p.j][p.k]/dissum;  
    }
  }

  void InitIntparray()
  {
    for (Point_iterator p(0,2); p.end(); ++p)
    {
      intpcell[p.i][p.j][p.k] = -1;
    }
  }
#endif  
};

struct IBCell
{
  int ci;
  /*touch point*/
  Pointxyz tp;
  /*Flow variables*/
  FlowVariables fv;
  double pre_grad;
  Pointxyz u_grad;

  HGCell hgc;

  int mynode;

  double yplus;
  double uplus;
  double ut;

  bool closetowall;
  int wallbox;

  double hg_yplus;
  double hg_vt;

  Pointxyz tangdir;
  double hg_viseddy;
  Pointxyz vn;
  Pointxyz patv;

  double angle_ib_to_pat;


};

struct CellMCoef
{
  Pointxyz var[3];
  CellMCoef()
  {}
//   CellMCoef(const double & kx, const double & ky, const double & kz,
//             const double & ix, const double & iy, const double & iz,
//             const double & zx, const double & zy, const double & zz)
//   {
//     var[0][0] = kx; var[0][1] = ky; var[0][2] = kz;
//     var[1][0] = ix; var[1][1] = iy; var[1][2] = iz;
//     var[2][0] = zx; var[2][1] = zy; var[2][2] = zz;
//   }
  inline Pointxyz & operator[](const int & i0)
  {
    return var[i0];
  }
  inline void zero()
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        var[i][j] = 0.0;
      }
    }
  }

  bool operator==(const CellMCoef & amc)
  {
    bool same = true;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        if (abs(var[i].xyz[j] - amc.var[i].xyz[j]) > 0.000001)
        {
          same = false;
        }
      }
    }
    return same;
  }

  void showdata(const string & mcstr)
  {
    printf("mcs are %f,%f,%f,%f,%f,%f,%f,%f,%f,tag %s\n",
      var[0][0],var[0][1],var[0][2],var[1][0],
      var[1][1],var[1][2],var[2][0],var[2][1],var[2][2],mcstr.c_str());
  }
//   inline double & at(const int & i0, const int & j0)
//   {
//     return var[i0][j0];
//   }
};

struct BoxCellGeom
{
  Pointxyz                 boxcenter;
  CellMCoef                keisa;
  double                   v;

  inline double & operator[](const int & i0)
  {
    return boxcenter[i0];
  }
};

struct Surfpatch
{
#ifdef TRI_ELEMENT  
  int corpt[3];
#endif
#ifdef QUA_ELEMENT
  int corpt[4];
#endif    
  /*patch center*/
  Pointxyz pc;
  /*normal vector*/
  Pointxyz nv;
  /*area*/
  double area;
  bool reverse_nv_dir;
  HGCell hgc;
  double ut;
  double hg_vt;
  double yplus;
  double wall_T;
  double Cp, Cf;
  int node;
  
  double filt_pre_input[FORDER];
  double filt_pre_output[FORDER];

  Surfpatch()
  {
    reverse_nv_dir = false;
    ut = -99.0;
    yplus = -99.0;
    wall_T = -99.0;
    Cp = -99.0;
    Cf = -99.0;
    node = -1;
  }

  inline double & operator[](const int & index)
  {
    return pc[index];
  }

  inline double patdis(Surfpatch & apat)
  {
    return pow(pow(apat[0]-pc[0],2)+pow(apat[1]-pc[1],2)+pow(apat[2]-pc[2],2),0.5);
  }

};

#endif

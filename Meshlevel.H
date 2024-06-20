//
//  Meshlevel.H
//  AMR
//
//  Created by wangzhuo on 2019/11/23.
//  Copyright © 2019 wangzhuo. All rights reserved.
//

#ifndef Meshlevel_h
#define Meshlevel_h

#include <vector>
#include "Box.H"
#include "AMRmpi.H"
#include "DataArray.H"
#include "FlowVariables.H"
#include "Face.H"
#include "GhostPair.H"
#include "Nodedata.H"
#include "MeshOp.H"

using namespace std;

#define RenewtheBoxFace(oldbox, newbox)\
  {\
    for (int i202 = 0; i202 < 3; ++i202)\
    {\
      for (int d204 = 0; d204 < 2; ++d204)\
      {\
        int f206 = m_box[oldbox].faces[i202][d204];\
        if (f206 > -1)\
        {\
          Assert(newbox > -2, "The new location of the box must be non-negative!!!",32);\
          m_face[f206][1-d204] = newbox;\
        }\
      }\
    }\
  }

#define EndConnection(bn)\
  {\
    for (Point_iterator p(0,3); p.end(); ++p)\
    {\
      int an0 = m_box[bn].neib[p.i][p.j][p.k];\
      if (an0 > -1)\
      {\
        int rvi = 2-p.i;\
        int rvj = 2-p.j;\
        int rvk = 2-p.k;\
        if (m_box[an0].neib[rvi][rvj][rvk] == bn)\
        {\
          m_box[an0].neib[rvi][rvj][rvk] = -1;\
        }\
      }\
      m_box[bn].neib[p.i][p.j][p.k] = -1;\
    }\
  }  

class Meshlevel
{
  friend class AMR;
  friend class Mesh;
  friend class Pltoutput;
  friend class NS_Solver;
  friend class BCValues;
  friend class TwoLevelOp;
  friend class Body;

public:
  Meshlevel()
  {}

  Meshlevel(vector<Box> & a_box, vector<mPoint> & a_point, vector<Face> & a_face, const int & dmgn)
  {
    m_box = DataArray<Box>(a_box, dmgn);
    m_point = DataArray<mPoint>(a_point, 0);
    m_face = DataArray<Face>(a_face, 0);
    GiveAFlag("Finish construct the Meshlevel for the back mesh!!!", 5);
  }

  void clearlevel()
  {
    m_box.cleararray();
    m_data.cleararray();
    m_geom.cleararray();
    m_point.cleararray();
    m_face.cleararray();
  }

  // void DataExchange();

  void Assignnewloc(const bool & holeflag);

  //void RenewtheBoxFace(const int & oldbox, const int & newbox);

  void EndMeshHoleConnection();

  //void EndConnection(const int & bn);

  void RemoveOldPoints(const int & ilevel);

  void RemoveOldFace(const int & ilevel);

  void RenewBoxFace(const int & ilevel);

  void CheckGhost(const int & ilevel);

  void TagRemovedBox();

  void RemoveConnections();

  void RenewFacesBox();

  void InitPair()
  {
    int bs = m_box.ps();
    int be = m_box.pe();
    for (int i = bs; i < be; ++i)
    {
      m_box[i].ptype = none;
    }
    bs = m_box.gps();
    be = m_box.gpe();
    for (int i = bs; i < be; ++i)
    {
      m_box[i].ptype = none;
    }
    MPI_Barrier(share_comm);
  }

  void NewFaceArea(Face & aface, const int & mylevel, const int & fd, const int & ibox, const int & f0)
  {
#ifndef IMPORT_MESH             
    ComptFaceArea(mylevel, aface);
    if (fd == 0)
    {
      ComptFaceCenter(aface, m_point[m_box[ibox].pts[f0][0][0]].xyz, 
                             m_point[m_box[ibox].pts[f0][1][0]].xyz, 
                             m_point[m_box[ibox].pts[f0][1][1]].xyz,
                             m_point[m_box[ibox].pts[f0][0][1]].xyz);
    }
    else if (fd == 1)
    {
      ComptFaceCenter(aface, m_point[m_box[ibox].pts[f0][0][0]].xyz, 
                             m_point[m_box[ibox].pts[f0][1][0]].xyz, 
                             m_point[m_box[ibox].pts[f0][1][1]].xyz,
                             m_point[m_box[ibox].pts[f0][0][1]].xyz);
    }
    else if (fd == 2)
    {
      ComptFaceCenter(aface, m_point[m_box[ibox].pts[0][0][f0]].xyz, 
                             m_point[m_box[ibox].pts[1][0][f0]].xyz, 
                             m_point[m_box[ibox].pts[1][1][f0]].xyz,
                             m_point[m_box[ibox].pts[0][1][f0]].xyz);
    }
#endif
#ifdef IMPORT_MESH
    if (fd == 0)
    {
      Assert(m_box[ibox].pts[f0][0][0] > -1, "Box point x1 must be non-negative!!!", 119);
      Assert(m_box[ibox].pts[f0][1][0] > -1, "Box point x2 must be non-negative!!!", 120);
      Assert(m_box[ibox].pts[f0][1][1] > -1, "Box point x3 must be non-negative!!!", 121);
      Assert(m_box[ibox].pts[f0][0][1] > -1, "Box point x4 must be non-negative!!!", 122);
      ComptFaceArea_FourPts(aface, m_point[m_box[ibox].pts[f0][0][0]].xyz, 
                                                    m_point[m_box[ibox].pts[f0][1][0]].xyz, 
                                                    m_point[m_box[ibox].pts[f0][1][1]].xyz,
                                                    m_point[m_box[ibox].pts[f0][0][1]].xyz);
    }
    else if (fd == 1)
    {
      Assert(m_box[ibox].pts[0][f0][0] > -1, "Box point y1 must be non-negative!!!", 119);
      Assert(m_box[ibox].pts[0][f0][1] > -1, "Box point y2 must be non-negative!!!", 120);
      Assert(m_box[ibox].pts[1][f0][1] > -1, "Box point y3 must be non-negative!!!", 121);
      Assert(m_box[ibox].pts[1][f0][0] > -1, "Box point y4 must be non-negative!!!", 122);
      ComptFaceArea_FourPts(aface, m_point[m_box[ibox].pts[0][f0][0]].xyz, 
                                                    m_point[m_box[ibox].pts[0][f0][1]].xyz, 
                                                    m_point[m_box[ibox].pts[1][f0][1]].xyz,
                                                    m_point[m_box[ibox].pts[1][f0][0]].xyz);
    }
    else if (fd == 2)
    {
      Assert(m_box[ibox].pts[0][0][f0] > -1, "Box point z1 must be non-negative!!!", 119);
      Assert(m_box[ibox].pts[1][0][f0] > -1, "Box point z2 must be non-negative!!!", 120);
      Assert(m_box[ibox].pts[1][1][f0] > -1, "Box point z3 must be non-negative!!!", 121);
      Assert(m_box[ibox].pts[0][1][f0] > -1, "Box point z4 must be non-negative!!!", 122);
      ComptFaceArea_FourPts(aface, m_point[m_box[ibox].pts[0][0][f0]].xyz, 
                                                    m_point[m_box[ibox].pts[1][0][f0]].xyz, 
                                                    m_point[m_box[ibox].pts[1][1][f0]].xyz,
                                                    m_point[m_box[ibox].pts[0][1][f0]].xyz);
    }
#endif
  }

  void IncludeNewFace(vector<Face> & nface)
  {
    int nfnum = nface.size();
    vector<int> nfnum_procs(sprocs);
    vector<int> nfnum_start(sprocs);
    int totnfnum = 0;
    int existface = m_face.realsize();
    MPI_Reduce(&nfnum, &totnfnum, 1, MPI_INT, MPI_SUM, 0, share_comm);
    MPI_Bcast(&totnfnum, 1, MPI_INT, 0, share_comm);
    MPI_Allgather(&nfnum, 1, MPI_INT, &nfnum_procs[0], 1, MPI_INT, share_comm);
    ArrayProcsStart(nfnum_procs, nfnum_start);
    PRINTFinLEVELRANK0("The total new face is %d", -1, totnfnum);
    for (int i = 0; i < nface.size(); ++i)
    {
      int i0 = i + nfnum_start[srank] + existface;
      for (int f0 = 0; f0 < 2; ++f0)
      {
        int b0 = nface[i][f0];
        int fn0 = nface[i].fnv;
#ifdef DEBUG
        Assert(b0 > -1, "The new face side box must be non-negative!!!", 161);
        if (m_box[b0].faces[fn0][1-f0] > 0)
        {
          printf("Box %d (%d,%d,%d) dir %d side %d is face %d current face is %d\n", b0,
            m_box[b0].ix(), m_box[b0].iy(), m_box[b0].iz(), fn0, 1-f0, m_box[b0].faces[fn0][1-f0], i0);
          MPI_Abort(MPI_COMM_WORLD, 165);
        }
#endif              
        m_box[b0].faces[fn0][1-f0] = i0;
      }
    }
    m_face.Addnew(nface);
    m_face.DirectlyReduceNew();
    MPI_Barrier(share_comm);
    //printf("new face number is %d\n", (int)nface.size());
  }

private:
  DataArray<Box>                                  m_box;
  DataArray<FlowVariables>                        m_data; 
  DataArray<BoxCellGeom>                          m_geom;
  DataArray<mPoint>                               m_point;
  DataArray<Face>                                 m_face;
  DataArray<NodePair>                             blockpair;
  DataArray<Domainghost> dghost;
  Nodedata                                        nd;
  /*For flux modification*/
  // DataArray<int>                            m_cneib;
  // DataArray<ShellFaceData<FlowVec> >              m_neibflux;
};

#endif
/* Meshlevel_h */

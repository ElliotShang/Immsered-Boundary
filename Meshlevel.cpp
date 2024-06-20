#include "Meshlevel.H"

	// void Meshlevel::DataExchange()
 //  {
 //    MPI_Win mywin;
 //    Assert(nd.recvdestcell.size()==nd.recvfv.size(), "Error in data exchange!!!", 5);
 //    //printf("in DataExchange!!! recvdestcell num %d recvfv size %d\n", (int)nd.recvdestcell.size(), (int)nd.recvfv.size());
 //    //double tk[6];
 //    //tk[0] = MPI_Wtime();
 //    int s0 = 0;
 //    for (int n0 = 0; n0 < nodenum; ++n0)
 //    {
 //      int nfvsize = nd.num2ranks[n0];
 //      //printf("node %d to node %d is %d recv is %d\n", node, n0, nfvsize, nd.numrecv[n0]);
 //      //int s0 = nd.d_cell_extrac[n0].ps();
 //      //Assert(nfvsize==nd.d_cell_extrac[n0].pe()-s0, "Error in data exchange!!!",14);
 //      // for (int i = 0; i < nfvsize; ++i)
 //      // {
 //      //   nd.flowmesg[n0][i] = m_data[nd.d_cell_extrac[n0][i+s0].localcell];
 //      // }
 //      for (int i = 0; i < nfvsize; ++i)
 //      {
 //        nd.flowmesg[i+s0] = m_data[nd.cell_extrac[n0][i].localcell];
 //      }
 //      s0 += nfvsize;
 //    }
 //    MPI_Barrier(MPI_COMM_WORLD);
 //    //printf("Get the flow mes...\n");
 //    int tot_send_num = nd.flowmesg.size();
 //    int have_recv_num = 0;
 //    //tk[1] = MPI_Wtime();   
 //    MPI_Win_create(&nd.flowmesg[0], tot_send_num*sizeof(FlowVariables), sizeof(FlowVariables), 
 //      MPI_INFO_NULL, nodecomm, &mywin);
 //    //MPI_Win_fence(0, mywin[i]);
 //    //MPI_Win_lock_all(0, nd.nodewin);
 //    MPI_Win_lock_all(0, mywin);
 //    // tk[2] = MPI_Wtime(); 
 //    //Assert(nd.recvfv.size() > 0,"Error in data exchange!!!", 35);    
 //    for (int i = 0; i < nodenum; ++i)
 //    {
 //      if (nd.numrecv[i] > 0)
 //      {
 //        //MPI_Get(&nd.recvfv[have_recv_num], nd.numrecv[i], MPI_FV, i, nd.getstart[i], nd.numrecv[i], MPI_FV, nd.nodewin);
 //        MPI_Get(&nd.recvfv[have_recv_num], nd.numrecv[i], MPI_FV, i, nd.getstart[i], nd.numrecv[i], MPI_FV, mywin);
 //        have_recv_num += nd.numrecv[i];
 //      }
 //    }
 //    //tk[3] = MPI_Wtime();
 //    MPI_Win_unlock_all(mywin);
 //    MPI_Win_free(&mywin);
 //    //MPI_Win_unlock_all(nd.nodewin);
 //    //GiveAFlag("Finish get remote data!!!", 5);
 //    //tk[4] = MPI_Wtime();   
 //    int recvsize = nd.recvfv.size();
 //    //printf("N%dR%d total recvsize is %d\n", node, srank, recvsize);
 //    for (int i = 0; i < recvsize; ++i)
 //    {
 //      int i0 = nd.recvdestcell[i];
 //      m_data[i0] = nd.recvfv[i];
 //      // if (m_box[i0].ix() == 241)
 //      // {
 //      //   printf("box (%d,%d,%d)\n", m_box[i0].ix(), m_box[i0].iy(), m_box[i0].iz());
 //      //   m_data[i0].showdata("recv fv");
 //      // }
 //    }
 //    //tk[5] = MPI_Wtime();
 //    MPI_Barrier(MPI_COMM_WORLD);
 //    //printf("1 %f 2 %f 3 %f 4 %f 5 %f\n", 
 //    //  tk[1]-tk[0], tk[2]-tk[1], tk[3]-tk[2], tk[4]-tk[3], tk[5]-tk[4]);
 //  }

void Meshlevel::Assignnewloc(const bool & holeflag)
{
  int lowc[3] = {0, 0, 1};
  int oldi0, newi0;
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
  if (holeflag)
  {
    for (int i = 0; i < m_box.hole.size(); ++i)
    {
      int i0 = m_box.hole[i];
      Assert(m_box[i0].neib[1][1][1] == -1, "Box array hole error!!!", 81);
    } 
    MPI_Barrier(MPI_COMM_WORLD);
    int d_hn = m_box.hole.size();
    MPI_Allreduce(MPI_IN_PLACE, &d_hn, 1, MPI_INT, MPI_SUM, share_comm);
    int bs, be;
    ArrayOrder_s(0, m_box.realsize(), bs, be, sprocs, srank);
    int d_hn0 = 0;
    for (int i = bs; i < be; ++i)
    {
      if (m_box[i].neib[1][1][1] == -1)
      {
        ++d_hn0;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &d_hn0, 1, MPI_INT, MPI_SUM, share_comm);
    if (d_hn0 != d_hn)
    {
      PRINTFinLEVEL("The hole number is %d but the removed box bunber is %d", -1, d_hn, d_hn0);
      MPI_Abort(MPI_COMM_WORLD, 103);
    }
  }
  else
  {
    int gnum1 = m_box.arrayghostnum();
    int gnum2 = 0;
    int bs, be;
    ArrayOrder_s(0, m_box.realsize(), bs, be, sprocs, srank);
    for (int i = bs; i < be; ++i)
    {
      if (m_box[i].neib[1][1][1] != i)
      {
        PRINTFinLEVEL("The array index is not right!!! Box %d (%d,%d,%d) index is %d",-1,i,
          m_box[i].ix(), m_box[i].iy(), m_box[i].iz(), m_box[i].neib[1][1][1]);
        MPI_Abort(MPI_COMM_WORLD,117);
      }
    }
    MPI_Barrier(share_comm);
    for (int i = 0; i < m_box.ghost_index.size(); ++i)
    {
      int i0 = m_box.ghost_index[i];
      if (m_box[i0].neib[1][1][1] != i0)
      {
        PRINTFinLEVEL("Ghost box %d (%d,%d,%d) index error index is %d!!!", -1, i0, m_box[i0].ix(), m_box[i0].iy(), m_box[i0].iz(), m_box[i0].neib[1][1][1]);
        MPI_Abort(MPI_COMM_WORLD, 113);
      }
      else
      {
        m_box[i0].neib[1][1][1] = -1;;
      }
    }
    MPI_Barrier(share_comm);
    for (int i = bs; i < be; ++i)
    {
      if (m_box[i].neib[1][1][1] == -1)
      {
        ++gnum2;
        m_box[i].neib[1][1][1] = i;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE ,&gnum2, 1, MPI_INT, MPI_SUM, share_comm);
    if (gnum1 != gnum2)
    {
      PRINTFinLEVEL("Sum of ghost index is %d but ghost number is %d", -1, gnum1, gnum2);
      MPI_Abort(MPI_COMM_WORLD, 138);
    }
    MPI_Barrier(share_comm);

  }
 #endif   
  m_box.holeplan(holeflag);
  vector<BoxlocChange> & loctag = m_box.locchange();
  MPI_Barrier(MPI_COMM_WORLD);
  vector<vector<vector<vector<int> > > > oldneib0(loctag.size(), 
    vector<vector<vector<int> > >(3, 
    vector<vector<int> >(3, 
    vector<int>(3, -1))));
  vector<vector<vector<vector<int> > > > oldneib1(loctag.size(), 
    vector<vector<vector<int> > >(3, 
    vector<vector<int> >(3, 
    vector<int>(3, -1))));
  for (int i = 0; i < loctag.size(); ++i)
  {
    for (Point_iterator p(0,3); p.end(); ++p)
    {
      oldneib0[i][p.i][p.j][p.k] = m_box[loctag[i].oldloc].neib[p.i][p.j][p.k];
    }
  }
  if (!holeflag)
  {
    for (int i = 0; i < loctag.size(); ++i)
    {
      for (Point_iterator p(0,3); p.end(); ++p)
      {
        oldneib1[i][p.i][p.j][p.k] = m_box[loctag[i].newloc].neib[p.i][p.j][p.k];
      }
    }
  }
  /*This barrier is necessary!!!*/
  MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Win_fence(0, m_box.arraywin());
  for (int i = 0; i < loctag.size(); ++i)
  {
    oldi0 = loctag[i].oldloc;
    newi0 = loctag[i].newloc;
    //printf("Node %d old %d new %d\n",node, oldi0, newi0);
    Assert(oldi0 != newi0, "Two boxes location is same!!!", 83);
    Assert(m_box[oldi0].neib[1][1][1] == oldi0, "The location of the box has been changed!!!", 83);
    Assert(newi0 > -1, "The new location of the box must be positive!!!", 88); 
    for (Point_iterator p(0,3); p.end(); ++p)
    {
#ifdef DEBUG
      if (holeflag)
      {
        if (m_box[newi0].neib[p.i][p.j][p.k] != -1)
        {
          printf("The hole box (%d,%d,%d) neib (%d,%d,%d) is %d not -1\n", 
            m_box[newi0].ix(),m_box[newi0].iy(),m_box[newi0].iz(),p.i,p.j,p.k,m_box[newi0].neib[p.i][p.j][p.k]);
          MPI_Abort(MPI_COMM_WORLD, 204);
        }
      }
#endif          
      int neib0 = oldneib0[i][p.i][p.j][p.k];
      int rvny = 2-p.j; int rvnz = 2-p.k;
      int rvnx = 2-p.i;
      if (neib0 > -1)
      {
#ifdef DEBUG        
        if (m_box[neib0].neib[rvnx][rvny][rvnz] != oldi0)
        {
          printf("Box %d neib (%d,%d,%d) is %d but its reverse neib is %d\n", oldi0,p.i,p.j,p.k,neib0,m_box[neib0].neib[rvnx][rvny][rvnz]);
        }
        Assert(m_box[neib0].neib[rvnx][rvny][rvnz] == oldi0, "The old neib error 1", 204);
#endif        
        m_box[neib0].neib[rvnx][rvny][rvnz] = newi0;
      }             
    }
    RenewtheBoxFace(oldi0, newi0);
  }
  if (!holeflag)
  {
    for (int i = 0; i < loctag.size(); ++i)
    {
      oldi0 = loctag[i].oldloc;
      newi0 = loctag[i].newloc;
      Assert(oldi0 != newi0, "Two ghost boxes location is same!!!", 115);
      Assert(m_box[newi0].neib[1][1][1] == newi0, "The new location of the ghost has been changed!!!", 114);
      Assert(newi0 > -1, "The new location of the box must be positive!!!", 118);
      Assert(oldi0 > -1, "The new location of the box must be positive!!!", 120);
      for (Point_iterator p(0,3); p.end(); ++p)
      {
        int neib0 = oldneib1[i][p.i][p.j][p.k];
        int rvny = 2-p.j;
        int rvnz = 2-p.k;
        int rvnx = 2-p.i;
        if (neib0 > -1)
        {
#ifdef DEBUG
          if (m_box[neib0].neib[rvnx][rvny][rvnz] != newi0)
          {
            printf("Box %d neib (%d,%d,%d) is %d but its reverse neib is %d\n", 
              newi0,p.i,p.j,p.k,neib0,m_box[neib0].neib[rvnx][rvny][rvnz]);
          }
          Assert(m_box[neib0].neib[rvnx][rvny][rvnz] == newi0, "The old neib error 2", 239);   
#endif
          m_box[neib0].neib[rvnx][rvny][rvnz] = oldi0;
        }
      }
      RenewtheBoxFace(newi0, oldi0);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void Meshlevel::TagRemovedBox()
{
  for (int i = 0; i < m_box.hole.size(); ++i)
  {
    int i0 = m_box.hole[i];
    Assert(m_box[i0].neib[1][1][1] != -1, "The hole box main index has also ready been set to -1!!!", 268);
    m_box[i0].neib[1][1][1] = -1;
  }
  MPI_Barrier(share_comm);
}

void Meshlevel::RemoveConnections()
{
  for (int i = 0; i < m_box.hole.size(); ++i)
    {
      int i0 = m_box.hole[i];
      RenewtheBoxFace(i0, -1);
      EndConnection(i0);
    }
    MPI_Barrier(share_comm);
}

  void Meshlevel::EndMeshHoleConnection()
  {

    for (int i = 0; i < m_box.hole.size(); ++i)
    {
      int i0 = m_box.hole[i];
      Assert(m_box[i0].neib[1][1][1] != -1, "The hole box main index has also ready been set to -1!!!", 268);
      m_box[i0].neib[1][1][1] = -1;
      // PRINTFinLEVEL("The hole box is (%d,%d,%d)",-1,m_box[i0].ix(),m_box[i0].iy(),m_box[i0].iz());
      RenewtheBoxFace(i0, -1);
      EndConnection(i0);
    }
    MPI_Barrier(share_comm);
  }

  void Meshlevel::RemoveOldPoints(const int & ilevel)
  {
    int ss[2] = {m_box.ps(), m_box.gps()};
    int se[2] = {m_box.pe(), m_box.gpe()};
    int mpps = m_point.ps();
    int mppe = m_point.pe();
    for (int i = mpps; i < mppe; ++i)
    {
      m_point[i].index = -1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ci = 0; ci < 2; ++ci)
    {
      for (int i = ss[ci]; i < se[ci]; ++i)
      {
        for (Point_iterator p(0,2); p.end(); ++p)
        {
          int ptn = m_box[i].pts[p.i][p.j][p.k];
          Assert(ptn > -1, "Negative point index in ExcludeOldPoints!!!", 752);
          Assert(ptn < m_point.size(), "Point not in range in ExcludeOldPoints!!!", 753);
          m_point[ptn].index = ptn;
        }
      }
    }
    Assert(m_point.Holenum() == 0, "Hole number is nonzero before removing!!!", 201);
    MPI_Barrier(share_comm);
    for (int i = mpps; i < mppe; ++i)
    {
      if (m_point[i].index < 0)
      {
        m_point.givehole(i);
      }
    }
    PRINTFinLEVELRANK0("Point number before CompressArray %d hole number %d", 
      ilevel, (int)m_point.size(), (int)m_point.Holenum());
    m_point.holeplan(true);
    for (int lc = 0; lc < m_point.newloctag.size(); ++lc)
    {
      int oldid = m_point.newloctag[lc].oldloc;
      int newid = m_point.newloctag[lc].newloc;
      m_point[oldid].index = newid;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ci = 0; ci < 2; ++ci)
    {
      for (int bn0 = ss[ci]; bn0 < se[ci]; ++bn0)
      {
        for (Point_iterator p(0,2); p.end(); ++p)
        {
          int pti = m_box[bn0].pts[p.i][p.j][p.k];
          if (pti < 0)
          {
            printf("A negative index for L%dB%d box (%d,%d,%d)\n", 
              ilevel,bn0,m_box[bn0].ix(),m_box[bn0].iy(),m_box[bn0].iz());
            MPI_Abort(MPI_COMM_WORLD, 856);
          }
          else
          {
            m_box[bn0].pts[p.i][p.j][p.k] = m_point[pti].index;
#ifdef DEBUG
            if (m_box[bn0].pts[p.i][p.j][p.k] < 0)
            {
              printf("A negative index for L%dB%d box (%d,%d,%d)\n", 
              ilevel,bn0,m_box[bn0].ix(),m_box[bn0].iy(),m_box[bn0].iz());
              MPI_Abort(MPI_COMM_WORLD,872);
            }
#endif          
          }
        }
      }
    }
    m_point.CompressArray();
    PRINTFinLEVELRANK0("Point number after CompressArray %d", ilevel, (int)m_point.size());
  }

  void Meshlevel::RemoveOldFace(const int & ilevel)
  {
    for (int i = m_face.ps(); i < m_face.pe(); ++i)
    {
      if (m_face[i][0] == -1 || m_box.isghost(m_face[i][0]))
      {
        if (m_face[i][1] == -1 || m_box.isghost(m_face[i][1]))
        {
          m_face.givehole(i);
          for (int f0 = 0; f0 < 2; ++f0)
          {
          	if (m_face[i][f0] > -1)
          	{
          		m_box[m_face[i][f0]].faces[m_face[i].fnv][1-f0] = -1;
          	}
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    m_face.holeplan(true);
    PRINTFinLEVELRANK0("Face number before CompressArray %d hole number %d", 
      ilevel, (int)m_face.size(), (int)m_face.Holenum());
    m_face.CompressArray();
    MPI_Barrier(MPI_COMM_WORLD);
    PRINTFinLEVELRANK0("Face number after CompressArray %d", ilevel, (int)m_face.size());
    RenewBoxFace(ilevel);
#ifdef DEBUG
    for (int i = m_box.gps(); i < m_box.gpe(); ++i)
    {
      for (int fd = 0; fd < 3; ++fd)
      {
        for (int fi = 0; fi < 2; ++fi)
        {
          if (m_box[i].faces[fd][fi] > -1)
          {
            int s0 = m_box[i].faces[fd][fi];
            if (!m_box.isnormal(m_face[s0][fi]))
            {
              PRINTFinLEVEL("Ghost box %d (%d,%d,%d) face (%d,%d) is %d not -1!!!", ilevel,
                i, m_box[i].ix(), m_box[i].iy(), m_box[i].iz(),
                fd, fi, m_box[i].faces[fd][fi]);
              MPI_Abort(MPI_COMM_WORLD, 393);
            }
          }
        }
      }
    }
#endif       
  }

  void Meshlevel::RenewBoxFace(const int & ilevel)
	{
    Point adp[3] = {Point(2,1,1),Point(1,2,1),Point(1,1,2)};
    int mfps = m_face.ps();
    int mfpe = m_face.pe();
		for (int i = mfps; i < mfpe; ++i)
		{
			for (int fn = 0; fn < 2; ++fn)
			{
				int bn = m_face[i][fn];
#ifdef DEBUG
        if (bn < 0)
        {
          Point adp(1,1,1);
          *(&adp[m_face[i].fnv]) = 2*fn;
          int fb = m_face[i][1-fn];
          int fb_n = m_box[fb].neib[adp[0]][adp[1]][adp[2]];
          PRINTFinLEVEL("Face %d side %d is box %d another side box is %d isghost %d"
            " neib to the face is %d the neib 1 face is %d neib 2 box face is %d neib 2 isghost %d", ilevel, 
            i, fn, bn, fb, m_box.isghost(fb), fb_n,
            m_box[fb].faces[m_face[i].fnv][fn],
            m_box[fb_n].faces[m_face[i].fnv][1-fn],
            m_box.isghost(fb_n));
        }
        Assert(bn>-1,"The face must have two boxes!!!",372);
#endif        
				m_box[bn].faces[m_face[i].fnv][1-fn] = i;
			}
#ifdef DEBUG
      Point & xy1 = m_box[m_face[i][0]].lowpt;
      Point & xy2 = m_box[m_face[i][1]].lowpt;
      int a10 = *(&xy1[m_face[i].fnv]);
      int a20 = *(&xy2[m_face[i].fnv]);
      if (a20 - a10 != 1 && !(periodic[m_face[i].fnv] && a20 + a10 + 1 == (highpt.xy[m_face[i].fnv] - lowpt.xy[m_face[i].fnv])*level_grid_ratio[ilevel][m_face[i].fnv]))
      {
        PRINTFinLEVEL("Face %d side 0 is box (%d,%d,%d) side 1 is box (%d,%d,%d)!!! in RenewBoxFace",ilevel,i,
          xy1[0],xy1[1],xy1[2],xy2[0],xy2[1],xy2[2]);
        MPI_Abort(MPI_COMM_WORLD,445);
      }
#endif      
		}
    //MPI_Win_fence(0, m_box.arraywin());
		MPI_Barrier(share_comm);
		//MPI_Win_fence(0, m_geom.arraywin());
		//MPI_Win_flush_all(m_geom.arraywin());	
	}

  void Meshlevel::RenewFacesBox()
  {
    MPI_Win_fence(0,m_face.arraywin());
    for (int i = m_face.ps(); i < m_face.pe(); ++i)
    {
      for (int fi = 0; fi < 2; ++fi)
      {
        int i0 = m_face[i][fi];
        if (i0 > -1)
        {
          m_face[i][fi] = m_box[i0].neib[1][1][1];
        }
      }
    }
    MPI_Win_fence(0,m_face.arraywin());
    MPI_Barrier(MPI_COMM_WORLD);
  }
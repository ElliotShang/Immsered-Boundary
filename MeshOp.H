#ifndef _MeshOp_H_
#define _MeshOp_H_
/*This function was used in the Mesh.H to compute the new cell geom*/
#define ComptCellParams(cellmc, ilevel, ibox)\
		cellmc.v = dh[ilevel][0]*dh[ilevel][1]*dh[ilevel][2];\
		cellmc.boxcenter.setvalue(0.0, 0.0, 0.0);\
		cellmc.boxcenter[0] = ((double)m_level[ilevel].m_box[ibox].ix()+0.5)*dh[ilevel][0];\
		cellmc.boxcenter[1] = ((double)m_level[ilevel].m_box[ibox].iy()+0.5)*dh[ilevel][1];\
		cellmc.boxcenter[2] = ((double)m_level[ilevel].m_box[ibox].iz()+0.5)*dh[ilevel][2];\
		cellmc.boxcenter.get_theta();\
		cellmc.keisa[0] = Pointxyz(1.0/dh[ilevel][0], 0.0, 0.0);\
		cellmc.keisa[1] = Pointxyz(0.0, 1.0/dh[ilevel][1], 0.0);\
		cellmc.keisa[2] = Pointxyz(0.0, 0.0, 1.0/dh[ilevel][2]);\
		// printf("Rank %d Level %d box %d center is %f,%f,%f\n", nrank, ilevel, ibox, \
		// 	cellmc.boxcenter[0], cellmc.boxcenter[1], cellmc.boxcenter[2]);

#define ComptFaceArea(ilevel, newface)\
			if (newface.fnv == 0)\
			{\
				newface.area = dh[ilevel][1]*dh[ilevel][2];\
				newface.keisa.setvalue(1.0, 0.0, 0.0);\
				newface.keisa1.setvalue(0.5, 0.0, 0.0);\
				newface.keisa2.setvalue(0.5, 0.0, 0.0);\
			}\
			else if (newface.fnv == 1)\
			{\
				newface.area = dh[ilevel][2]*dh[ilevel][0];\
				newface.keisa.setvalue(0.0, 1.0, 0.0);\
				newface.keisa1.setvalue(0.0, 0.5, 0.0);\
				newface.keisa2.setvalue(0.0, 0.5, 0.0);\
			}\
			else if (newface.fnv == 2)\
			{\
				newface.area = dh[ilevel][0]*dh[ilevel][1];\
				newface.keisa.setvalue(0.0, 0.0, 1.0);\
				newface.keisa1.setvalue(0.0, 0.0, 0.5);\
				newface.keisa2.setvalue(0.0, 0.0, 0.5);\
			}

#define ComptFaceCenter(newface, pt0, pt1, pt2, pt3)\
			newface.fc1 = (pt1+pt0+pt3)/3.0;\
			newface.fc2 = (pt1+pt2+pt3)/3.0

#define ComptFaceArea_FourPts(newface, pt0, pt1, pt2, pt3)\
			Pointxyz pt66 = pt1 - pt0;\
			Pointxyz pt67 = pt3 - pt0;\
			Pointxyz pt68 = pt1 - pt2;\
			Pointxyz pt69 = pt3 - pt2;\
			ComptFaceCenter(newface, pt0, pt1, pt2, pt3);\
			newface.keisa1 = pt66.cross(pt67)*0.5;\
			newface.keisa2 = pt69.cross(pt68)*0.5;\
			newface.keisa = newface.keisa1+newface.keisa2;\
			newface.area = newface.keisa.length();\
			newface.keisa /= newface.area;\
			newface.keisa1 /= newface.area;\
			newface.keisa2 /= newface.area

#define ComptFaceNmv_ThreePts(nmv, pt0, pt1, pt2)\
			Pointxyz pt66 = pt1 - pt0;\
			Pointxyz pt67 = pt2 - pt0;\
			nmv = pt66.cross(pt67);\
			nmv.normalize()

// #define ComptCellParams_Pts_Mesh(ilevel,ibox)\
// 			m_level[ilevel].m_geom[ibox].boxcenter.setvalue(0.0,0.0,0.0);\
// 			for (Point_iterator p45(0,2); p45.end(); ++p45)\
// 			{\
// 				int pttt0 = m_level[ilevel].m_box[ibox].pts[p45.i][p45.j][p45.k];\
// 				Assert(pttt0 > -1, "The box point index must be non-negative when computing the mesh geom!!!",49);\
// 				m_level[ilevel].m_geom[ibox].boxcenter += m_level[ilevel].m_point[pttt0].xyz*0.125;\
// 			}\
// 			Pointxyz p2p44 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
// 											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][1]].xyz;\
// 			Pointxyz p2p46 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][1]].xyz -\
// 											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz;\
// 			Pointxyz p2p48 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
// 											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][0]].xyz;\
// 			Face aface90;\
// 			Pointxyz fareakeisa71[3][2];\
// 			for (int fd62 = 0; fd62 < 3; ++fd62)\
// 			{\
// 				for (int fi64 = 0; fi64 < 2; ++fi64)\
// 				{\
// 					Face * aimface;\
// 					int face66 = m_level[ilevel].m_box[ibox].faces[fd62][fi64];\
// 					if (face66 > -1)\
// 					{\
// 						aimface = &m_level[ilevel].m_face[face66];\
// 					}\
// 					else\
// 					{\
// 						if (fd62 == 0)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][0][0]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][1][0]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][1][1]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][0][1]].xyz);\
// 						}\
// 						else if (fd62 == 1)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][fi64][0]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][fi64][1]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][fi64][1]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][fi64][0]].xyz);\
// 						}\
// 						else if (fd62 == 2)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][fi64]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][fi64]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][fi64]].xyz,\
// 															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][1][fi64]].xyz);\
// 						}\
// 						aimface = &aface90;\
// 					}\
// 					fareakeisa71[fd62][fi64] = aimface->keisa*aimface->area;\
// 				}\
// 			}\
// 			m_level[ilevel].m_geom[ibox].v = (abs(fareakeisa71[0][0].dot(p2p44))+\
// 																				abs(fareakeisa71[1][1].dot(p2p46))+\
// 																				abs(fareakeisa71[2][0].dot(p2p48)))/3.0;\
// 			for (int di = 0; di < 3; ++di)\
// 			{\
// 				m_level[ilevel].m_geom[ibox].keisa[di] = \
// 					(fareakeisa71[di][0] + fareakeisa71[di][1])*0.5/m_level[ilevel].m_geom[ibox].v;\
// 			}

// #define ComptCellParams_Pts_AMR(ilevel,ibox)\
// 			m_mesh.m_level[ilevel].m_geom[ibox].boxcenter.setvalue(0.0,0.0,0.0);\
// 			for (Point_iterator p45(0,2); p45.end(); ++p45)\
// 			{\
// 				int pttt0 = m_level[ilevel].m_box[ibox].pts[p45.i][p45.j][p45.k];\
// 				m_mesh.m_level[ilevel].m_geom[ibox].boxcenter += m_mesh.m_level[ilevel].m_point[pttt0].xyz*0.125;\
// 			}\
// 			Pointxyz p2p44 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
// 											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][1]].xyz;\
// 			Pointxyz p2p46 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][1]].xyz -\
// 											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz;\
// 			Pointxyz p2p48 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
// 											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][0]].xyz;\
// 			Face aface90;\
// 			Pointxyz fareakeisa71[3][2];\
// 			for (int fd62 = 0; fd62 < 3; ++fd62)\
// 			{\
// 				for (int fi64 = 0; fi64 < 2; ++fi64)\
// 				{\
// 					Face * aimface;\
// 					int face66 = m_mesh.m_level[ilevel].m_box[ibox].faces[fd62][fi64];\
// 					if (face66 > -1)\
// 					{\
// 						aimface = &m_mesh.m_level[ilevel].m_face[face66];\
// 					}\
// 					else\
// 					{\
// 						if (fd62 == 0)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][0][0]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][1][0]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][1][1]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][0][1]].xyz);\
// 						}\
// 						else if (fd62 == 1)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][fi64][0]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][fi64][1]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][fi64][1]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][fi64][0]].xyz);\
// 						}\
// 						else if (fd62 == 2)\
// 						{\
// 							ComptFaceArea_FourPts(aface90,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][0][fi64]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][0][fi64]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][1][fi64]].xyz,\
// 																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][1][fi64]].xyz);\
// 						}\
// 						aimface = &aface90;\
// 					}\
// 					fareakeisa71[fd62][fi64] = aimface->keisa*aimface->area;\
// 				}\
// 			}\
// 			m_mesh.m_level[ilevel].m_geom[ibox].v = (abs(fareakeisa71[0][0].dot(p2p44))+\
// 																							 abs(fareakeisa71[1][1].dot(p2p46))+\
// 																							 abs(fareakeisa71[2][0].dot(p2p48)))/3.0;\
// 			for (int di = 0; di < 3; ++di)\
// 			{\
// 				m_mesh.m_level[ilevel].m_geom[ibox].keisa[di] = \
// 					(fareakeisa71[di][0] + fareakeisa71[di][1])*0.5/m_mesh.m_level[ilevel].m_geom[ibox].v;\
// 			}

#define ComptCellParams_Pts_Mesh(ilevel,ibox)\
			m_level[ilevel].m_geom[ibox].boxcenter.setvalue(0.0,0.0,0.0);\
			for (Point_iterator p45(0,2); p45.end(); ++p45)\
			{\
				int pttt0 = m_level[ilevel].m_box[ibox].pts[p45.i][p45.j][p45.k];\
				Assert(pttt0 > -1, "The box point index must be non-negative when computing the mesh geom!!!",49);\
				m_level[ilevel].m_geom[ibox].boxcenter += m_level[ilevel].m_point[pttt0].xyz*0.125;\
			}\
			m_level[ilevel].m_geom[ibox].boxcenter.get_theta();\
			Pointxyz p2p44 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][1]].xyz;\
			Pointxyz p2p46 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][1]].xyz -\
											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz;\
			Pointxyz p2p48 = m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
											 m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][0]].xyz;\
			Face aface90;\
			Pointxyz fareakeisa71[3][2];\
			for (int fd62 = 0; fd62 < 3; ++fd62)\
			{\
				for (int fi64 = 0; fi64 < 2; ++fi64)\
				{\
					if (fd62 == 0)\
					{\
						ComptFaceArea_FourPts(aface90,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][0][0]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][1][0]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][1][1]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[fi64][0][1]].xyz);\
					}\
					else if (fd62 == 1)\
					{\
						ComptFaceArea_FourPts(aface90,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][fi64][0]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][fi64][1]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][fi64][1]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][fi64][0]].xyz);\
					}\
					else if (fd62 == 2)\
					{\
						ComptFaceArea_FourPts(aface90,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][fi64]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][fi64]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][fi64]].xyz,\
															m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][1][fi64]].xyz);\
					}\
					fareakeisa71[fd62][fi64] = aface90.keisa*aface90.area;\
				}\
			}\
			m_level[ilevel].m_geom[ibox].v = (abs(fareakeisa71[0][0].dot(p2p44))+\
																				abs(fareakeisa71[1][1].dot(p2p46))+\
																				abs(fareakeisa71[2][0].dot(p2p48)))/3.0;\
			for (int di = 0; di < 3; ++di)\
			{\
				m_level[ilevel].m_geom[ibox].keisa[di] = \
					(fareakeisa71[di][0] + fareakeisa71[di][1])*0.5/m_level[ilevel].m_geom[ibox].v;\
			}

#define ComptCellParams_Pts_AMR(ilevel,ibox)\
			m_mesh.m_level[ilevel].m_geom[ibox].boxcenter.setvalue(0.0,0.0,0.0);\
			for (Point_iterator p45(0,2); p45.end(); ++p45)\
			{\
				int pttt0 = m_level[ilevel].m_box[ibox].pts[p45.i][p45.j][p45.k];\
				m_mesh.m_level[ilevel].m_geom[ibox].boxcenter += m_mesh.m_level[ilevel].m_point[pttt0].xyz*0.125;\
			}\
			m_mesh.m_level[ilevel].m_geom[ibox].boxcenter.get_theta();\
			Pointxyz p2p44 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[0][0][1]].xyz;\
			Pointxyz p2p46 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][1][1]].xyz -\
											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz;\
			Pointxyz p2p48 = m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][1]].xyz -\
											 m_mesh.m_level[ilevel].m_point[m_level[ilevel].m_box[ibox].pts[1][0][0]].xyz;\
			Face aface90;\
			Pointxyz fareakeisa71[3][2];\
			for (int fd62 = 0; fd62 < 3; ++fd62)\
			{\
				for (int fi64 = 0; fi64 < 2; ++fi64)\
				{\
					if (fd62 == 0)\
					{\
						ComptFaceArea_FourPts(aface90,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][0][0]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][1][0]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][1][1]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[fi64][0][1]].xyz);\
					}\
					else if (fd62 == 1)\
					{\
						ComptFaceArea_FourPts(aface90,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][fi64][0]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][fi64][1]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][fi64][1]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][fi64][0]].xyz);\
					}\
					else if (fd62 == 2)\
					{\
						ComptFaceArea_FourPts(aface90,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][0][fi64]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][0][fi64]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[1][1][fi64]].xyz,\
																		m_mesh.m_level[ilevel].m_point[m_mesh.m_level[ilevel].m_box[ibox].pts[0][1][fi64]].xyz);\
					}\
					fareakeisa71[fd62][fi64] = aface90.keisa*aface90.area;\
				}\
			}\
			m_mesh.m_level[ilevel].m_geom[ibox].v = (abs(fareakeisa71[0][0].dot(p2p44))+\
																							 abs(fareakeisa71[1][1].dot(p2p46))+\
																							 abs(fareakeisa71[2][0].dot(p2p48)))/3.0;\
			for (int di = 0; di < 3; ++di)\
			{\
				m_mesh.m_level[ilevel].m_geom[ibox].keisa[di] = \
					(fareakeisa71[di][0] + fareakeisa71[di][1])*0.5/m_mesh.m_level[ilevel].m_geom[ibox].v;\
			}			

#define ComptPtxyz_Refine(ilevel,mombox,newptxyz,sonx,sony,sonz,ptx,pty,ptz)\
			int spt73x = sonx+ptx;\
			int spt74y = sony+pty;\
			int spt75z = sonz+ptz;\
			int in76pt = (pt_range_e[spt73x]-pt_range_s[spt73x])*\
						 (pt_range_e[spt74y]-pt_range_s[spt74y])*\
						 (pt_range_e[spt75z]-pt_range_s[spt75z]);\
			double pt_ratio_77 = 1.0/double(in76pt);\
			for (int ni78 = pt_range_s[spt73x]; ni78 < pt_range_e[spt73x]; ++ni78)\
			{\
				for (int ni80 = pt_range_s[spt74y]; ni80 < pt_range_e[spt74y]; ++ni80)\
				{\
					for (int ni82 = pt_range_s[spt75z]; ni82 < pt_range_e[spt75z]; ++ni82)\
					{\
						int mypt84 = m_mesh.m_level[ilevel].m_box[mombox].pts[ni78][ni80][ni82];\
						newptxyz += m_mesh.m_level[ilevel].m_point[mypt84].xyz*pt_ratio_77;\
					}\
				}\
			}

#define ComptPtxyz_Refine_2d(ilevel,mombox,newptxyz,sonx,sony,sonz,ptx,pty,ptz)\
			int spt73x = sonx+ptx;\
			int spt74y = sony+pty;\
			int in76pt = (pt_range_e[spt73x]-pt_range_s[spt73x])*\
									 (pt_range_e[spt74y]-pt_range_s[spt74y]);\
			double pt_ratio_77 = 1.0/double(in76pt);\
			for (int ni78 = pt_range_s[spt73x]; ni78 < pt_range_e[spt73x]; ++ni78)\
			{\
				for (int ni80 = pt_range_s[spt74y]; ni80 < pt_range_e[spt74y]; ++ni80)\
				{\
					int mypt84 = m_mesh.m_level[ilevel].m_box[mombox].pts[ni78][ni80][ptz];\
					newptxyz += m_mesh.m_level[ilevel].m_point[mypt84].xyz*pt_ratio_77;\
				}\
			}			

#define ComptRefinePtxyz(momlevel,momindex,sonxyz,sonx,sony,sonz,px,py,pz)\
			int leftpt124 = m_mesh.m_level[momlevel].m_box[momindex].pts[0][0][0];\
			sonxyz[0] = m_mesh.m_level[momlevel].m_point[leftpt124][0]+dh[momlevel+1][0]*double(sonx+px);\
			sonxyz[1] = m_mesh.m_level[momlevel].m_point[leftpt124][1]+dh[momlevel+1][1]*double(sony+py);\
			sonxyz[2] = m_mesh.m_level[momlevel].m_point[leftpt124][2]+dh[momlevel+1][2]*double(sonz+pz)

#define ComptRefinePtxyz_2d(momlevel,momindex,sonxyz,sonx,sony,sonz,px,py,pz)\
			int leftpt124 = m_mesh.m_level[momlevel].m_box[momindex].pts[0][0][0];\
			sonxyz[0] = m_mesh.m_level[momlevel].m_point[leftpt124][0]+dh[momlevel+1][0]*double(sonx+px);\
			sonxyz[1] = m_mesh.m_level[momlevel].m_point[leftpt124][1]+dh[momlevel+1][1]*double(sony+py);\
			sonxyz[2] = m_mesh.m_level[momlevel].m_point[leftpt124][2]+dh[momlevel][2]*double(pz)	

#define ComptDerefinePtxyz(flevel,findex,cxyz,px,py,pz)\
			int fi130 = m_mesh.m_level[flevel].m_box[findex].neib[px+1][py+1][pz+1];\
			int pt131 = m_mesh.m_level[flevel].m_box[fi130].pts[px][py][pz];\
			Assert(pt131 > -1, "A negative point index when ComptDerefinePtxyz!!!", 500);\
			cxyz[0] = m_mesh.m_level[flevel].m_point[pt131][0];\
			cxyz[1] = m_mesh.m_level[flevel].m_point[pt131][1];\
			cxyz[2] = m_mesh.m_level[flevel].m_point[pt131][2]

#define ComptDerefinePtxyz_2d(flevel,findex,cxyz,px,py,pz)\
			int fi130 = m_mesh.m_level[flevel].m_box[findex].neib[px+1][py+1][1];\
			int pt131 = m_mesh.m_level[flevel].m_box[fi130].pts[px][py][pz];\
			Assert(pt131 > -1, "A negative point index when ComptDerefinePtxyz!!!", 500);\
			cxyz[0] = m_mesh.m_level[flevel].m_point[pt131][0];\
			cxyz[1] = m_mesh.m_level[flevel].m_point[pt131][1];\
			cxyz[2] = m_mesh.m_level[flevel].m_point[pt131][2]			
		
#endif
#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>
#include <cmath>
#include <mpi.h>
using namespace std;

struct Point
{

	Point ()
	{
		xy[0] = 0; xy[1] = 0; xy[2] = 0;
	}

	Point (int i, int j, int k)
	{
		xy[0] = i;
		xy[1] = j;
		xy[2] = k;
	}

	void define (int &i, int &j, int &k)
	{
		xy[0] = i;
		xy[1] = j;
		xy[2] = k;
	}

	void changex (int &i)
	{
		xy[0] = i;
	}

	Point & copy()
	{
		return *this;
	}

	Point operator+ (const Point& a_p)
	{
		Point t_p;
		t_p.xy[0] = this->xy[0] + a_p.xy[0];
		t_p.xy[1] = this->xy[1] + a_p.xy[1];
		t_p.xy[2] = this->xy[2] + a_p.xy[2];
		return t_p;
	}

	int & operator[](const int & i0)
	{
		return xy[i0];
	}

	Point operator- (const Point& a_p)
	{
		Point t_p;
		t_p.xy[0] = this->xy[0] - a_p.xy[0];
		t_p.xy[1] = this->xy[1] - a_p.xy[1];
		t_p.xy[2] = this->xy[2] - a_p.xy[2];
		return t_p;
	}

	bool operator> (const Point& a_p)
	{
		if (this->xy[0] > a_p.xy[0] &&
			this->xy[1] > a_p.xy[1] &&
			this->xy[2] > a_p.xy[2])					
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool operator< (const Point& a_p)
	{
		if (this->xy[0] < a_p.xy[0] &&
			this->xy[1] < a_p.xy[1] &&
			this->xy[2] < a_p.xy[2])
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline const int & ix() const
	{
		return xy[0];
	}
	inline const int & iy() const
	{
		return xy[1];
	}
	inline const int & iz() const
	{
		return xy[2];
	}

	bool operator==(const Point & ap)
	{
		if (xy[0] == ap.xy[0] && xy[1] == ap.xy[1] && xy[2] == ap.xy[2])
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	int xy[3];
	
};

struct Pointxyz
{
	double xyz[3];
	double pt_theta;

	Pointxyz()
	{
		xyz[0] = 0.0;
		xyz[1] = 0.0;
		xyz[2] = 0.0;
	}

	Pointxyz(const double & x0, const double & y0, const double & z0)
	{
		xyz[0] = x0;
		xyz[1] = y0;
		xyz[2] = z0;
	}

	inline void setvalue(const double & x0, const double & y0, const double & z0)
	{
		xyz[0] = x0;
		xyz[1] = y0;
		xyz[2] = z0;
	}
	inline double & operator[](const int & index)
	{
		return xyz[index];
	}

	inline void zero()
	{
		xyz[0] = 0.0;
		xyz[1] = 0.0;
		xyz[2] = 0.0;
	}

	inline void get_theta()
	{
		pt_theta = atan2(xyz[2], xyz[1]);
	}

	inline Pointxyz operator+(const Pointxyz & apt)
	{
		return Pointxyz(xyz[0]+apt.xyz[0], xyz[1]+apt.xyz[1], xyz[2]+apt.xyz[2]);
	}

	inline Pointxyz operator*(const double & a0)
	{
		return Pointxyz(xyz[0]*a0, xyz[1]*a0, xyz[2]*a0);
	}

	inline Pointxyz operator/(const double & a0)
	{
		return Pointxyz(xyz[0]/a0, xyz[1]/a0, xyz[2]/a0);
	}

	inline void operator/=(const double & a0)
	{
		xyz[0] /= a0;
		xyz[1] /= a0;
		xyz[2] /= a0;
	}

	inline void normalize()
	{
		double length = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
		//printf("The vector x %d y %d z %d\n", xyz[0], xyz[0], xyz[1]);
		xyz[0] /= length;
		xyz[1] /= length;
		xyz[2] /= length;
#ifdef DEBUG
		if (isnan(xyz[0]))
		{
			printf("NaN erro when normalize a double vector in dir x length is %f!!!\n", length);
			MPI_Abort(MPI_COMM_WORLD, 43);
		}
		if (isnan(xyz[1]))
		{
			printf("NaN erro when normalize a double vector in dir y length is %f!!!\n", length);
			MPI_Abort(MPI_COMM_WORLD, 49);
		}
		if (isnan(xyz[2]))
		{
			printf("NaN erro when normalize a double vector in dir z length is %f!!!\n", length);
			MPI_Abort(MPI_COMM_WORLD, 54);
		}
#endif		
	}

	inline double length()
	{
		return sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
	}

	inline double length_2d()
	{
		return sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	}

	inline Pointxyz operator-(const Pointxyz & p0)
	{
		Pointxyz cc;
		cc[0] = xyz[0]-p0.xyz[0];
		cc[1] = xyz[1]-p0.xyz[1];
		cc[2] = xyz[2]-p0.xyz[2];
		return cc;
	}
	// inline double operator*(const Pointxyz & p0)
	// {
	// 	Pointxyz cc;
	// 	cc[0] = this->xyz[1]*p0.xyz[2]-this->xyz[2]*p0.xyz[1];
	// 	cc[1] = this->xyz[2]*p0.xyz[0]-this->xyz[0]*p0.xyz[2];
	// 	cc[2] = this->xyz[0]*p0.xyz[1]-this->xyz[1]*p0.xyz[0];
	// 	return cc.length();
	// }

	inline Pointxyz cross(const Pointxyz & p0)
	{
		Pointxyz cc;
		cc[0] = xyz[1]*p0.xyz[2]-xyz[2]*p0.xyz[1];
		cc[1] = xyz[2]*p0.xyz[0]-xyz[0]*p0.xyz[2];
		cc[2] = xyz[0]*p0.xyz[1]-xyz[1]*p0.xyz[0];
		return cc;
	}

	inline double dot(const Pointxyz & p0)
	{
		return xyz[0]*p0.xyz[0]+xyz[1]*p0.xyz[1]+xyz[2]*p0.xyz[2];
	}


	inline double & at(const int & i0)
	{
		return xyz[i0];
	}

	inline Pointxyz & operator+=(const Pointxyz & apt)
	{
		xyz[0] += apt.xyz[0];
		xyz[1] += apt.xyz[1];
		xyz[2] += apt.xyz[2];
		return *this;
	}

	inline Pointxyz & operator*=(const double & d0)
	{
		xyz[0] *= d0;
		xyz[1] *= d0;
		xyz[2] *= d0;
		return *this;
	}

	inline void rotate_x(const double & rot_angle)
	{
		double y0 = xyz[1]*cos(rot_angle) - xyz[2]*sin(rot_angle);
		double z0 = xyz[2]*cos(rot_angle) + xyz[1]*sin(rot_angle);
		xyz[1] = y0;
		xyz[2] = z0;
	}

	inline void rotate_z(const double & rot_angle)
	{
		double x0 = xyz[0]*cos(rot_angle) - xyz[1]*sin(rot_angle);
		double y0 = xyz[1]*cos(rot_angle) + xyz[0]*sin(rot_angle);
		xyz[0] = x0;
		xyz[1] = y0;
	}
	
};

struct mPoint
{
	int index;
	Pointxyz xyz;

	inline double & operator[](const int & i0)
	{
		return xyz[i0];
	}
};




// const Point zeropt(0,0,0);
// const Point ptthree(3,3,3);

#endif
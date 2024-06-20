#ifndef _Slopelimiter_H_
#define _Slopelimiter_H_

#define minvalue(a0, a, b)\
		if (a < b)\
		{\
			a0 = a;\
		}\
		else\
		{\
			a0 = b;\
		}

#define Superbee(qr, lvalue)\
		lvalue = max(min(2.0*qr, 1.0), min(qr, 2.0))

#define VanLeer(qr, lvalue)\
		lvalue = 2.0*qr/(1+qr)	

#define ALimiter(qr,limitvalue) VanLeer(qr,limitvalue)

#define SlopeCorrection_Left(q_i, qi, qi1)\
		dqi = qi - q_i;\
		dqii = qi1 - qi;\
		qr = dqii/dqi;\
		if (abs(qr) > 1.0e10 || isnan(qr))\
		{\
			if (dqii > 0.0) qr = 1000000.0;\
			else qr = -1.0;\
		}\
		if (qr <= 0.0)\
		{\
			qr = 0.0;\
		}\
		ALimiter(qr, limitvalue);\
		Assert(!isnan(limitvalue), "Limit value left nan error!!!", 5);\
		Assert(limitvalue > -0.0000001 && limitvalue < 2.0000001, "Limit value left range error!!!", 5)

#define SlopeCorrection_Right(q_i, qi, qi1)\
		dqi = qi - q_i;\
		dqii = qi1 - qi;\
		qr = dqi/dqii;\
		if (abs(qr) > 1.0e10 || isnan(qr))\
		{\
			if (dqi > 0.0) qr = 1000000.0;\
			else qr = -1.0;\
		}\
		if (qr <= 0.0)\
		{\
			qr = 0.0;\
		}\
		ALimiter(qr, limitvalue);\
		Assert(!isnan(limitvalue), "Limit value right nan error!!!", 5);\
		if (limitvalue < -0.0000001 || limitvalue > 2.0000001)\
		{\
			printf("Limit value right range error!!! Right limitvalue %f dqi %f dqii %f qr %f\n", limitvalue, dqi, dqii, qr);\
			MPI_Abort(MPI_COMM_WORLD, 57);\
		}

class Slopelimiter
{
public:	
	Slopelimiter()
	{}

	inline void ApplyLimiter(const double & qr, double & phivalue)
	{
		minmod(qr, phivalue);
#ifdef DEBUG		
		if (isnan(phivalue))
		{
			printf("Limiter give a NaN!!!\n");
			MPI_Abort(MPI_COMM_WORLD, 560);
		}
#endif
	}

	inline void doubleminmod(const double & qr, double & phivalue)
	{		
		minvalue(phivalue, 2.0*qr, 2.0);
		minvalue(phivalue, (1.0+qr)/2.0, phivalue);
	}

	inline void minmod(const double & qr, double & phivalue)
	{
		minvalue(phivalue, qr, 1.0);
	}

	inline void superbee(const double & qr, double & phivalue)
	{
		phivalue = max(min(2.0*qr, 1.0), min(qr, 2.0));
	}

};

#endif
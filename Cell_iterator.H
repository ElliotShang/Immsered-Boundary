#ifndef _Cell_iterator_H_
#define _Cell_iterator_H_

class Cell_iterator
{
private:
	Mesh & m0;
	int preloc;
	const int & level0;
	const int & cbloc;
	const int & pi;
	const int & pj;
	const int & pk;
	int is, js, ks;
	static int nbse[3][2];	

public:
	Cell_iterator(Mesh & amesh, const int & alevel0, const int & acbloc, const int & i0, const int & j0, const int & k0):
	m0(amesh),level0(alevel0),cbloc(acbloc),pi(i0),pj(j0),pk(k0)
	{
		is = nbse[i0][0];
		js = nbse[j0][0];
		ks = nbse[k0][0];
	}

	// Cell_iterator(Mesh & amesh, const int & alevel0, const int & acbloc, const int & i0, const int & j0, const int & k0):
	// m0(amesh),cbloc(*abox),pi(i0),pj(j0),pk(k0)
	// {
	// 	is = nbse[i0][0];
	// 	js = nbse[j0][0];
	// 	ks = nbse[k0][0];
	// }

	inline Cell_iterator & operator++()
	{
		++ks;
		if (ks == nbse[pk][1])
		{
			ks = nbse[pk][0];
			++js;
			if (js == nbse[pj][1])
			{
				js = nbse[pj][0];
				++is;
			}
		}
		return *this;
	}

	inline bool end()
	{
		if (is == nbse[pi][1])
		{
			return false;
		}
		else
		{
			preloc = m0.BoxNeib(level0, cbloc, is, js, ks);
			return true;
		}
	}

	inline const int & index()
	{
		return preloc;
	}

};

#endif
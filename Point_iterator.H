#ifndef _Point_iterator_H_
#define _Point_iterator_H_

class Point_iterator
{
public:
	Point_iterator()
	{}

	Point_iterator(int start0, int end0):
	start(start0), endside(end0)
	{
		i = start0;
		j = start0;
		k = start0;
	}

	Point_iterator(const int & i0, const int & j0, const int & k0)
	{
		i = i0;
		j = j0;
		k = k0;
	}

	int i;
	int j;
	int k;

	inline Point_iterator & operator ++()
	{
		++k;
		if (k == endside)
		{
			k = start;
			++j;
			if (j == endside)
			{
				j = start;
				++i;
			}
		}
		return *this;
	}

	bool inrange(int low, int high)
	{
		if (i >= low && i <= high &&
			j >= low && j <= high &&
			k >= low && k <= high)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool operator!=(const Point_iterator & ap)
	{
		if (i != ap.i || j != ap.j || k != ap.k)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool operator < (int a)
	{
		if (i < a && j < a && k < a)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool operator > (int a)
	{
		if (i > a && j > a && k > a)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool end()
	{
		if (i == endside)			
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline int & operator[](const int & ai)
	{
		return *(&i+ai);
	}

private:
	int start;
	int endside;
};

class Four_iterator
{
public:
	Four_iterator();

	Four_iterator(int start0, int end0):
	start(start0), endside(end0), p3(start0, end0)
	{
		i = start0;
		j = start0;
		k = start0;
		s = start0;
	}

	int i, j, k, s;

	Four_iterator & operator ++()
	{
		s += 1;
		if (s == endside)
		{
			++p3;
			s = start;
			i = p3.i;
			j = p3.j;
			k = p3.k; 
		}
		return *this;
	}

	bool end()
	{
		if (i == endside)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

private:
	int start, endside;
	Point_iterator p3;
	
};

class Point_iterator_2d
{
public:
	Point_iterator_2d()
	{}

	Point_iterator_2d(int start0, int end0):
	start(start0), endside(end0)
	{
		i = start0;
		j = start0;
		k = 0;
	}

	int i;
	int j;
	int k;

	inline Point_iterator_2d & operator ++()
	{
		++j;
		if (j == endside)
		{
			j = start;
			++i;
		}
		return *this;
	}

	inline bool end()
	{
		if (i == endside)			
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline int & operator[](const int & ai)
	{
		return *(&i+ai);
	}

private:
	int start;
	int endside;
};
#endif
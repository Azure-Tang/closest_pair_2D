// Mohamad Al Mustapha
// Implementation of Divide-and-Conquer algorithm to find closest pair of points in XY plane

#ifndef CLOSEST_PAIR_H
#define CLOSEST_PAIR_H

#include <vector>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <fstream>
using namespace std;

class Point
{
public:
	friend std::ostream& operator<<(std::ostream& os, const Point& pt)
	{
		os << "(" << pt.x << "," << pt.y << ")";
		return os;
	}

	inline bool byX(const Point& p)
	{
		return x < p.x;
	}

	inline bool byY(const Point& p)
	{
		return y < p.y;
	}

	// coordinates
	double x;
	double y;
	// original index of the point
	size_t index;
};

class Pair
{
public:
	Pair() :index1(NULL), index2(NULL), distance(DBL_MAX){};

	size_t index1;
	size_t index2;
	double distance;

	friend bool operator< (const Pair & a, const Pair & b)
	{
		if (a.distance < b.distance)
			return true;
		return false;
	}
};

class ClosestPair
{
public:
	ClosestPair(){};

	ClosestPair& read_file(char* filepath)
	{
		ifstream file(filepath);
		if (!file)
		{
			cout << "Error opening file" << endl;
			return *this;
		}

		double x, y;
		Point p;

		//read abscissa of first point to make sure the file is not empty
		file >> x;
		for (size_t i = 0; !file.eof(); ++i)
		{
			//Format: XX.XX XX.XX
			file >> y;

			p.x = x;
			p.y = y;
			p.index = i;
			points.push_back(p);
			file >> x;
		}

		file.close();
		return *this;
	}

	//a function to generate points with random coordinates within a rectangle 
	void generate_random_points(double bottom_x, double top_x, double bottom_y, double top_y, size_t num)
	{
		srand(time(0));

		points.clear();
		points = vector<Point>(num);

		for (size_t i = 0; i < num; ++i)
		{
			points[i].x = random_double(bottom_x, top_x);
			points[i].y = random_double(bottom_y, top_y);
			points[i].index = i;
			cout << points[i] << endl;
		}
	}

	// a function which finds the closest pair of points in theta n squared
	Pair brute_force_closest_pair(double &time)
	{
		if (points.size() < 2)
			return Pair();

		long double sum;
		clock_t start_s = clock();

		Pair p = _brute_force_closest_pair(points);
		cout << "Minimum Distance between " << points[p.index1] << " and " << points[p.index2] << " at " << sqrt(p.distance) << endl;

		clock_t stop_s = clock();
		sum = ((long double)(stop_s - start_s) / CLOCKS_PER_SEC) / 100000;
		time = sum;
		return p;
	}

	Pair DAC_closest_pair(double& time)
	{
		if (points.size() < 2)
			return Pair();

		long double sum;
		clock_t start_s = clock();

		vector<Point> p1(points.begin(), points.end());
		vector<Point> p2(points.begin(), points.end());

		//sort p1 by increasing x coordinate
		sort(p1, &Point::byX);
		//sort p2 by increasing y coordinate
		sort(p2, &Point::byY);

		Pair p = _recursive_closest_pair(p1, p2);
		cout << "Minimum Distance between " << points[p.index1] << " and " << points[p.index2] << " at " << sqrt(p.distance) << endl;
		clock_t stop_s = clock();
		sum = ((double)(stop_s - start_s) / CLOCKS_PER_SEC) / 100000;
		time = sum;
		return p;
	}

	size_t size()
	{
		return points.size();
	}

	~ClosestPair()
	{
	}

private:
	vector<Point> points;

	Pair _recursive_closest_pair(vector<Point>& X, vector<Point>& Y)
	{
		if (X.size() <= 3)
		{
			return _brute_force_closest_pair(X);
		}
		else
		{
			size_t size_R = X.size() / 2;
			size_t size_L = X.size() - size_R;
			double L = X[size_L - 1].x;

			//split the arrays
			vector<Point> XL(size_L), XR(size_R), YL, YR;

			vector<Point>::iterator it = X.begin() + size_L;
			XL = vector<Point>(X.begin(), it);
			XR = vector<Point>(it, X.end());

			size_t size_Y = Y.size();

			for (size_t i = 0; i < size_Y; ++i)
			{
				// if Y[i] belongs to PL
				if (Y[i].x <= L)
					YL.push_back(Y[i]);
				else
					YR.push_back(Y[i]);
			}

			//call the function recursively on the left and right arrays
			Pair GL = _recursive_closest_pair(XL, YL);
			Pair GR = _recursive_closest_pair(XR, YR);

			//find the minimum between both
			Pair G;
			if (GL < GR)
				G = GL;
			else G = GR;

			// find the points in the zone around the median
			vector<Point> Y_prime;
			for (size_t i = 0; i < size_Y; ++i)
			{
				double dist = Y[i].x - L;
				if (dist*dist <= G.distance)
					Y_prime.push_back(Y[i]);
			}

			size_t size_Y_prime = Y_prime.size();
			Pair G_prime;

			if (size_Y_prime < 2)
				return G;

			for (size_t i = 0; i < size_Y_prime - 1; ++i)
			{
				for (size_t j = i + 1; j < i + 7 && j < size_Y_prime; ++j)
				{
					double dist = distance_squared(Y_prime[i], Y_prime[j]);
					if (dist < G_prime.distance)
					{
						G_prime.distance = dist;
						G_prime.index1 = Y_prime[i].index;
						G_prime.index2 = Y_prime[j].index;
					}
				}
			}

			if (G_prime < G)
				return G_prime;
			return G;
		}
	}

	inline Pair _brute_force_closest_pair(const vector<Point>& pts)
	{
		size_t size = pts.size();
		Pair pair;

		if (size <= 1)
			return pair;

		for (size_t i = 0; i < size - 1; ++i)
		{
			for (size_t j = i + 1; j < size; ++j)
			{
				double dist = distance_squared(pts[i], pts[j]);
				if (dist < pair.distance)
				{
					pair.distance = dist;
					pair.index1 = pts[i].index;
					pair.index2 = pts[j].index;
				}
			}
		}

		return pair;
	}

	// a function to generate a random double number between bottom and top
	double random_double(double bottom, double top)
	{
		double r = double(rand()) / RAND_MAX;
		return bottom + r * (top - bottom);
	}

	//a function to compute the squared distance between two points
	double distance_squared(const Point& p1, const Point& p2)
	{
		double x = p1.x - p2.x;
		double y = p1.y - p2.y;
		return (x*x + y*y);
	}

	// sorting using binary max heap with complexity: O(n log n)
	// adapted from Dr. Florin Balasa
	void sort(vector<Point>& A, bool (Point::*fptr)(const Point&))
	{
		size_t crtSize = A.size();
		// build a binary max heap from A
		for (int i = crtSize / 2 - 1; i >= 0; i--)
		{
			Point x = A[i];
			size_t hole = i;
			size_t child;
			// percolate the hole down
			while ((child = hole * 2 + 1) < crtSize)
			{
				if (child + 1 < crtSize && (A[child].*fptr)(A[child + 1]))
					child++;
				if ((x.*fptr)(A[child]))
				{
					A[hole] = A[child];
					hole = child;
				}
				else break;
			}
			A[hole] = x;
		}

		// delete the max heap element A[0] from the heap
		while (crtSize > 1)
		{
			Point x = A[--crtSize];
			A[crtSize] = A[0];
			size_t hole = 0;
			size_t child;
			// re-adjust the max heap
			while ((child = hole * 2 + 1) < crtSize)
			{
				if (child + 1 < crtSize && (A[child].*fptr)(A[child + 1]))
					child++;
				if ((x.*fptr)(A[child]))
				{
					A[hole] = A[child];
					hole = child;
				}
				else break;
			}
			A[hole] = x;
		}
	}
};

#endif
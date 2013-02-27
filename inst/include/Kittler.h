/**********************************************************************************
* TITLE:        Kittler.hpp
*
* PURPOSE:      PERFORMS KITTLER'S AUTOMATIC THRESHOLDING ALGORITH .
*
* AUTHOR:       Rave Harpaz
*               Pattern Recognition Laboratory
*               Department of Computer Science
*               The Graduate Center
*               The City University of New York
*               365 Fifth Avenue, New York, New York 10016
*               e-mail: rave_harpaz@yahoo.com
*
* DATE:         03/25/2004
*
*
* VERSION:      1.00
*
* LANGUAGE:     C++
*
* SYSTEM:       Red Hat Linux 9
*
* COMPILER:     GNU g++
*
* REFERENCES:   J. Kittler & J. Illingworth: "Minimum Error Thresholding"
                Pattern Recognition, Vol 19, nr 1. 1986, pp. 41-47.
*
* REVISIONS:
*
* Copyright 2002 Pattern Recognition Laboratory, The City University of New York
*********************************************************************************/
#ifndef KITTLER_H
#define KITTLER_H

#include <vector>
#include <deque>

using namespace std;

static const double SMALL_DEPTH = 1.0e-5;
static const double SMALL_NUMBER = 1.0e-5;

class Kittler
{
public:
	Kittler():Depth(0), Discriminability(0), GlobalMin(0), Threshold(0) {};
	bool FindThreshold(const vector<double>& H, double RHmin, double RHmax);
	double GetDepth() {
		return Depth;
	}
	double GetDiscrim() {
		return Discriminability;
	}
	double GetThreshold() {
		return Threshold;
	}
private:
	deque<bool> MarkMinima( const vector<double> &J);
	bool FindGlobalMin(const deque<bool> &IsMin, const vector<double> &J);
	
	// separation statistics
	double Depth;                                           // largest deifference in criterion function
	double Discriminability;
	double GlobalMin;
	double Threshold;

};


#endif

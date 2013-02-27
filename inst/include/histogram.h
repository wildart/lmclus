#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <iostream> 

using namespace std;


/* class Histogram
 * ---------------
 * used to create 2D frequency histograms with constant bin width.
 * bin width is either computed using a heuristic or by specifying the number of bins desired.
 */
class Histogram{
	public:
  		Histogram(){};
  		void createHistogram ( const vector<double> &data, int NumBins );
		void createHistogram ( const vector<double> &data, double portion );
  		~Histogram(){};
  		void insert_data(vector<double> &data);                    // insert all data into histogram
  		int get_element_bin(double data)const;                     // return the bin number associated with the data
  		double get_min_H () const {return Min_H;}
  		double get_max_H () const {return Max_H;}
  		vector<double> get_H() const {return H_n;}
  		void display(const double thres );
  		void ShowBins(const double thres);
  		void Normalize();                                          // convert histogram to relative frequancy histogram
  		friend ostream & operator<<(ostream &out, const Histogram &h);
	private:
  		vector<double> H_n;           // the histogram, points per bin
  		unsigned int Size_pt;                  // number of points in histogram
  		int Max_Bin_Portion;          // max number of points per bin
  		double Min_pt;                // smallest data point 
  		double Max_pt;                // largest data point
  		double Min_H;                 // left bound in histogram
  		double Max_H;                 // right bound in histogram
  		double Range_pt;              // range of points
  		double Range_H;               // range of histogram
  		double Bin_Width;             
  		unsigned int Num_Of_Bins;       
};

# endif

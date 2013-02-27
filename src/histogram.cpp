#include <cmath>
#include <iterator>
#include <algorithm>
#include "histogram.h"

/* createHistogram
 * ---------------
 * Generate an empty histogram setting the parameters for the histogram's structure based
 * on the data to be stored in the histogram.
 * The bin-width of the histogram is determined by the formula: bin-Width=min{ Z_n+i - Z_i },i.e.
 * find the smallest difference between i successive points.
 * (where Z_n is a point with index n, and i is the max number of points we allow
 * to be stored in a single bin).
 */
void Histogram::createHistogram ( const vector<double> &data, double portion )
{
  	vector<double> Sorted=data;
  	sort(Sorted.begin(),Sorted.end());                      // sort the data

  	Size_pt=Sorted.size();
  	Max_Bin_Portion=static_cast<int>(Size_pt * portion);    // max number of points per bin
  	Min_pt=Sorted[0];                                       // smallest point from  the sorted data 
  	Min_H=Min_pt;                                           // left bound of histogram
  
  	Max_pt=Sorted[Size_pt-1];                               // largest point from the sorted data 
	Range_pt=Max_pt - Min_pt;                               // range of data
 
	//cout << "init hist: "<< Size_pt << ",  " << Max_Bin_Portion << endl;
 
  	// find minimum width using the above formula
  	Bin_Width=Range_pt + 1;                                 // upper bound for bin width
	if (Max_Bin_Portion > 0)
  	for (int i=Size_pt-1; i-Max_Bin_Portion+1 >= 0; i--)
  	{ 
		int j = i-Max_Bin_Portion+1;
		double width=Sorted[i]-Sorted[j];   // distance between 'MaxBinPortion' of consecutive points
		if (width<Bin_Width && width!=0)
			Bin_Width=width;
  	}
	
	//cout << "find minimum" << endl;
  
  	// determine the number of bins in the histogram, based on the range of
  	// the histogram and the bin-width.
  	Num_Of_Bins=static_cast<int> (Range_pt/Bin_Width)+1;         
  	Max_H=Min_H + (Bin_Width * Num_Of_Bins);
  	Range_H=Max_H - Min_H; 

  	
  	Bin_Width=Range_pt/static_cast<double>(Num_Of_Bins-1);

	// deal with special cases
  	if(Num_Of_Bins<=2)
  	{
 		//cerr<<"n="<<Num_Of_Bins<<" R="<<Range_pt<<" bw="<<Bin_Width<<endl;
		Num_Of_Bins=30;
    		Bin_Width=Range_pt/static_cast<double>(Num_Of_Bins-1);
		Max_H=Min_H + (Bin_Width * Num_Of_Bins);
    		Range_H=Max_H - Min_H;
		//exit(1);
  	}
	//cout << "deal with special cases" << endl;
  	
	if(Num_Of_Bins > static_cast<int>(ceil(static_cast<double>(Size_pt)/20)) )
  	{
		//cerr<<"correcting bin size, range="<<Range_pt<<" binwidth="<<Bin_Width<<"  number of bins="<<Num_Of_Bins<<endl;
    	Num_Of_Bins=static_cast<int>(ceil(static_cast<double>(Size_pt)/20));
		if(Num_Of_Bins<2)
			Num_Of_Bins=10;
    	Bin_Width=Range_pt/static_cast<double>(Num_Of_Bins-1);
    	Max_H=Min_H + (Bin_Width * Num_Of_Bins);
    	Range_H=Max_H - Min_H;
    	//cerr<<" corrected binwidth="<<Bin_Width<<" Num_Of_Bins="<<Num_Of_Bins<<endl;
  	}	
  
	//cerr<<"n="<<Num_Of_Bins<<" R="<<Range_pt<<" bw="<<Bin_Width<<endl;
  	vector<double> His(Num_Of_Bins,0);                     // initialize all bin counts to zero
  	H_n=His;
}


/* createHistogram
 * ---------------
 * Generate an empty constant size histogram.
 */
void Histogram::createHistogram ( const vector<double> &data, int NumBins )
{
  	vector<double> Sorted=data;
  	sort(Sorted.begin(),Sorted.end());                      // sort the data

  	Size_pt=Sorted.size();
  	Min_pt=Sorted[0];                                       // smallest point from  the sorted data 
  	Min_H=Min_pt;   
  	Max_pt=Sorted[Size_pt-1];                               // largest point from the sorted data 
  	Range_pt=Max_pt - Min_pt;                               // range of data
 
    
  	Num_Of_Bins=NumBins;
  	Bin_Width=Range_pt/static_cast<double>(Num_Of_Bins-1);
  	Max_H=Min_H + (Bin_Width * Num_Of_Bins);
  	Range_H=Max_H - Min_H; 

  
  	vector<double> His(Num_Of_Bins,0);                     // initialize all bin count to zero
  	H_n=His;

	return;
}



/* get_element_bin
 * ---------------
 * return the bin number associated with a data point.
 */
int Histogram::get_element_bin(double data) const
{
  	// associate point with bin by using the formula: BinIndex=(D_n - left_bound)/BinWidth.
  	// this gives the index of the bin that a falls into. 
  	int BinIndex;
  	BinIndex=static_cast<int>( (data - Min_pt) / Bin_Width );
  
  	return BinIndex;
}


/* insert_data
 * -----------
 * insert entire data into histogram.
 */
void Histogram::insert_data(vector<double> &data)
{
  	// associate point with bin by using the formula: BinIndex=(D_n - left_bound)/BinWidth.
  	// this gives the index of the bin that a falls into.  
  	unsigned int BinIndex;
  	
  	for (unsigned int i=0; i< data.size(); i++)
  	{    
		BinIndex=static_cast<int>( (data[i] - Min_pt) / Bin_Width );
     	if (BinIndex<0 || BinIndex >=Num_Of_Bins)
		{ 
			cout<<"histogram bin index out of range"<<endl;
     			cout<<"data="<<data[i]<<" min pt="<<Min_pt<<" bw="<<Bin_Width; 
     			cout<<" bin index "<<BinIndex<<endl; 
     			display(0); 
     			//exit(1);
		}
		else
			H_n[BinIndex]++;                                           // increment bin counter
  	}
  	return;
}	



/* ShowBins
 * ----------
 * display bin counts for each bin in the histogram and highlight(color escape sequnce C-q+[)
 * a bin associated with a point(values are doubles).
 */
void Histogram::ShowBins(const double thres)
{    	
  	unsigned int index=get_element_bin(thres);    
  
  	for(unsigned int i=0; i < H_n.size(); i++)
  	{
    		if (i==index  )
      			cout<<"[31m"<<H_n[i]<<"[0m ";
    		else
      			cout<<H_n[i]<<" ";
       	}
  	cout<<endl;

	return;
}

/* Normalize
 * ----------
 * normalize histogram to a relative frequency histogram. Used by kittler's threshold algorithm.
 */
void Histogram::Normalize()
{
	for(unsigned int i=0;i<Num_Of_Bins;i++)
		H_n[i]=H_n[i] / static_cast<double>(Size_pt);

	return;
}





/* display
 * -------
 * used to display histogram debugging info with color highlighting.
 */
void Histogram::display(const double thres)
{
  	cout<<"data size="<<Size_pt<<endl;
  	cout<<"max number of points per bin="<<Max_Bin_Portion<<endl;
  	cout<<"bin Width="<<Bin_Width<< endl;
  	cout<<"number of bins="<<Num_Of_Bins<<endl;
  	cout<<"range of data="<<Min_pt<<"<-->"<<Max_pt<<" ( "<<Range_pt<<" )"<<endl;
  	cout<<"range of histogram="<<Min_H<<"<-->"<<Max_H<<" ( "<<Range_H<<" )"<<endl;
  	cout<<"histogram:"<<endl;
  	ShowBins(thres);
  	cout<<endl;

	return;
}

/* operator<<
 * ----------
 * overload the insertion operator to enable output of histogram.
 * (used by the description class to output model/description info)
 */
ostream & operator<<(ostream &out, const Histogram &h)
{
  	out<<h.Size_pt<<endl;
  	out<<h.Max_Bin_Portion<<endl;
  	out<<h.Bin_Width<< endl;
  	out<<h.Num_Of_Bins<<endl;
  	out<<h.Min_pt<<endl;
  	out<<h.Max_pt<<endl;
  	out<<h.Min_H<<endl;
  	out<<h.Max_H<<endl;
  
  	for(unsigned int i=0; i < h.H_n.size(); i++)
    		out<<h.H_n[i]<<" ";
  	out<<endl;

  	return out;
}

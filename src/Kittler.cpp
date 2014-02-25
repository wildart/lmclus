/* Copyright 2005-2014 Pattern Recognition Laboratory, The City University of New York
 * =====================================================================================
 *
 *       Filename:  Kittler.cpp
 *
 *    Description:  Kittler's minimum error thresholding algorithm
 *
 *        Version:  1.0
 *       Compiler:  gcc 4.7.2
 *
 *         Author:  Rave Harpaz, Art Wild (wildart@gmail.com) 
 *        Company:  Patter Recognition Lab, GC CUNY
 *
 * =====================================================================================
 */
#include <cmath>
#include <deque>
#include <vector>
#include "Kittler.hpp"

/* FindThreshold
 * -------------
 * using kittler's minimum error thresholding algorithm, find the threshold, and separation statistics of a histogram.
 * input: H-histogram in frequancy or relative frequancy format ( use relative frequancy for normalized separation statistics).
 *        RHmin, RHmax - range of values in histogram.
 * output: false - histogram is unimodal and cannot be thresholded, true - histrogram is bimodal can be thresholded
 */
bool Kittler::FindThreshold(const arma::vec &H, double RHmin, double RHmax) {
    int N = H.n_elem;
    vector<double> P1(N, 0);   // P1[i]- prior to the left of threshold i
    vector<double> P2(N, 0);   // P2[i]- prior to the right of threshold i
    vector<double> Mu1(N, 0);
    vector<double> Mu2(N, 0);
    vector<double> Var1(N, 0);
    vector<double> Var2(N, 0);

    // stopping cases , threshold i threshold separtes bin i from the rest
    P1[0] = H[0];
    P2[N-2] = H[N-1];
    Mu1[0] = 0;
    Mu2[N-2] = (H[N-1] == 0 ? 0 : (N-1));
    Var1[0] = 0;
    Var2[N-2] = 0;

    // recursive defintions
    for (int i = 1, j = N-3; i <= N-2 ; i++, j--) {
        P1[i] = P1[i-1] + H[i];
        if ( P1[i] != 0 ) {
            Mu1[i] = ((Mu1[i-1] * P1[i-1]) + (i * H[i])) / P1[i];
            Var1[i]= (P1[i-1] *
                        (Var1[i-1] + (Mu1[i-1]-Mu1[i]) * (Mu1[i-1]-Mu1[i])) +
                        H[i] * (i - Mu1[i]) * (i - Mu1[i]) ) / P1[i];
        } else {
            Mu1[i] = 0, Var1[i] = 0;
        }

        P2[j] = P2[j+1] + H[j+1];
        if ( P2[j] != 0 ) {
            Mu2[j] = ((Mu2[j+1] * P2[j+1]) + ((j+1) * H[j+1])) / P2[j];
            Var2[j]= (P2[j+1] *
                        (Var2[j+1] + (Mu2[j+1]-Mu2[j]) * (Mu2[j+1]-Mu2[j])) +
                        H[j+1] * (j+1 - Mu2[j]) * (j+1 - Mu2[j]) ) / P2[j];
        } else {
            Mu2[j] = 0, Var2[j] = 0;
        }
    }

    // compute criterion function
    vector<double> J(N-1, 0);

    for (int T = 0; T < N-1 ; T++)
        // if (P1[T] && P2[T] && Var1[T]>SMALL_NUMBER && Var2[T]>SMALL_NUMBER)
            J[T] = 1 + 2*(P1[T]*log(sqrt(Var1[T])) + P2[T]*log(sqrt(Var2[T])))
                    - 2*(P1[T]*log(P1[T]) + P2[T]*log(P2[T]));
        // else
        //    J[T] = -100000000;               // undefined

    CriterionFunc = J;

    // mark all local minima of the criterion function J
    deque<bool> IsMin = MarkMinima(J);

    if (FindGlobalMin(IsMin, J) == false) {  // no minimum, unimode histogram
        return false;
    } else {  // compute separation statistics
        MinIndex = static_cast<int> (GlobalMin);
        Threshold = RHmin + ((GlobalMin+1) * (RHmax - RHmin) / N);
        Discriminability = (fabs(Mu1[MinIndex]-Mu2[MinIndex]))/
                            (sqrt(Var1[MinIndex]+Var2[MinIndex]));
        return true;
    }
}

/* MarkMinima
 * ----------
 * Mark local minima of criterion function values.
 * input: J-criterion function values for each bin.
 * output: positions of local minima for values for criterion fucntion.
 */
deque<bool> Kittler::MarkMinima(const vector<double> &J) {
    int N = J.size();
    deque<bool> IsMin(N, false);

    if (N-1 < 1) return IsMin;

    int MinCount = 0;
    double DifPrev = J[1]-J[0], DifCur = 0;

    for (int i = 1; i < N-1; i++) {
        DifCur = J[i+1]-J[i];

        if (DifPrev <= 0 && DifCur >= 0) {
            IsMin[i] = true;
            ++MinCount;
        }
        DifPrev = DifCur;
    }

    return IsMin;
}

/* FindGlobalMin
 * -------------
 * find global minima of criterion funtion if one exists.
 * input: J-criterion function values for each bin.
 *        IsMin-positions of local minima for values for criterion fucntion.
 * output: does a a global minima exits for values for criterion fucntion.
 */
bool Kittler::FindGlobalMin(const deque<bool> &IsMin, const vector<double> &J) {
    // find first minimum
    unsigned int LeftMin = 0;
    while (LeftMin < IsMin.size() && IsMin[LeftMin] == false) LeftMin++;

    if (LeftMin == IsMin.size())
        return false;                       // no minimum, unimode histogram

    while (LeftMin < IsMin.size()) {
        // Dedect flat
        unsigned int RightMin = LeftMin;
        double LocalMin;
        while (IsMin[RightMin] == true && RightMin < IsMin.size()) RightMin++;
        LocalMin = static_cast<double>(LeftMin + RightMin-1)/2.0;

        // Monotonically ascend to the left
        unsigned int LeftHeight = static_cast<int>(LocalMin);
        while (LeftHeight >=1 && J[LeftHeight-1] >= J[LeftHeight])
            LeftHeight--;

        // Monotonically ascend to the right
        unsigned int RightHeight = static_cast<int>(LocalMin);
        while (RightHeight < IsMin.size()-1 &&
                J[RightHeight] <= J[RightHeight+1])
                RightHeight++;

        // choose minimum of height values and compute depth
        double LocalDepth = 0;
        if (J[LeftHeight] < J[RightHeight])
            LocalDepth = J[LeftHeight]-J[static_cast<int>(LocalMin)];
        else
            LocalDepth = J[RightHeight]-J[static_cast<int>(LocalMin)];

        if (LocalDepth > Depth) {
            Depth = LocalDepth;
            GlobalMin = LocalMin;
        }

        // get next minima
        LeftMin = RightMin;
        while (LeftMin < IsMin.size() && IsMin[LeftMin] == false) LeftMin++;
    }


    if (Depth < SMALL_DEPTH)
        return false;   // no minimum, unimode histogram
    else
        return true;
}

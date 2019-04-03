#pragma once

#include <map>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstddef>             // for size_t
#include <cmath>               // for fabs

#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>


namespace beamlinereco {
    template <class T>
    struct hit_t {
        unsigned int  DigitizerChannel;
        uint32_t      Timestamp;
        T             TStartInNanoSec;
        T             TPeakInNanoSec;
        T             IntegratedChargeInADCTimeTicks;
        T             IntegratedChargeInADCNanoSec;
        T             AmplitudeInMiliVolt;
        T             AmplitudeInADC;
        T             RiseTimeInNanoSec;
        T             FallTimeInNanoSec;
        bool          IsContained;
    };

    enum CFDParams {
        kADCNBits,
        kADCDynamicRange,
        kADCOffset,
        kTimeSamplingInterval,
        kNSamplingPoints,
        kIsWaveformNegativePolarity,
        kCFDThreshold,
        kRawHitFinderThresholdInNoiseSigma,
        kShortRawHitIgnoringDurationInTicks,
        kConsecutiveHitSeperationDurationInTicks,
        kGSFilter,
        kGSFilterWindow,
        kGSFilterDegree,
        kIntergratedWindowFixed,
        kIntergratedWindowLowerLimitIndex,
        kIntergratedWindowUpperLimitIndex,
    };

    std::vector<double> sg_smooth(const std::vector<double> &v, const int w, const int deg);
    std::vector<double> sg_derivative(const std::vector<double> &v, const int w, const int deg, const double h = 1.0);

    class SGSmoothing {
    public:
        static void Smooth(size_t nsamples, double* in, double* out, const int w, const int deg);
        static void Derivative(size_t nsamples, double* in, double* out, const int w, const int deg, const double h = 1.0);
    };

    template <class T>
    class CFDHitFinder {
    public:
        CFDHitFinder() {
            _HitCollection.clear();
            _RawHitLogicNanosec.clear();
            _PedestalInADC = 0;
            _NoiseSigmaInADC = 0;
            _CFDThresholdInADC = 0;
            _WaveformADCNanosec.clear();
            _ChannelNumber = 993;
            _Timestamp = 0;
            _NonfilterWaveformADCNanosec.clear();
            _CFDParamSet.clear();
        }
        virtual ~CFDHitFinder() = default;

        
        void SetWaveform(std::vector<uint16_t>& waveform, unsigned int channelNo, uint32_t timestamp);
        void SetParams(const std::map<CFDParams, double>& paramSet);

        inline const std::map<T, hit_t<T> >& GetHitCollection() const { return _HitCollection; }
        inline const T& GetPedestal() const { return _PedestalInADC; }
        inline const T& GetNoiseSigma() const { return _NoiseSigmaInADC; }
        inline const T& GetCFDThreshold() const { return _CFDThresholdInADC; }

        virtual void Go();

    private:
        
        void SetChannel(unsigned int channelNo);
        void SetTimestamp(uint32_t timestamp);
        void SetParam(CFDParams param, T value);

    private:
        unsigned int   _ChannelNumber;
        uint32_t       _Timestamp;
        std::vector<T> _WaveformADCNanosec;
        std::vector<T> _NonfilterWaveformADCNanosec;
        std::vector<bool> _RawHitLogicNanosec;
        T _PedestalInADC;
        T _NoiseSigmaInADC;
        T _CFDThresholdInADC;

        std::map<CFDParams, T> _CFDParamSet;
        std::map<T, hit_t<T> > _HitCollection;

    private:
        virtual void FindPedestal();
        virtual void FindRawHitLogic();
        virtual void FindCFDHits();
        virtual bool BackwardFindingOfHitStart(size_t hitPeakTimeIndex, T hitPeakValue, T& hitStartTimeIndex, T& hitRiseTimeInIndex);
        virtual bool ForwardFindingOfHitFallTime(size_t hitPeakTimeIndex, T& hitFallTimeInIndex);
        virtual T    IntegrateWaveformInADC(size_t hitStartTimeIndex, size_t hitEndTimeIndex);

        virtual inline void Reset() {
            _HitCollection.clear();
            _RawHitLogicNanosec.clear();
            _PedestalInADC = 0;
            _NoiseSigmaInADC = 0;
            _CFDThresholdInADC = 0;
            _WaveformADCNanosec.clear();
            _ChannelNumber = 993;
            _Timestamp = 0;
            _NonfilterWaveformADCNanosec.clear();
        }
    };
}

template struct beamlinereco::hit_t<double>;
template class beamlinereco::CFDHitFinder<double>;

template<class T>
void beamlinereco::CFDHitFinder<T>::SetWaveform(std::vector<uint16_t>& waveform, unsigned int channelNo, uint32_t timestamp) {
    Reset();

    SetChannel(channelNo);
    SetTimestamp(timestamp);

    double polaritySignFactor = _CFDParamSet[kIsWaveformNegativePolarity] ? 1. : -1.;

    T* nonfilterwaveform = (T*) malloc(_CFDParamSet[kNSamplingPoints] * sizeof(T));
    T* filterwaveform    = (T*) malloc(_CFDParamSet[kNSamplingPoints] * sizeof(T));

    for (auto itr = waveform.begin(); itr != waveform.end(); itr++) {
        _NonfilterWaveformADCNanosec.push_back((T)(*itr) * polaritySignFactor);
        nonfilterwaveform[itr - waveform.begin()] = (T)(*itr) * polaritySignFactor;
    }

    if (_CFDParamSet[kGSFilter]) {
        SGSmoothing::Smooth(_CFDParamSet[kNSamplingPoints],
                            nonfilterwaveform,
                            filterwaveform,
                            _CFDParamSet[kGSFilterWindow],
                            _CFDParamSet[kGSFilterDegree]);
    } else {
        for (auto itr = waveform.begin(); itr != waveform.end(); itr++) {
            filterwaveform[itr - waveform.begin()] = (T)(*itr) * polaritySignFactor;
        }
    }

    for (size_t i = 0; i < _CFDParamSet[kNSamplingPoints]; i++) {
        _WaveformADCNanosec.push_back(*(filterwaveform + i));
    }
}

template<class T>
void beamlinereco::CFDHitFinder<T>::SetParam(CFDParams param, T value) {
    _CFDParamSet[param] = value;
}

template<class T>
void beamlinereco::CFDHitFinder<T>::SetParams(const std::map<CFDParams, double>& paramSet) {
    for (auto & itr : paramSet) {
        SetParam(itr.first, itr.second);
    }
}

template<class T>
void beamlinereco::CFDHitFinder<T>::SetChannel(unsigned int channelNo) {
    _ChannelNumber = channelNo;
}

template<class T>
void beamlinereco::CFDHitFinder<T>::SetTimestamp(uint32_t ts) {
    _Timestamp = ts;
}

template<class T>
void beamlinereco::CFDHitFinder<T>::FindPedestal() {
/* 2nd way
    TH1D* pedHist = new TH1D("", "", 10000, 0, 5000);
    for (size_t i = 700; i < 1024; i++) {
        pedHist->Fill(_NonfilterWaveformADCNanosec.at(i));
    }

    _PedestalInADC   = (T) pedHist->GetMean();
    _NoiseSigmaInADC = (T) pedHist->GetStdDev();

    delete (pedHist);
*/

/* 3rd way*/
    T NoiseBand = 1E9;
    T OldPedestal = 1E9;
    T Pedestal  = 0;

    while (TMath::Abs(Pedestal - OldPedestal) >= 2) {
        T average = 0;
        T stdDev = 0;
        size_t countSamples = 0;
        for (typename std::vector<T>::iterator itr = _NonfilterWaveformADCNanosec.begin(); itr != _NonfilterWaveformADCNanosec.end(); itr++) {
            if ((*itr <= Pedestal + NoiseBand) && (*itr >= Pedestal - NoiseBand)) {
                average += *itr;
                countSamples++;
            }
        }
        average = average / (T)countSamples;
        for (typename std::vector<T>::iterator itr = _NonfilterWaveformADCNanosec.begin(); itr != _NonfilterWaveformADCNanosec.end(); itr++) {
            if ((*itr <= Pedestal + NoiseBand) && (*itr >= Pedestal - NoiseBand)) {
                stdDev += TMath::Power((*itr) - average, 2);
            }
        }

        stdDev = TMath::Sqrt(stdDev / (T)countSamples);
        OldPedestal = Pedestal;
        Pedestal = average;
        NoiseBand = stdDev;
    }

    _PedestalInADC = Pedestal;
    _NoiseSigmaInADC = NoiseBand;
/**/
}

template<class T>
void beamlinereco::CFDHitFinder<T>::FindRawHitLogic() {
    T rawHitThreshold = (_PedestalInADC) - ((_NoiseSigmaInADC) * _CFDParamSet[kRawHitFinderThresholdInNoiseSigma]);
    //std::cout << "RawHitThreshold: " << rawHitThreshold << std::endl;
    for (size_t i = 0; i < _WaveformADCNanosec.size(); i++) {
        _RawHitLogicNanosec.push_back(_WaveformADCNanosec.at(i) < rawHitThreshold);
    }

    for (size_t i = 1; i < _RawHitLogicNanosec.size() - 1; i++) {
        if ((not _RawHitLogicNanosec.at(i - 1)) and _RawHitLogicNanosec.at(i)) {
            size_t countDuration = 0;
            for (size_t j = i; j < _RawHitLogicNanosec.size() - 1; j++) {
                if (_RawHitLogicNanosec.at(j)) countDuration++;
                if (not(_RawHitLogicNanosec.at(j + 1))) break;
            }

            if (countDuration < _CFDParamSet[kShortRawHitIgnoringDurationInTicks]) {
                for (size_t j = i; j < i + countDuration; j++) {
                    _RawHitLogicNanosec.at(j) = false;
                }
            }
        }
    }

    for (size_t i = 1; i < _RawHitLogicNanosec.size() - 1; i++) {
        if (_RawHitLogicNanosec.at(i - 1) && not(_RawHitLogicNanosec.at(i))) {
            size_t countDuration = 0;
            for (size_t j = i; j < _RawHitLogicNanosec.size() - 1; j++) {
                if (not(_RawHitLogicNanosec.at(j))) countDuration++;
                if (_RawHitLogicNanosec.at(j + 1)) break;
            }

            if (countDuration < _CFDParamSet[kConsecutiveHitSeperationDurationInTicks]) {
                for (size_t j = i; j < i + countDuration; j++) {
                    _RawHitLogicNanosec.at(j) = true;
                }
            }
        }
    }

}

template<class T>
void beamlinereco::CFDHitFinder<T>::FindCFDHits() {
    std::vector<T> hitPeakValueVector; hitPeakValueVector.clear();
    std::vector<size_t> hitPeakTimeIndexVector; hitPeakTimeIndexVector.clear();
    std::vector<T> hitStartTimeIndexVector; hitStartTimeIndexVector.clear();
    std::vector<T> hitRiseTimeInIndexVector; hitRiseTimeInIndexVector.clear();
    std::vector<T> hitFallTimeInIndexVector; hitFallTimeInIndexVector.clear();
    std::vector<bool> isHitContainedVector; isHitContainedVector.clear();

    bool logicFallingPart = true;
    for (size_t i = 1; i < _RawHitLogicNanosec.size() - 1; i++) {
        T hitPeakValue = 1E10;
        size_t hitPeakTimeIndex = 0;
        if (_RawHitLogicNanosec.at(i) && !(_RawHitLogicNanosec.at(i - 1)) && logicFallingPart) {
            logicFallingPart = false;
            for (size_t j = i; j < _RawHitLogicNanosec.size() - 1; j++) {
                hitPeakValue = hitPeakValue > _WaveformADCNanosec.at(j) ? (T) _WaveformADCNanosec.at(j) : hitPeakValue;
                hitPeakTimeIndex = hitPeakValue == _WaveformADCNanosec.at(j) ? j : hitPeakTimeIndex;
                if (_RawHitLogicNanosec.at(j) && !(_RawHitLogicNanosec.at(j + 1)) && !logicFallingPart) {
                    logicFallingPart = true;
                    T hitStartTimeIndex = (T)hitPeakTimeIndex;
                    T hitRiseTimeInIndex = (T)hitPeakTimeIndex;
                    T hitFallTimeInIndex = (T)hitPeakTimeIndex;
                    bool containedFromRisingEdge = BackwardFindingOfHitStart(hitPeakTimeIndex, hitPeakValue, hitStartTimeIndex, hitRiseTimeInIndex);
                    bool containedFromFallingEdge = ForwardFindingOfHitFallTime(hitPeakTimeIndex, hitFallTimeInIndex);

                    hitPeakValueVector.push_back(hitPeakValue);
                    hitPeakTimeIndexVector.push_back(hitPeakTimeIndex);
                    hitStartTimeIndexVector.push_back(hitStartTimeIndex);
                    hitRiseTimeInIndexVector.push_back(hitRiseTimeInIndex);
                    hitFallTimeInIndexVector.push_back(hitFallTimeInIndex);
                    isHitContainedVector.push_back(containedFromRisingEdge && containedFromFallingEdge);
                    break;
                }
            }
        }
    }

    for (size_t i = 0; i < hitStartTimeIndexVector.size(); i++) {
        hit_t<double> aHit;
        aHit.Timestamp                       = _Timestamp;
        aHit.DigitizerChannel                = _ChannelNumber;
        aHit.TStartInNanoSec                 = hitStartTimeIndexVector.at(i) * _CFDParamSet[kTimeSamplingInterval];
        aHit.TPeakInNanoSec                  = hitPeakTimeIndexVector.at(i) * _CFDParamSet[kTimeSamplingInterval];
        aHit.AmplitudeInMiliVolt             = hitPeakValueVector.at(i) * (_CFDParamSet[kADCDynamicRange] / (TMath::Power(2, _CFDParamSet[kADCNBits]) - 1));
        aHit.AmplitudeInADC                  = hitPeakValueVector.at(i);
        aHit.RiseTimeInNanoSec               = hitRiseTimeInIndexVector.at(i) * _CFDParamSet[kTimeSamplingInterval];
        aHit.FallTimeInNanoSec   	         = hitFallTimeInIndexVector.at(i) * _CFDParamSet[kTimeSamplingInterval];
        aHit.IsContained                     = isHitContainedVector.at(i);
        aHit.IntegratedChargeInADCTimeTicks  = IntegrateWaveformInADC(hitPeakTimeIndexVector.at(i) - hitRiseTimeInIndexVector.at(i),
                                                                      hitPeakTimeIndexVector.at(i) + hitFallTimeInIndexVector.at(i));
        aHit.IntegratedChargeInADCNanoSec    = aHit.IntegratedChargeInADCTimeTicks * _CFDParamSet[kTimeSamplingInterval];

        _HitCollection[hitStartTimeIndexVector.at(i)] = aHit;
    }
}

template<class T>
T beamlinereco::CFDHitFinder<T>::IntegrateWaveformInADC(size_t hitStartTimeIndex, size_t hitEndTimeIndex) {

    T integration = 0;

    if (!_CFDParamSet[kIntergratedWindowFixed]) {
        hitStartTimeIndex = hitStartTimeIndex < 0 ? 0 : hitStartTimeIndex;
        hitEndTimeIndex   = hitEndTimeIndex > (_CFDParamSet[kNSamplingPoints] - 1) ? _CFDParamSet[kNSamplingPoints] - 1 : hitEndTimeIndex;
        for (size_t i = hitStartTimeIndex; i <= hitEndTimeIndex; i++) {
            integration = integration + (_WaveformADCNanosec.at(i) - _PedestalInADC);
        }
    } else {
        for (size_t i = _CFDParamSet[kIntergratedWindowLowerLimitIndex]; i <= _CFDParamSet[kIntergratedWindowUpperLimitIndex]; i++) {
            integration = integration + (_WaveformADCNanosec.at(i) - _PedestalInADC);
        }
    }

    return integration;
}

template<class T>
bool beamlinereco::CFDHitFinder<T>::BackwardFindingOfHitStart(size_t hitPeakTimeIndex, T hitPeakValue, T& hitStartTimeIndex, T& hitRiseTimeInIndex) {
    bool isThisHitContainedFromTheRisingEdge = true;
    _CFDThresholdInADC = _PedestalInADC - _CFDParamSet[kCFDThreshold] * (_PedestalInADC - hitPeakValue);

    size_t tmp_index = hitPeakTimeIndex;
    while (_WaveformADCNanosec.at(tmp_index) < _CFDThresholdInADC) {
        if (tmp_index == 2 && _WaveformADCNanosec.at(tmp_index - 1) < _CFDThresholdInADC) {
            isThisHitContainedFromTheRisingEdge = false;
            break;
        }
        tmp_index = tmp_index - 1;
    }

    T V1 = _WaveformADCNanosec.at(tmp_index);
    T V2 = _WaveformADCNanosec.at(tmp_index - 1);
    size_t t1 = tmp_index;
    size_t t2 = tmp_index - 1;
    T slope = (V2 - V1) / ((T)t2 - (T)t1);
    hitStartTimeIndex = (_CFDThresholdInADC - V1) / slope + (T)t1;

    if (!isThisHitContainedFromTheRisingEdge) {
        hitRiseTimeInIndex = (T) hitPeakTimeIndex - (T) tmp_index;
    } else {
        tmp_index = hitPeakTimeIndex;
        while (_WaveformADCNanosec.at(tmp_index) <= _PedestalInADC) {
            if (tmp_index == 2 && _WaveformADCNanosec.at(tmp_index - 1) < _PedestalInADC) {
                isThisHitContainedFromTheRisingEdge = false;
                break;
            }
            tmp_index = tmp_index - 1;
        }
        hitRiseTimeInIndex = (T) hitPeakTimeIndex - (T) tmp_index;
    }

    return isThisHitContainedFromTheRisingEdge;
}

template<class T>
bool beamlinereco::CFDHitFinder<T>::ForwardFindingOfHitFallTime(size_t hitPeakTimeIndex, T& hitFallTimeInIndex) {
    bool isThisHitContainedFromTheFallingEdge = true;
    size_t tmp_index = hitPeakTimeIndex;
    while (_WaveformADCNanosec.at(tmp_index) <= _PedestalInADC) {
        if (tmp_index == _CFDParamSet[kNSamplingPoints] - 2 && _WaveformADCNanosec.at(tmp_index + 1) < _PedestalInADC) {
            isThisHitContainedFromTheFallingEdge = false;
            break;
        }
        tmp_index = tmp_index + 1;
    }

    hitFallTimeInIndex = (T)tmp_index - (T)hitPeakTimeIndex;
    return isThisHitContainedFromTheFallingEdge;
}

template<class T>
void beamlinereco::CFDHitFinder<T>::Go() {
    if (_WaveformADCNanosec.size() == 0) return;

    FindPedestal(); // Get pedestal and noise band sigma
    FindRawHitLogic(); // From pedestal and noise band, get rawHitThreshold and RawHitLogic vector
    FindCFDHits(); // From RawHitLogic vector, get hitPeakValue and hitPeakTimeIndex. From hitPeakValue and pedestal, get cfdHitThresholdValue. Go back from hitPeakTimeIndex to the point when
}






























// GOLAY_SAVITZKY ALGORITHM

//! default convergence
static const double TINY_FLOAT = 1.0e-300;

//! comfortable array of doubles
using float_vect = std::vector<double>;
//! comfortable array of ints;
using int_vect = std::vector<int>;

/*! matrix class.
 *
 * This is a matrix class derived from a vector of float_vects.  Note that
 * the matrix elements indexed [row][column] with indices starting at 0 (c
 * style). Also note that because of its design looping through rows should
 * be faster than looping through columns.
 *
 * \brief two dimensional floating point array
 */
class float_mat : public std::vector<float_vect> {
private:
    //! disable the default constructor
    explicit float_mat() {};
    //! disable assignment operator until it is implemented.
    float_mat &operator =(const float_mat &) { return *this; };
public:
    //! constructor with sizes
    float_mat(const size_t rows, const size_t cols, const double def=0.0);
    //! copy constructor for matrix
    float_mat(const float_mat &m);
    //! copy constructor for vector
    float_mat(const float_vect &v);

    //! use default destructor
    // ~float_mat() {};

    //! get size
    size_t nr_rows(void) const { return size(); };
    //! get size
    size_t nr_cols(void) const { return front().size(); };
};



// constructor with sizes
float_mat::float_mat(const size_t rows,const size_t cols,const double defval)
        : std::vector<float_vect>(rows) {

    for (unsigned int i = 0; i < rows; ++i) {
        (*this)[i].resize(cols, defval);
    }
    if ((rows < 1) || (cols < 1)) {
        char buffer[1024];

        sprintf(buffer, "cannot build matrix with %d rows and %d columns\n",
                (int)rows, (int)cols);
        //sgs_error(buffer);
    }
}

// copy constructor for matrix
float_mat::float_mat(const float_mat &m) : std::vector<float_vect>(m.size()) {

    float_mat::iterator inew = begin();
    float_mat::const_iterator iold = m.begin();
    for (/* empty */; iold < m.end(); ++inew, ++iold) {
        const size_t oldsz = iold->size();
        inew->resize(oldsz);
        const float_vect oldvec(*iold);
        *inew = oldvec;
    }
}

// copy constructor for vector
float_mat::float_mat(const float_vect &v)
        : std::vector<float_vect>(1) {

    const size_t oldsz = v.size();
    front().resize(oldsz);
    front() = v;
}

//////////////////////
// Helper functions //
//////////////////////

//! permute() orders the rows of A to match the integers in the index array.
void permute(float_mat &A, int_vect &idx)
{
    int_vect i(idx.size());
    unsigned int j,k;

    for (j = 0; j < A.nr_rows(); ++j) {
        i[j] = j;
    }

    // loop over permuted indices
    for (j = 0; j < A.nr_rows(); ++j) {
        if (i[j] != idx[j]) {

            // search only the remaining indices
            for (k = j+1; k < A.nr_rows(); ++k) {
                if (i[k] ==idx[j]) {
                    std::swap(A[j],A[k]); // swap the rows and
                    i[k] = i[j];     // the elements of
                    i[j] = idx[j];   // the ordered index.
                    break; // next j
                }
            }
        }
    }
}

/*! \brief Implicit partial pivoting.
 *
 * The function looks for pivot element only in rows below the current
 * element, A[idx[row]][column], then swaps that row with the current one in
 * the index map. The algorithm is for implicit pivoting (i.e., the pivot is
 * chosen as if the max coefficient in each row is set to 1) based on the
 * scaling information in the vector scale. The map of swapped indices is
 * recorded in swp. The return value is +1 or -1 depending on whether the
 * number of row swaps was even or odd respectively. */
static int partial_pivot(float_mat &A, const size_t row, const size_t col,
                         float_vect &scale, int_vect &idx, double tol)
{
    if (tol <= 0.0)
        tol = TINY_FLOAT;

    int swapNum = 1;

    // default pivot is the current position, [row,col]
    unsigned int pivot = row;
    double piv_elem = fabs(A[idx[row]][col]) * scale[idx[row]];

    // loop over possible pivots below current
    unsigned int j;
    for (j = row + 1; j < A.nr_rows(); ++j) {

        const double tmp = fabs(A[idx[j]][col]) * scale[idx[j]];

        // if this elem is larger, then it becomes the pivot
        if (tmp > piv_elem) {
            pivot = j;
            piv_elem = tmp;
        }
    }

#if 0
    if(piv_elem < tol) {
      //sgs_error("partial_pivot(): Zero pivot encountered.\n")
#endif

    if(pivot > row) {           // bring the pivot to the diagonal
        j = idx[row];           // reorder swap array
        idx[row] = idx[pivot];
        idx[pivot] = j;
        swapNum = -swapNum;     // keeping track of odd or even swap
    }
    return swapNum;
}

/*! \brief Perform backward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is upper
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the lower triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_backsubst(float_mat &A, float_mat &a, bool diag=false)
{
    int r,c;
    unsigned int k;

    for (r = (A.nr_rows() - 1); r >= 0; --r) {
        for (c = (A.nr_cols() - 1); c > r; --c) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if(!diag) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Perform forward substitution.
 *
 * Solves the system of equations A*b=a, ASSUMING that A is lower
 * triangular. If diag==1, then the diagonal elements are additionally
 * assumed to be 1.  Note that the upper triangular elements are never
 * checked, so this function is valid to use after a LU-decomposition in
 * place.  A is not modified, and the solution, b, is returned in a. */
static void lu_forwsubst(float_mat &A, float_mat &a, bool diag=true)
{
    unsigned int r, k, c;
    for (r = 0;r < A.nr_rows(); ++r) {
        for(c = 0; c < r; ++c) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] -= A[r][c] * a[c][k];
            }
        }
        if(!diag) {
            for (k = 0; k < A.nr_cols(); ++k) {
                a[r][k] /= A[r][r];
            }
        }
    }
}

/*! \brief Performs LU factorization in place.
 *
 * This is Crout's algorithm (cf., Num. Rec. in C, Section 2.3).  The map of
 * swapped indeces is recorded in idx. The return value is +1 or -1
 * depending on whether the number of row swaps was even or odd
 * respectively.  idx must be preinitialized to a valid set of indices
 * (e.g., {1,2, ... ,A.nr_rows()}). */
static int lu_factorize(float_mat &A, int_vect &idx, double tol=TINY_FLOAT)
{
    if ( tol <= 0.0)
        tol = TINY_FLOAT;

    if ((A.nr_rows() == 0) || (A.nr_rows() != A.nr_cols())) {
        //sgs_error("lu_factorize(): cannot handle empty "
        //           "or nonsquare matrices.\n");

        return 0;
    }

    float_vect scale(A.nr_rows());  // implicit pivot scaling
    unsigned int i;
    int j;
    for (i = 0; i < A.nr_rows(); ++i) {
        double maxval = 0.0;
        for (j = 0; j < (int)A.nr_cols(); ++j) {
            if (fabs(A[i][j]) > maxval)
                maxval = fabs(A[i][j]);
        }
        if (maxval == 0.0) {
            //sgs_error("lu_factorize(): zero pivot found.\n");
            return 0;
        }
        scale[i] = 1.0 / maxval;
    }

    int swapNum = 1;
    unsigned int c,r;
    for (c = 0; c < A.nr_cols() ; ++c) {            // loop over columns
        swapNum *= partial_pivot(A, c, c, scale, idx, tol); // bring pivot to diagonal
        for(r = 0; r < A.nr_rows(); ++r) {      //  loop over rows
            int lim = (r < c) ? r : c;
            for (j = 0; j < lim; ++j) {
                A[idx[r]][c] -= A[idx[r]][j] * A[idx[j]][c];
            }
            if (r > c)
                A[idx[r]][c] /= A[idx[c]][c];
        }
    }
    permute(A,idx);
    return swapNum;
}

/*! \brief Solve a system of linear equations.
 * Solves the inhomogeneous matrix problem with lu-decomposition. Note that
 * inversion may be accomplished by setting a to the identity_matrix. */
static float_mat lin_solve(const float_mat &A, const float_mat &a,
                           double tol=TINY_FLOAT)
{
    float_mat B(A);
    float_mat b(a);
    int_vect idx(B.nr_rows());
    unsigned int j;

    for (j = 0; j < B.nr_rows(); ++j) {
        idx[j] = j;  // init row swap label array
    }
    lu_factorize(B,idx,tol); // get the lu-decomp.
    permute(b,idx);          // sort the inhomogeneity to match the lu-decomp
    lu_forwsubst(B,b);       // solve the forward problem
    lu_backsubst(B,b);       // solve the backward problem
    return b;
}

///////////////////////
// related functions //
///////////////////////

//! Returns the inverse of a matrix using LU-decomposition.
static float_mat invert(const float_mat &A)
{
    const int n = A.size();
    float_mat E(n, n, 0.0);
    float_mat B(A);
    int i;

    for (i = 0; i < n; ++i) {
        E[i][i] = 1.0;
    }

    return lin_solve(B, E);
}

//! returns the transposed matrix.
static float_mat transpose(const float_mat &a)
{
    float_mat res(a.nr_cols(), a.nr_rows());
    unsigned int i,j;

    for (i = 0; i < a.nr_rows(); ++i) {
        for (j = 0; j < a.nr_cols(); ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

//! matrix multiplication.
float_mat operator *(const float_mat &a, const float_mat &b)
{
    float_mat res(a.nr_rows(), b.nr_cols());
    if (a.nr_cols() != b.nr_rows()) {
        //sgs_error("incompatible matrices in multiplication\n");
        return res;
    }

    unsigned int i,j,k;

    for (i = 0; i < a.nr_rows(); ++i) {
        for (j = 0; j < b.nr_cols(); ++j) {
            double sum(0.0);
            for (k = 0; k < a.nr_cols(); ++k) {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}


//! calculate savitzky golay coefficients.
static float_vect sg_coeff(const float_vect &b, const size_t deg)
{
    const size_t rows(b.size());
    const size_t cols(deg + 1);
    float_mat A(rows, cols);
    float_vect res(rows);

    // generate input matrix for least squares fit
    unsigned int i,j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            A[i][j] = pow(double(i), double(j));
        }
    }

    float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));

    for (i = 0; i < b.size(); ++i) {
        res[i] = c[0][0];
        for (j = 1; j <= deg; ++j) {
            res[i] += c[j][0] * pow(double(i), double(j));
        }
    }
    return res;
}

/*! \brief savitzky golay smoothing.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.  The needed coefficients are
 * generated dynamically by doing a least squares fit on a "symmetric" unit
 * vector of size 2w+1, e.g. for w=2 b=(0,0,1,0,0). evaluating the polynome
 * yields the sg-coefficients.  at the border non symmectric vectors b are
 * used. */
float_vect beamlinereco::sg_smooth(const float_vect &v, const int width, const int deg)
{
    float_vect res(v.size(), 0.0);
    if ((width < 1) || (deg < 0) || ((int)v.size() < (2 * width + 2))) {
        //sgs_error("sgsmooth: parameter error.\n");
        return res;
    }

    const int window = 2 * width + 1;
    const int endidx = v.size() - 1;

    // do a regular sliding window average
    int i,j;
    if (deg == 0) {
        // handle border cases first because we need different coefficients
#if defined(_OPENMP)
#pragma omp parallel for private(i,j) schedule(static)
#endif
        for (i = 0; i < (int)width; ++i) {
            const double scale = 1.0/double(i+1);
            const float_vect c1(width, scale);
            for (j = 0; j <= i; ++j) {
                res[i]          += c1[j] * v[j];
                res[endidx - i] += c1[j] * v[endidx - j];
            }
        }

        // now loop over rest of data. reusing the "symmetric" coefficients.
        const double scale = 1.0/double(window);
        const  float_vect c2(window, scale);
#if defined(_OPENMP)
#pragma omp parallel for private(i,j) schedule(static)
#endif
        for (i = 0; i <= (int) (v.size() - window); ++i) {
            for (j = 0; j < (int)window; ++j) {
                res[i + width] += c2[j] * v[i + j];
            }
        }
        return res;
    }

    // handle border cases first because we need different coefficients
#if defined(_OPENMP)
#pragma omp parallel for private(i,j) schedule(static)
#endif
    for (i = 0; i < (int) width; ++i) {
        float_vect b1(window, 0.0);
        b1[i] = 1.0;

        const float_vect c1(sg_coeff(b1, deg));
        for (j = 0; j < (int) window; ++j) {
            res[i]          += c1[j] * v[j];
            res[endidx - i] += c1[j] * v[endidx - j];
        }
    }

    // now loop over rest of data. reusing the "symmetric" coefficients.
    float_vect b2(window, 0.0);
    b2[width] = 1.0;
    const float_vect c2(sg_coeff(b2, deg));

#if defined(_OPENMP)
#pragma omp parallel for private(i,j) schedule(static)
#endif
    for (i = 0; i <= (int) (v.size() - window); ++i) {
        for (j = 0; j < (int) window; ++j) {
            res[i + width] += c2[j] * v[i + j];
        }
    }
    return res;
}

/*! least squares fit a polynome of degree 'deg' to data in 'b'.
 *  then calculate the first derivative and return it. */
static float_vect lsqr_fprime(const float_vect &b, const int deg)
{
    const int rows(b.size());
    const int cols(deg + 1);
    float_mat A(rows, cols);
    float_vect res(rows);

    // generate input matrix for least squares fit
    int i,j;
    for (i = 0; i < (int) rows; ++i) {
        for (j = 0; j < (int) cols; ++j) {
            A[i][j] = pow(double(i), double(j));
        }
    }

    float_mat c(invert(transpose(A) * A) * (transpose(A) * transpose(b)));

    for (i = 0; i < (int) b.size(); ++i) {
        res[i] = c[1][0];
        for (j = 1; j < (int) deg; ++j) {
            res[i] += c[j + 1][0] * double(j+1)
                      * pow(double(i), double(j));
        }
    }
    return res;
}

/*! \brief savitzky golay smoothed numerical derivative.
 *
 * This method means fitting a polynome of degree 'deg' to a sliding window
 * of width 2w+1 throughout the data.
 *
 * In contrast to the sg_smooth function we do a brute force attempt by
 * always fitting the data to a polynome of degree 'deg' and using the
 * result. */
float_vect beamlinereco::sg_derivative(const float_vect &v, const int width,
                                       const int deg, const double h)
{
    float_vect res(v.size(), 0.0);
    if ((width < 1) || (deg < 1) || ((int) v.size() < (2 * width + 2))) {
        //sgs_error("sgsderiv: parameter error.\n");
        return res;
    }

    const int window = 2 * width + 1;

    // handle border cases first because we do not repeat the fit
    // lower part
    float_vect b(window, 0.0);
    int i,j;

    for (i = 0; i < (int) window; ++i) {
        b[i] = v[i] / h;
    }
    const float_vect c(lsqr_fprime(b, deg));
    for (j = 0; j <= (int) width; ++j) {
        res[j] = c[j];
    }
    // upper part. direction of fit is reversed
    for (i = 0; i < (int) window; ++i) {
        b[i] = v[v.size() - 1 - i] / h;
    }
    const float_vect d(lsqr_fprime(b, deg));
    for (i = 0; i <= (int) width; ++i) {
        res[v.size() - 1 - i] = -d[i];
    }

    // now loop over rest of data. wasting a lot of least squares calcs
    // since we only use the middle value.
#if defined(_OPENMP)
#pragma omp parallel for private(i,j) schedule(static)
#endif
    for (i = 1; i < (int) (v.size() - window); ++i) {
        for (j = 0; j < (int) window; ++j) {
            b[j] = v[i + j] / h;
        }
        res[i + width] = lsqr_fprime(b, deg)[width];
    }
    return res;
}

void beamlinereco::SGSmoothing::Smooth(size_t nsamples, double *in, double *out, const int w, const int deg) {
    std::vector<double> in_vec;
    for (size_t idx = 0; idx < nsamples; idx++) {
        in_vec.push_back(*(in + idx));
    }

    std::vector<double> out_vec = sg_smooth(in_vec, w, deg);

    for (size_t idx = 0; idx < nsamples; idx++) {
        *(out + idx) = out_vec.at(idx);
    }
}

void beamlinereco::SGSmoothing::Derivative(size_t nsamples, double *in, double *out, const int w, const int deg, const double h) {
    std::vector<double> in_vec;
    for (size_t idx = 0; idx < nsamples; idx++) {
        in_vec.push_back(*(in + idx));
    }

    std::vector<double> out_vec = sg_derivative(in_vec, w, deg, h);

    for (size_t idx = 0; idx < nsamples; idx++) {
        *(out + idx) = out_vec.at(idx);
    }
}




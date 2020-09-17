#ifndef _WFDAnalysis3_
#define _WFDAnalysis3_

#include <iterator>
#include <cassert>
#include <algorithm>
#include <memory>

// ---- List of function ---- //
//
// template <class InputIterator>
// double GetBaseline(InputIterator first, InputIterator last)
//
// template <class InputIterator, class OutputIterator>
// void GetCalibratedWaveform(InputIterator first1, InputIterator last1,
//                            OutputIterator result,
//                            double scale,
// 			      int baseline_begin, int baseline_end);
//
// template <class InputIterator>
// double GetAmplitude(InputIterator first, InputIterator last);
//



class WFDAnalysis{
 public:
  WFDAnalysis(){};
  ~WFDAnalysis(){};
  template <class InputIterator>
    double GetBaseline(InputIterator first, int sampling_num);
      
  template <class InputIterator, class OutputIterator>
    void GetCalibratedWaveform(InputIterator first1, InputIterator last1,
			       OutputIterator result,
			       double scale,
			       int baseline_begin, int baseline_end);

  template <class InputIterator>
    double GetAmplitude(InputIterator first, InputIterator last);
 private:
  std::shared_ptr<TGraph> g_amp;
    
};
  
  template <class InputIterator>
    double WFDAnalysis::GetBaseline(InputIterator first, int sampling_num){
    double acc = 0;
    auto last = first + sampling_num;
    std::size_t size = std::distance(first, first + sampling_num);
    assert(size > 0);
    while(first != last)
      acc += *first++;
    return acc / size;
  }

  template <class InputIterator, class OutputIterator>
    void WFDAnalysis::GetCalibratedWaveform(InputIterator first1, InputIterator last1,
					    OutputIterator result,
					    double scale,
					    int baseline_begin, int baseline_end){
    double baseline = GetBaseline(first1 + baseline_begin,
				  baseline_end - baseline_begin + 1);
    std::transform(first1, last1, result, [&](decltype(*first1) x){return (x - baseline) * scale; });
    return;
  }

  template <class InputIterator>
    double WFDAnalysis::GetAmplitude(InputIterator first, InputIterator last){
    std::size_t x = std::distance(first, last);
    assert(x > 0);
    auto it_min = std::min_element(first + 1, last);
    auto i_min  = std::distance(first, it_min);
    // -- Fit function -- //
    std::shared_ptr<TF1> pol2(new TF1("pol2", "[2] * (x - [1])^2 + [0]") );
    pol2 -> SetParameter(1, i_min);
    pol2 -> SetParameter(0, *it_min);
    // -- TGraph -- //
    g_amp = std::shared_ptr<TGraph>(new TGraph() );
    int ith_point = 0;
    for(auto it = first ; it != last ; it++){
      g_amp -> SetPoint(ith_point, ith_point, *it);
      ith_point++;
    }//end of <it>
    g_amp -> Fit("pol2", "Q", "", i_min - 3, i_min + 3);
    return -1 * ( pol2 -> GetParameter(0) );
  }

#endif

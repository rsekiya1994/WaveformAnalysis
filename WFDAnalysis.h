#ifndef _WFDAnalysis3_
#define _WFDAnalysis3_

#include <iterator>
#include <cassert>
#include <algorithm>
#include <memory>
#include <TAxis.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>

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
  WFDAnalysis(double sampling_frequency);
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
  
  template <class CalibratedWaveformInputIterator>
  double CFD(CalibratedWaveformInputIterator first,
	     CalibratedWaveformInputIterator last,
	     double constant, int index_delay);
  
  template <class InputIterator>
  void ShowRawWaveform(InputIterator first,
		       InputIterator last ,
		       std::string add = "",
		       std::string ns_scale = ""); // add = "add", ns_scale = "ns"
  
  template <class InputIterator>
  void ShowCalibratedWaveform(InputIterator first,
			      InputIterator last ,
			      std::string add = "");// add = "add"

  void ShowAmplitudeFitting(std::string add = "");
  void ShowCFDWaveform(std::string add = "");
  
  private:
  double frequency_;  // GHz
  double ns_factor_;
  std::vector<std::shared_ptr<TCanvas> > c_draw;
  std::shared_ptr<TGraph> g_raw;
  std::shared_ptr<TGraph> g_calib;
  std::shared_ptr<TGraph> g_amp;
  std::shared_ptr<TF1> pol2_amplitude;
  std::shared_ptr<TGraph> g_cfd;
  std::shared_ptr<TF1> pol2_cfd;
};

WFDAnalysis::WFDAnalysis(double sampling_frequency){
  frequency_ = sampling_frequency;
  ns_factor_  = 1. / frequency_;
}

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
  pol2_amplitude = std::make_shared<TF1>("pol2_amp", "[2] * (x - [1])^2 + [0]");
  pol2_amplitude -> SetParameter(1, i_min);
  pol2_amplitude -> SetParameter(0, *it_min);
  // -- TGraph -- //
  g_amp = std::shared_ptr<TGraph>(new TGraph() );
  int ith_point = 0;
  for(auto it = first ; it != last ; it++){
    g_amp -> SetPoint(ith_point, ith_point, *it);
    ith_point++;
  }//end of <it>
  g_amp -> Fit("pol2_amp", "Q", "", i_min - 3, i_min + 3);
  return -1 * ( pol2_amplitude -> GetParameter(0) );
}

template <class CalibratedWaveformInputIterator>
double WFDAnalysis::CFD(CalibratedWaveformInputIterator first,
	   CalibratedWaveformInputIterator last,
	   double constant, int index_delay){
  // ---- Size investigation and some preparation ---- //
  int size = (int)(last - first);
  std::vector<double> wfd_cfd(size);
  std::vector<double> x(size);
  // ---- Make CFD waveform ---- //
  for(int i = 0 ; i < size ; i++){
    if( i <= index_delay ){
      wfd_cfd[i] = 0;
    }else{
      wfd_cfd[i] = ( *(first + i - index_delay) ) - ( *(first + i) ) * constant;
    }// if
    x[i] = i * ns_factor_;
  }// i
  
  // ---- CFD Graph Definition ---- //
  g_cfd    = std::make_shared<TGraph>(size, x.data(), wfd_cfd.data() );
  pol2_cfd = std::make_shared<TF1>("pol2_cfd", "[2] * (x - [1])^2 + [0]");
  
  // ---- Search Minimum Point ---- //
  auto it_min = std::min_element(wfd_cfd.begin(), wfd_cfd.end());
  auto i_min  = it_min - wfd_cfd.begin();
  
  // ---- Search Maximum Point ---- //
  auto it_max = std::max_element(wfd_cfd.begin(), wfd_cfd.end());
  auto i_max  = it_max - wfd_cfd.begin();

  // ---- Exception Judgement ---- //
  // --> Value of it_max should be positive
  // --> and
  // --> it_max should come before it_min
 
  if( *it_max < 0 ){
    return - 11000;
  }
  if( (it_max > it_min) ){
    return - 10000;
  }
  
  // ---- Search Zero Cross Point ---- //
  auto it_zerocross = std::find_if(it_max, it_min, [](double const &y){ return y <= 0; } );
  auto i_zerocross  = it_zerocross - wfd_cfd.begin();
  
  // ---- Fitting near zerocross points ---- //
  pol2_cfd -> SetParameter(1, i_max * ns_factor_);
  pol2_cfd -> SetRange(i_max * ns_factor_, (i_zerocross + 1) * ns_factor_);
  g_cfd -> Fit("pol2_cfd", "Q", "", i_max * ns_factor_, (i_zerocross + 1) * ns_factor_);
  
  // ---- Get Fitting Result ---- //
  double p0 = pol2_cfd -> GetParameter(0);
  double p1 = pol2_cfd -> GetParameter(1);
  double p2 = pol2_cfd -> GetParameter(2);
  
  double D = - p0 / p2;
  // ---- D (= b^2 - 4ac) should be positive ---- //
  if( D < 0 ){
    return -12000;
  }
  // ---- solve 2nd order polynomial ---- //
  double ans1 = p1 + std::sqrt(D);
  double ans2 = p1 - std::sqrt(D);
  double distance_zerocross1 = std::fabs(ans1 - i_zerocross * ns_factor_);
  double distance_zerocross2 = std::fabs(ans2 - i_zerocross * ns_factor_);

  if(distance_zerocross1 <= distance_zerocross2){
    return ans1;
  }else{
    return ans2;
  }
  return -13000;
}

template <class InputIterator>
void WFDAnalysis::ShowRawWaveform(InputIterator first,
				  InputIterator last,
				  std::string add,
				  std::string ns_scale){
  if(add == "add"){
    c_draw.push_back(nullptr);
  }
  if(c_draw.size() == 0){
    c_draw.push_back(nullptr);
  }
  int c_ith = (int)c_draw.size() - 1;
  c_draw[c_ith] = std::make_shared<TCanvas>("c_raw", "c_raw", 500, 400);
  g_raw = std::make_shared<TGraph>();
  int ith_point = 0;
  double scale = 1;
  if(ns_scale == "ns"){
     scale = ns_factor_;
  }
  for(auto it = first ; it != last ; it++){
    g_raw -> SetPoint(ith_point, ith_point * scale, *it);
    ++ith_point;
  }
  g_raw -> SetMarkerStyle(20);
  g_raw -> SetMarkerSize(0.5);
  g_raw -> SetTitle("Raw Waveform");
  if(ns_scale == "ns"){
     g_raw -> GetXaxis() -> SetTitle("time [ns]");
  }else{
    g_raw -> GetXaxis() -> SetTitle("sampling number");
  }
  c_draw[c_ith] -> Draw();
  g_raw -> GetYaxis() -> SetTitle("");
  g_raw -> Draw("AP");
}

template <class InputIterator>
void WFDAnalysis::ShowCalibratedWaveform(InputIterator first,
					 InputIterator last ,
					 std::string add){
  if(add == "add"){
    c_draw.push_back(nullptr);
  }
  if(c_draw.size() == 0){
    c_draw.push_back(nullptr);
  }
  int c_ith = (int)c_draw.size() - 1;
  c_draw[c_ith] = std::make_shared<TCanvas>("c_calib", "c_calib", 500, 400);
  g_calib = std::make_shared<TGraph>();
  int ith_point = 0;
  for(auto it = first ; it != last ; it++){
    g_calib -> SetPoint(ith_point, ith_point * ns_factor_, *it);
    ++ith_point;
  }
  g_calib -> SetMarkerStyle(20);
  g_calib -> SetMarkerSize(0.5);
  g_calib -> SetTitle("Calibrated Waveform");
  g_calib -> GetXaxis() -> SetTitle("time [ns]");
  g_calib -> GetYaxis() -> SetTitle("voltage [mV]");
  g_calib -> Draw("AP");
}

void WFDAnalysis::ShowAmplitudeFitting(std::string add){
  if(add == "add"){
    c_draw.push_back(nullptr);
  }
  if(c_draw.size() == 0){
    c_draw.push_back(nullptr);
  }
  int c_ith = (int)c_draw.size() - 1;
  c_draw[c_ith] = std::make_shared<TCanvas>("c_amp_fit", "c_amp_fit", 500, 400);
  g_amp -> SetMarkerStyle(20);
  g_amp -> SetMarkerSize(0.5);
  g_amp -> SetTitle("Waveform with Amplitude Fitting");
  g_amp -> GetXaxis() -> SetTitle("time [ns]");
  g_amp -> GetYaxis() -> SetTitle("");
  g_amp -> Draw("AP");
  pol2_amplitude -> Draw("same");
}

void WFDAnalysis::ShowCFDWaveform(std::string add){
  if(add == "add"){
    c_draw.push_back(nullptr);
  }
  if(c_draw.size() == 0){
    c_draw.push_back(nullptr);
  }
  int c_ith = (int)c_draw.size() - 1;
  c_draw[c_ith] = std::make_shared<TCanvas>("c_cfd", "c_cfd", 500, 400);
  g_cfd -> SetMarkerStyle(20);
  g_cfd -> SetMarkerSize(0.5);
  g_cfd -> SetTitle("CFD Waveform");
  g_cfd -> GetXaxis() -> SetTitle("time [ns]");
  g_cfd -> GetYaxis() -> SetTitle("");
  g_cfd -> Draw("AP");
  pol2_cfd -> Draw("same");
}

#endif

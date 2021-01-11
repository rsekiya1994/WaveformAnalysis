# Waveform Data Analysis
Some library to analyze waveform data


The class WFDAnalysis code runs upon ROOT developed by CERN.
Also it requires c++11.
This is a header only library.

Definition
--
```c++
namespace sky {
  class WFDAnalysis;
}
```

Member functions
--
```c++
WFDAnalysis(int sampling_frequency);
```
The constructor. 
The unit of the parameter `sampling_frequency` is [GHz].
For example, the sampling rate is 2.5 GHz, you should give 2.5 to the argument.

---
```c++
~WFDAnalysis();
```
The destructor.

---
```c++
template <class InputIterator>
double GetBaseline(InputIterator first, int sampling_num);
```
This function calculates a baseline of an input waveform. It takes average of the continuous `sampling_num` points of values from the `first` points.

---

# Xenon_Generic_Calibration
Matlab script and example files for analyzing hyperpolarized 129Xe calibration data 

## Authorship
Author: Peter J. Niedbalski 

Contact: pniedbalski@kumc.edu

Written: June 21, 2021

## Main Function: calibration_analysis.m
This function is written to analyze dissolved and gas FIDs acquired using a hyperpolarized 129Xe calibration spectroscopy sequence.

Several values are hardcoded in this function to correspond to the parameters recommended by the 129Xe MRI clinical trials consortium. These include:
- Dwell Time
- Echo Time
- Flip Angle
- Dissolved Frequency Offset

As written, this function is designed for data acquired at 3T. However, values for 1.5T are included as comments.

Fitting is performed in the time domain using code provided in:

Robertson SH, et al. Uncovering a third dissolved-phase 129Xe resonance in the human lung: Quantifying spectroscopic features in healthy subjects and patients with idiopathic pulmonary fibrosis. Magn Reson Med 2017;78(4):1306-1315.

### Inputs
- disfids: all dissolved FIDs, shaped as (NPts x NFIDs)
- gasfids: all gas FIDs, shaped as (NPts x NFIDs)

### Outputs: 
There is a single structure output. This structure has fields:
- Freq_Offset: The true gas center frequency in relation to the set frequency (in Hz) - i.e. to get the center frequency for future scans, set frequency to Current_Frequency+Freq_Offset
- Set2Act_Flip: Ratio of prescribed to actual flip angle
- TE90: The true TE90
- RBC2TP: Ratio of xenon RBC signal to tissue/plasma signal.

## Usage
Different scanner platforms and different sequence implementations will have different data types and different formats for output data. As such, we leave it up to each individual site to utilize the appropriate code for reading in data and sorting out the dissolved and the gaseous FIDs.

Call the function by providing dissolved and gas FIDs. As calculations are performed, a figure is written to the screen. This figure shows the values that need to be set based on this calibration sequence:
- Frequency offset: The frequency (in Hz) of the gas-phase signal in relation to the currently set center frequency. Future scans should be set to "Current_Frequency + Frequency Offset"
- Prescribed to Actual Flip Angle Ratio: To ensure that the scanner plays out excitation pulses with the proper flip angle, either flip angle or scanner reference voltage needs to be multiplied by this value.
- TE90: Echo time to be set for 1-point Dixon Gas Exchange Imaging
- RBC/TP Ratio: Ratio of RBC signal to Tissue/Plasma signal. Used for analysis of Gas Exchange data. Also, useful for real-time quality control (RBC/TP should be between 0.1 and 0.7 - if not, could indicate poor fitting).

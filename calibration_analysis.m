function Output_Struct = calibration_analysis(disfids,gasfids)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate Gas Center Frequency, TE90, RBC/TP ratio, and
% Actual to applied flip angle
%
% Note that values for Dwell Time, Echo Time, Flip Angle, Dissolved
% Frequency Offset are hardcoded (at the top of the code) to match the
% recommendations provided in the 129Xe MRI Clinical Trials Consortium
% Position Paper for imaging at 3T(Niedbalski et al. Submitted to Magnetic Resonance in
% Medicine)
%
% Fitting is performed in the time domain using code provided in 
% Robertson SH, et al. Uncovering a third dissolved-phase 129Xe resonance 
% in the human lung: Quantifying spectroscopic features in healthy subjects 
% and patients with idiopathic pulmonary fibrosis. Magn Reson Med 2017;78(4):1306-1315.
%
% Different scanner platforms and different sequence implementations will
% have different data types and different formats for output data. As such,
% we leave it up to each individual site to utilize the appropriate code
% for reading in data and sorting out the dissolved and the gaseous FIDs.
%
% Inputs: disfids - all dissolved FIDs, shaped as (NPts x NFIDs)
%         gasfids - all gas FIDs, shaped as (NPts x NFIDs)
%
% Outputs: Output_Struct - structure with fields:
%                *Freq_Offset - The true gas center frequency in relation
%                to the set frequency - i.e. to get the center frequency
%                for future scans, set frequency to
%                Current_Frequency+Freq_Offset.
%                *Set2Act_Flip - Ratio of prescribed to actual flip angle
%                *TE90 - True TE90
%                *RBC2TP - RBC/TP ratio
%
% Author: Peter J. Niedbalski 
% Contact: pniedbalski@kumc.edu
% Written: June 21, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters according to 129Xe CTC Position Paper
%Dwell Time in s
dwell_time = 39e-6/2;
%From dwell time, get time vector
t = double((0:(length(disfids)-1))*dwell_time');
%Echo Time in ms
te = 0.45;
% te = 0.8 %For 1.5T
%Flip angle in degrees
FA = 20;
%Dissolved Offset in Hz
Dis_Freq = 7430;
% Dis_Freq = 3715; %For 1.5T

%% Fit gas spectra for center frequency 
% only fit first gas spectrum
gas1 = gasfids(:,1);
% Set initial guesses
area_guess = 1e-4;
freq_guess = 0;
fwhm_guess = 30;
phase_guess = 0;
line_broadening = 0;
zeroPadSize = 10000;
% Fit gas FID in the time domain
gasfitObj = Spectroscopy.NMR_TimeFit(gas1,t,area_guess,freq_guess,fwhm_guess,phase_guess,line_broadening,zeroPadSize);
gasfitObj.fitTimeDomainSignal();
% From gas fit object, get offset frequency
Output_Struct.Freq_Offset = gasfitObj.freq(1);

%% Use time-domain max amplitude for fitting flip angle
max_amp = max(abs(gasfids));
% Equation 1 from position paper
gas_decay_fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1)+coefs(3);
guess(1) = max(max_amp);
guess(2) = 15*pi/180; %Guess 15 degrees
guess(3) = 0;

xdata = 1:length(max_amp);
ydata = max_amp;
% Fit decay of gas signal - just looking for flip angle here, so can ignore
% most outputs
fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,~,~,~,~,~,~]  = lsqcurvefit(gas_decay_fitfunct,guess,xdata,ydata,[],[],fitoptions);
% Write to output
Act_Flip = fitparams(2)*180/pi;
Output_Struct.Set2Act_Flip = FA/Act_Flip;
%% Display a Figure
Cal_Fig = figure('Name','Calibration Output');
set(Cal_Fig,'color','white','Units','normalized','Position',[0.1 0.1 0.8 0.5]);
subplot(1,2,1)
plot(xdata,ydata,'bo','MarkerFaceColor','b');
hold on;
plot(xdata,gas_decay_fitfunct(fitparams,xdata),'-r');
hold off
xlabel('FID Number');
ylabel('Magnitude');
title('Decay of Gas Signal')
set(gca,'FontSize',14);
annotation('textbox',[0.25 0.7 0.2 0.2],'String',{['Gas Frequency Offset = ' num2str(gasfitObj.freq(1),'%.0f') ' Hz'];['Prescribed Flip Angle = ' num2str(FA,'%.0f') '\circ'];['Actual Flip Angle = ' num2str(Act_Flip,'%.1f') '\circ'];['Ratio Prescribed/Actual = ' num2str(Output_Struct.Set2Act_Flip,'%.3f')]},'FitBoxToText','on','FontSize',14,'BackgroundColor','white')

%% Fit Dissolved Spectra
% To get adequate signal, average dissolved spectra together
disfid_avg = mean(disfids,2);
% Set initial guesses for fitting - Each spectrum should have 3 components
% - RBC, Tissue/Plasma, Gas - so pass 3 vectors of guesses. index 1 is
% guess for RBC, index 2 for Tissue/plasma, index 3 for gas
area_guess = [1 1 1];
freq_guess = [0 -700 -7400];
% freq_guess = [0 -350 -3700]; %for 1.5T
fwhm_guess = [250 200 30];
fwhmG_guess = [0 200 0];
phase_guess = [0 0 0];
line_broadening = 0;
zeroPadSize = length(t); 

disfitObj = Spectroscopy.NMR_TimeFit_v(disfid_avg,t,area_guess,freq_guess,fwhm_guess,fwhmG_guess,phase_guess,line_broadening,zeroPadSize); 
disfitObj = disfitObj.fitTimeDomainSignal();
%Get Final fitting in Frequency Domain for ease of viewing
disfitFinal = dwell_time*fftshift(fft(disfitObj.calcComponentTimeDomainSignal(t),[],1),1);

deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-barrier freq difference
deltaTe90= (90-deltaPhase)/(360*deltaF); % how far off are we?
Output_Struct.TE90 = te + deltaTe90; % in msec

Output_Struct.RBC2TP = disfitObj.area(1)/disfitObj.area(2); 

%% Finish second half of figure.
subplot(1,2,2)
plot(disfitObj.f, abs(disfitObj.spectralDomainSignal),'k')%mag
hold('on');
plot(disfitObj.f,abs(sum(disfitFinal,2)),'Color',[0 0 1 0.33],'LineWidth',3)%mag
hold('off');
title('Dissolved Phase Spectrum and Fit')
set(gca,'Xdir','reverse','FontSize',14)
xlabel('Frequency (Hz)')
ylabel('NMR Signal (a.u.)')
xlim([-10000 5000]);
annotation('textbox',[0.79 0.7 0.2 0.2],'String',{['TE90 = ' num2str(Output_Struct.TE90,'%.2f') ' ms'];['RBC/TP = ' num2str(disfitObj.area(1)/disfitObj.area(2),'%.2f')]},'FitBoxToText','on','FontSize',14,'BackgroundColor','white')

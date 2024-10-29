%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================== SetUpL7_4PA_Delay_BModeToggle.m =====================

% Last modified: 2/14/2020
% Author: Vinoin Devpaul Vincely
% Adapted from Verasonics example script by Jonah Harmon

% Defaults to interleaved B-mode and photoacoustic imaging. Captures one
% full B-mode frame using multiangled plane wave transmits, and then
% captures a number of PA frames to generate one coherently compounded
% final image (number of acqs to compound is specified by ne; can set to 1
% if higher frame rate is desired). Users can subsequently compound
% additional frames in post-processing by adding IQ frames.

% A conditional event sequence is specified to allow for GUI-initated
% saving of IQ data, both B-mode and PA. The script will save 
% Resource.RcvBuffer(1).numFrames IQ frames, with each PA IQ frame being
% the composite of ne acqs. So, if you want to effectively compound 200 PA
% acquisitions, set ne = 5, and Resource.RcvBuffer(1).numFrames = 40. You
% can then read the separate files into Matlab, add the IQ data from all 40
% frames, and have basically the same end result as ne = 200. Just depends
% on what time resolution you need, as after they are compounded, you can't
% go back. Setting ne low and Resource.RcvBuffer(1).numFrames high
% preserves more data.

% A toggle button allows the user to save the full RF frame in the same mat
% file as the PA IQ data. This RF frame includes all B-mode acquisitions
% and PA acqs in the same 2D matrix. To facilitate separation of the
% different acquisitions and isolation of PA from B-mode data, a few
% parameters are saved alongside the RF and IQ data. The start sample, end
% sample, and length of each PA acquisition (in samples) are saved. So, to
% get only PA data, you can index RF, e.g. PADat =
% RF(PAStartSample:PAEndSample, :); and can subsequently separate out PA
% acqs for beamforming.

% As an overview of the default imaging sequence:
%   1) Transmits na angled plane waves to generate one coherently
%      compounded B-mode frame
%   2) Waits for an input trigger from the laser flash lamp
%   3) Does nothing (noop) while the laser pumps, syncs the hardware and
%      software sequencers; user can specify wait time in microseconds
%   4) Sends external trigger to Q-switch to fire laser, immediately starts
%      listening for PA signal
%   5) Repeat 2-4 ne times for full PA acquisition
%   6) Repeat 1-5 Resource.RcvBuffer(1).numFrames times to fill RF buffer
%   7) Transfer RF to host computer
%   8) Recon, display, jump back to event 1



% NOTES PRESERVED FROM VERASONICS EXAMPLE SCRIPT:

% When running in simulation (with oneway=0), the media points are imaged 
% as though the medium moves between B-mode and PA acquisitions, and thus 
% are slightly offset from each other to better visualize each result.

%   - For the PA acquisition, 'ne' T/R events with output trigger for an 
%       external source (laser) are performed, and use receive-only 
%       reconstruction (when oneway=1).

%       * The sequence waits for an input trigger (e.g. from laser flash 
%           lamp), then pauses while the laser pumps (PR.flash2Qdelay 
%           microsecs), then fires the Q-switch using a trigger out that 
%           also begins receive acquisition. This approach results in very 
%           low jitter between the laser pulse and the ultrasound.

%       * Multiple PA acquisitions are reconstructed and then coherently 
%           averaged by accumulation in the Inter Buffer. Incoherent 
%           averaging is possible by summing in the Image buffer instead, 
%           and relaxes jitter requirements. Though not implemented here, 
%           the RF data could be accumulated in hardware prior to transfer
%           using the Receive.mode attribute; doing so also requires 
%           precise synchronization with the laser.

%       * The receive-only reconstruction mode is not currently supported 
%           in simulation, so to test the PA code we simulate by turning on
%           the transmitters for a conventional T/R acquisition and
%           reconstruction. Toggling between modes is easily done using the
%           'oneway' parameter (=0 turn on transmitters, =1 turn off
%           transmitters for real PA acquisition) below.

%       * Two TGC structures are defined, but only TGC(1) is set up to be 
%           controlled by the GUI sliders. To use the TGC slider controls
%           for the PA acquisitions, use TGC(1) for PA and TGC(2) for 2D.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Move to and activate directory; run with most recent update
% vsrootv4

clear all;
close all;

%% Specify parameters for B-mode and PA - user defined
% Upsample 
upsample = 1; 

% Split to allow some flexibility, but these can be the same if desired
P(1).startDepth = 0;       % Acquisition start depth (in wavelengths)
P(1).endDepth = 165;       % Acquisition end depth (in wavelengths)
P(2).startDepth = 0;       % Acquisition start depth (in wavelengths)
P(2).endDepth = 165;       % Acquisition end depth may be different for PA

% Set PA parameters
oneway = 1; % 0 for simulate mode, 1 for hardware operation

% Time from flash lamp trigIn to VSX trigOut to Q-switch (@ time = 0);
% needs to be 75-230 us accounting for 100 us pulse width q-switch trigIn
% 175 us seems to be the sweet spot! MUST USE FUNCTION GENERATOR FOR TRIG,
% 100 us PULSE WIDTH IS NECESSARY
signalPathDelay = 0; % us; delay from routing through function generator
trigPulseWidth = 0; % us; QSwitch trigger pulse width
PR.trigOutDelay = signalPathDelay + trigPulseWidth; % Time to wait from trig out to listening
PR.flash2Qdelay = 210; % microseconds between trigger input and start of acquisition (which outputs a trigger pulse at time=0)

PR.ne = 5;             % Num of acqs in PA ensemble, coherently compounded in inter buffer
PR.numPAFrames = 4;   % Number of frames captured per burst

%% No need for user input - these will remain same for all operation
% For GUI control over saving data
PR.numImagingEvents = 0;               % For tracking correct startEvent
PR.currentLoc = 1;                     % Current imaging location
newDir = strrep(datestr(now),':','-'); % Current datetime on setup
newDir = strrep(newDir,' ','_');       % Remove spaces
PR.initTime = newDir;                  % Directory for saving data
PR.savingRFToggle = 0;                 % Flag for saving PA RF
PR.tempRF = [];                        % Temp storage for PA RF data
PR.tempBModeIData = [];                % Temp storage for BMode I data
PR.tempBModeQData = [];                % Temp storage for BMode Q data

% Set number of angles for coherent compounding, B-mode
PR.na = 12;     
if PR.na>1
    dtheta2D = (36*pi/180)/(PR.na-1); % Change in angle
else
    dtheta2D = 0; % For single plane wave transmits; much faster acqs
end

% Define start angle
if fix(PR.na/2) == PR.na/2   % If PR.na is even
    startAngle = (-(fix(PR.na/2) - 1) - 0.5)*dtheta2D;
else
    startAngle = -fix(PR.na/2)*dtheta2D;
end

% Basically never used
PA_Angle = 0;       % Angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 10;        % PA PRF in Hz. Only used if not using flash lamp trigger.
                    
% Enable or disable transmit for simulate or true imaging respectively
if oneway == 1
    disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***')
else
    disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
end

%% ========== Specify system parameters
Resource.Parameters.numTransmit = 128;     % Number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % Number of receive channels.
Resource.Parameters.speedOfSound = 1540;   % Set speed of sound in m/sec before calling computeTrans
Resource.Parameters.connector = 2;         % Use L7-4 with second connector
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 2;    % Operate with hardware
% Resource.Parameters.simulateMode = 1;      % Operate in simulate mode

%% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm';             % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % computeTrans is used for known transducers.
nElem = Trans.numelements;
Trans.maxHighVoltage = 50;      % Set a reasonable high voltage limit.

%% Specify PData structure arrays.
% 2D Bmoe PData structure
PData(1).PDelta(1) = 1.0;
PData(1).PDelta(3) = 0.5;
PData(1).Size(1,1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3));   % rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(1,3) = 1;  % Single image page
PData(1).Origin = [-Trans.spacing*63.5,0,P(1).startDepth]; % x,y,z of uppr lft crnr.

% PA PData structure
PData(2).PDelta(1) = 1.0;
PData(2).PDelta(3) = 0.5;
PData(2).Size(1,1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3));   % PA rows
PData(2).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % PA cols
PData(2).Size(1,3) = 1;  % Single image page
PData(2).Origin = [-Trans.spacing*63.5,0,P(2).startDepth]; % x,y,z of upper lft crnr.

%% Specify Media object and point displacement function
% Media.MP(1,:) = [0,0,P(1).startDepth+(P(1).endDepth-P(1).startDepth)/2,1.0];
% Media.attenuation = -0.5;

pt1;
Media.function = 'movePoints';

%% Specify Resources.
% Upsampling 
% warning('\nUpsampling = %d! Save post-collection workspace for reference ...\n', upsample); 

% RcvBuffer(1) is for both 2D and PA acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*2048*(PR.na + PR.ne); % Needs to be larger than intended data length
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = PR.numPAFrames;   % Frame(s) allocated for RF acqusitions.
PR.imageNum = 1-Resource.RcvBuffer(1).numFrames;  % Set to not save burst on startup

% InterBuffer(1) is for 2D Bmode reconstructions.
Resource.InterBuffer(1).numFrames = 1;  % One intermediate frame needed for 2D.

% InterBuffer(2) is for PA reconstructions.
Resource.InterBuffer(2).numFrames = 1;  % One intermediate frame needed for PA.

% ImageBuffer(1) is for 2D Bmode image.
Resource.ImageBuffer(1).datatype = 'double';  % Image buffer for 2D
Resource.ImageBuffer(1).numFrames = 20;

% ImageBuffer(2) is for PA image.
Resource.ImageBuffer(2).datatype = 'double';  % Image buffer for PA
Resource.ImageBuffer(2).numFrames = Resource.ImageBuffer(1).numFrames;

%% Specify Display Windows
% Separate Display Window, B-mode
Resource.DisplayWindow(1).Title = 'B-mode';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics'; % This is fine if no ROI controls needed
Resource.DisplayWindow(1).numFrames = 10;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Separate Display Window, PA
Resource.DisplayWindow(2).Title = 'Photoacoustic Signal';
Resource.DisplayWindow(2).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [750,(ScrnSize(4)-(DwHeight+150))/2, ...  % right of B-mode
                                      DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics'; % This is fine if no ROI controls needed
Resource.DisplayWindow(2).numFrames = 10;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = hot(256);

%% Specify Transmit waveforms structure
% 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% PA transmit waveform - not used in actual imaging, for sim only
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,6,1];

%% Specify Transmit beams structure
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, PR.na+1); % PR.na TXs for 2D + 1 for PA

% Assign angle for each transmit, compute delays
for n = 1:PR.na  
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% PA TX struct - all zero apodization if doing actual imaging
TX(PR.na+1).waveform = 2;
if oneway
    TX(PR.na+1).Apod = zeros(1,Trans.numelements); % Specifies receive only beamformer
else
    % This is the conventional T/R condition and invokes the default beamformer
    TX(PR.na+1).Apod = ones(1,Trans.numelements);  
end

% Assign steer angle, compute delays
TX(PR.na+1).Steer = [PA_Angle,0.0];            % Only used for sim mode
TX(PR.na+1).Delay = computeTXDelays(TX(PR.na+1)); % Only used for sim mode

%% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 50;

% This allows one to use different transmit profile for PA. Can't really
% see a case where this is used unless moving to something like
% sonophotoacoustics, but no reason not to leave it here.
TPC(2).name = 'PA';
TPC(2).maxHighVoltage = 35;

%% Analog front end gain settings.
% Sets the gain level for the initial low noise analog amplifier
RcvProfile(1).LnaGain = 18;   % 15, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

RcvProfile(2).LnaGain = 24;   % Additional amplification of PA signal
RcvProfile(2).condition = 'immediate';

%% Specify Receive structure arrays.
% We need to acquire all the 2D and PA data within a single RcvBuffer frame. 
% This allows the transfer-to-host DMA after each frame to transfer a large
% amount of data, improving throughput.

% We need PR.na Receives for a 2D frame and PR.ne Receives for a PA frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthPA = sqrt(P(2).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);    % wavelengths in a 128 sample block for 4 smpls per wave round trip.
% wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
                        'TGC', 1, ...     % TGC(1) is tied to the GUI sliders; B-mode
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (PR.na+PR.ne)*Resource.RcvBuffer(1).numFrames);
               
% Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (PR.na + PR.ne)*(i-1); % k keeps track of Receive index increment per frame.
    
    % Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1; % Move points before doing ensemble of different angle plane waves
    
    % Acquisitions for 2D
    for j = 1:PR.na
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    
    % PA acquisitions
    for j = (PR.na+1):(PR.na+PR.ne)
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;   % PA acqNums continue after 2D
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
        Receive(j+k).TGC = 2;      % TGC(1) is tied to the GUI sliders
        % move points between 2D and PA to see difference in simulation
        if j==PR.na+1, Receive(j+k).callMediaFunc = 1; end 
    end
end

% Save some info for separating out RF frame components
PR.PAStart = PR.na*2*4*Receive(1).endDepth + 1; % First sample after B-mode acqs
lastRcv = (PR.na+PR.ne)*Resource.RcvBuffer(1).numFrames;
PR.PAAcqLength = 2*4*(Receive(lastRcv).endDepth-Receive(lastRcv).startDepth);  % Find length of each PA acq
PR.PAEnd = PR.ne * PR.PAAcqLength + PR.PAStart;   % For finding last PA sample in RF frame

%% Specify TGC Waveform structures.
% 2D TGC
TGC(1).CntrlPts = [1023 1023 1023 1023 1023 1023 1023 1023];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

% PA TGC
TGC(2).CntrlPts = round([1023 1023 1023 1023 1023 1023 1023 1023]./2); % User needs to adjust gain appropriately
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

%% Specify Recon structure arrays.
% We need two Recon structures, one for 2D, one for PA. These will be 
% referenced in the same event, so that they will use the same (most 
% recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 0), 1, 2);
           
% Set Recon values for 2D frame.
Recon(1).RINums(1:PR.na) = 1:PR.na;  % PR.na ReconInfos needed for PR.na angles
k = PR.na + 1;

% Set Recon values for PA ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums(1:PR.ne) = k:(k+PR.ne-1);  % 'PR.ne' ReconInfos needed for PA ensemble.

%% Define ReconInfo structures.
% For 2D, we need PR.na ReconInfo structures for PR.na steering angles.
% For PA, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, PR.na + PR.ne);
               
% ReconInfos for 2D frame.
ReconInfo(1).mode = 'replaceIQ';
for j = 1:PR.na
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
if PR.na > 1
    % For multiangle transmit, add last IQ frame, replace intensity
    ReconInfo(PR.na).mode = 'accumIQ_replaceIntensity';
else
    % For single plane wave, just replace intensity
    ReconInfo(PR.na).mode = 'replaceIntensity';
end

% ReconInfos for PA ensemble.
k = PR.na;
for j = 1:PR.ne
    if j==1, ReconInfo(k+j).mode = 'replaceIQ'; end
    ReconInfo(k+j).txnum = PR.na + 1;
    ReconInfo(k+j).rcvnum = PR.na + j;
end
if PR.ne > 1
    ReconInfo(PR.na+PR.ne).mode = 'accumIQ_replaceIntensity';
else
    ReconInfo(PR.na+PR.ne).mode = 'replaceIntensity';
end

%% Specify Process structure arrays.
% Define display parameters here for use in UI control later
cpers = 80; % Persistence
rejlv = 10; % Reject level

% B-mode
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',10,...     % reject level
                         'persistMethod','simple',...
                         'persistLevel',50,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',... % 'lowerHalf' for split colormap
                         'display',1,...     % Display image after processing
                         'displayWindow',1};

% PA
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',rejlv,...  % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',cpers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',... % 'upperHalf' for split colormap
                         'display',1,...     % Display image after processing
                         'displayWindow',2};
                     
% Plot PA RF data 
Process(6).classname = 'External'; 
Process(6).method = 'PlotPA_RF'; 
Process(6).Parameters = {'srcbuffer', 'receive', ... % name of the buffer to process 
                        'srcbufnum', 1, ...
                        'srcframenum', -1, ... %% Process the most recent frame
                        'dstbuffer', 'none'};
                     
% Save B-mode IQ frame
Process(3).classname = 'External';
Process(3).method = 'saveBModeIQ';
Process(3).Parameters = {'srcbuffer', 'inter', 'srcbufnum', 1, ... % Grab 2D IQ frame
                         'srcframenum', 1, 'srcpagenum', 1, ...    % 1 frame
                         'dstbuffer', 'none'};                     % Not sending data anywhere

% Save PA IQ frame
Process(4).classname = 'External';
Process(4).method = 'savePAIQ';
Process(4).Parameters = {'srcbuffer', 'inter', 'srcbufnum', 2, ... % Grab PA IQ frame
                         'srcframenum', 1, 'srcpagenum', 1, ...    % 1 frame
                         'dstbuffer', 'none'};                     % Not sending data anywhere
                     
% Save RF frame
Process(5).classname = 'External';
Process(5).method = 'saveRFData';
Process(5).Parameters = {'srcbuffer', 'receive', 'srcbufnum', 1, ... % Grab RF frame
                         'srcframenum', -1, 'dstbuffer', 'none'};    % Not sending data anywhere
                     
%% Specify SeqControl structure arrays.
% Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200; % 5000 Hz between B-mode angles

% Change to Profile 2 (PA); not used since no PA transmit pulse
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;

% Time between 2D acquisition and PA ensemble. Set to allow time for Rcv profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 7000; % time in usec

% PRF for PA ensemble; not used since triggering VSX with flash lamp
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(PA_PRF*1e-06)); % (10 msecs for PA_PRF=100 Hz)

% Change to Profile 1 (2D); not used
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;

% Time between PA and next 2D acquisition. Set to allow time for Rcv profile change.
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = 7000; % time in usec

% Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;

% Set receive profile - for selecting amplifier gain pre-ADC
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;

% Output trigger; 1 us, 3.3V, falling edge
SeqControl(10).command = 'triggerOut';

% Input trigger - timeout range is 1:255 in 250 msec steps; 0 means timeout
% disabled; may have to change condition depending on trig from flash lamp
SeqControl(11).command = 'triggerIn';
SeqControl(11).condition = 'Trigger_1_Rising'; % Trigger input 1, enable with rising edge
SeqControl(11).argument = 2; % 500 msec timeout delay
    
% Noop delay between trigger in and start of acquisition
SeqControl(12).command = 'noop';
SeqControl(12).argument = fix(PR.flash2Qdelay)*5; % noop counts are in 0.2 microsec increments

% Sync command
SeqControl(13).command = 'sync';

% Noop to account for trigger pulse width
SeqControl(14).command = 'noop';
SeqControl(14).argument = fix(PR.trigOutDelay)*5;

% Conditional branch to allow GUI toggling of B-mode only, for higher frame
% rate while finding target tissue, and B-mode + PA for actually localizing
% before burst acquisitions. Set event to jump to later.
SeqControl(15).command = 'cBranch';
SeqControl(15).condition = 'bFlag';

% Conditional branch to toggle back to B-mode only; basically conditional
% jump back to Event 1
SeqControl(16).command = 'cBranch';
SeqControl(16).condition = 'bFlag';
SeqControl(16).argument = 1;

% Set to 17; defining an additional jump SeqControl later in script
nsc = 18;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

%% Specify default, B-mode only imaging sequence, higher frame rate
n = 1;

% Need to sync hardware and software sequencers since we aren't really sure
% which is doing what at the branch point
Event(n).info = 'sync after branch';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 13;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:PR.na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (PR.na+PR.ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [3, nsc];   % Longer TTNA between frames, transfer frame
    SeqControl(nsc).command = 'transferToHost'; % Transfer frame to host buffer
      nsc = nsc+1;
    
    Event(n).info = 'recon and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor  % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
    
    Event(n).info = 'Conditional Jump (cbranch)';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 15;
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;
n = n+1;
% keyboard
%% Specify B-mode + PA sequence, toggled with GUI control
% This is the start of the cBranch Flash image acquisition and processing
% events
SeqControl(15).argument = n;  % branch Event for the cBranch command
SeqControl(17).command = 'jump';
SeqControl(17).argument = n;  % Define jump to allow continuous B-mode + PA

% Need to sync hardware and software sequencers since we aren't really sure
% which is doing what at the branch point
Event(n).info = 'sync after branch';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 13;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:PR.na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (PR.na+PR.ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [3,9];   % replace last 2D acquisition Event's seqControl (longer TTNA and new RCV profile)

    % Acquire PA ensemble.
    for j = (PR.na+1):(PR.na+PR.ne)
        % Wait for input trigger from flash lamp firing
        Event(n).info = 'Wait for Trigger IN';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 11;
        n = n+1;

        % Pause for optical buildup
        Event(n).info = 'noop and sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [12,13];
        n = n+1;
        
        % Fire trigger - KEEP ONLY IF USING FUNCTION GENERATOR TO EMIT 100
        % US PULSE WIDTH TRIG
        Event(n).info = 'trigger out';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 10;
        n = n+1;
        
        % Wait for signal path, trigger width - KEEP ONLY IF USING FUNCTION
        % GENERATOR
        Event(n).info = 'noop for trigger pulse width';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 14;
        n = n+1;

        % Send trigger output at start of every PA acquisition to fire Q-switch
        % CHANGE SEQCONTROL TO 10 IF NOT USING FUNCTION GENERATOR
        Event(n).info = 'Acquire PA event';
        Event(n).tx = PR.na+1;
        Event(n).rcv = (PR.na+PR.ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 0;
        n = n+1;
    end
    Event(n-1).seqControl = [6,8]; % replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change

    Event(n).info = 'Plot RF data'; 
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 6; % External Function, plot RF 
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % Transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'PA image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
    
    Event(n).info = 'Conditional Jump (cbranch)';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 16;
    n = n+1;
end

Event(n).info = 'Jump back to B-mode + PA';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 17;
n = n+1;

% For setting correct start event
PR.numImagingEvents = n-1;
Resource.Parameters.startEvent = PR.numImagingEvents+1;  % Start with burst

%% Specify conditional sequence - acq, plus save IQ frame
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:PR.na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (PR.na+PR.ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [3,9];   % replace last 2D acquisition Event's seqControl (longer TTNA and new RCV profile)

    % Acquire PA ensemble.
    for j = (PR.na+1):(PR.na+PR.ne)
        % Wait for input trigger from flash lamp firing
        Event(n).info = 'Wait for Trigger IN';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 11;
        n = n+1;

        % Pause for optical buildup
        Event(n).info = 'noop and sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [12,13];
        n = n+1;
        
        % Fire trigger - KEEP ONLY IF USING FUNCTION GENERATOR TO EMIT 100
        % US PULSE WIDTH TRIG
        Event(n).info = 'trigger out';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 10;
        n = n+1;
        
        % Wait for signal path, trigger width - KEEP ONLY IF USING FUNCTION
        % GENERATOR
        Event(n).info = 'noop for trigger pulse width';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 14;
        n = n+1;

        % Send trigger output at start of every PA acquisition to fire Q-switch
        % CHANGE SEQCONTROL TO 10 IF NOT USING FUNCTION GENERATOR
        Event(n).info = 'Acquire PA event';
        Event(n).tx = PR.na+1;
        Event(n).rcv = (PR.na+PR.ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 0;
        n = n+1;
    end
    Event(n-1).seqControl = [10,6,8]; % replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change
    
    Event(n).info = 'Plot RF data'; 
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 6; % External Function, plot RF 
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % Transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'Recon B-mode and PA';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'PA image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;
    
    % Process event - grab RF
    Event(n).info = 'Grab RF';
    Event(n).tx = 0;            % No transmit
    Event(n).rcv = 0;           % No receive
    Event(n).recon = 0;         % No recon
    Event(n).process = 5;       % External function, grab RF frame
    Event(n).seqControl = 0;    % No need for SeqControl
    n = n+1;
    
    % Process event - grab IQ, B-mode
    Event(n).info = 'Grab IQ, B-mode';
    Event(n).tx = 0;            % No transmit
    Event(n).rcv = 0;           % No receive
    Event(n).recon = 0;         % No recon
    Event(n).process = 3;       % External function, save B-mode IQ frame
    Event(n).seqControl = 0;    % No need for SeqControl
    n = n+1;
    
    % Process event - save all desired data
    Event(n).info = 'Save data';
    Event(n).tx = 0;            % No transmit
    Event(n).rcv = 0;           % No receive
    Event(n).recon = 0;         % No recon
    Event(n).process = 4;       % External function, save PA IQ frame
    Event(n).seqControl = 0;    % No need for SeqControl
    n = n+1;
end

Event(n).info = 'Jump back to B-mode + PA';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 17;

%% User specified UI Control Elements
% 2D and PA acq, save IQ data
UI(1).Control = {'UserB1','Style','VsPushButton','Label','Save Burst',...
                 'Label','SaveImageSequence'};
UI(1).Callback = text2cell('%BurstAcqCallback');

% Move to next imaging location
UI(2).Control = {'UserB2','Style','VsPushButton','Label','Next Burst',...
                 'Label','NextImagingLocation'};
UI(2).Callback = text2cell('%NextLocCallback');

% Sensitivity Cutoff
UI(3).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(3).Callback = text2cell('%SensCutCallback');

% Color Persistence Slider
UI(4).Control = {'UserB6','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,cpers],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ColPersCallback');

% Reject Level Slider
UI(5).Control = {'UserB5','Style','VsSlider','Label','Thresholding','SliderMinMaxVal',[0,100,rejlv],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(5).Callback = text2cell('%RejLevelCallback');

% Save PA RF (on) or only save IQ (off)
UI(6).Control = {'UserB3','Style','VsToggleButton','Label','Save RF...',...
                 'Label','saverf'};
UI(6).Callback = text2cell('%RFFrameToggle');

% B-mode only / B-mode + PA cBranch toggle
UI(7).Control = {'UserC2','Style','VsPushButton','Label','Switch Modes...'};
UI(7).Callback = text2cell('%CBranchButtonCallback');

% RF channel for display 
nr = Resource.Parameters.numRcvChannels; 
UI(8).Control = {'UserC8','Style','VsSlider','SliderMinMaxVal',[1,nr,64], 'Label', 'RF Channel'};
UI(8).Callback = {'assignin(''base'',''myPlotChn1'',round(UIValue))'};

%% Specify external functions
EF(1).Function = text2cell('%SaveIQ2D');
EF(2).Function = text2cell('%SaveIQPA');
EF(3).Function = text2cell('%SaveRFFrame');
EF(4).Function = vsv.seq.function.ExFunctionDef('PlotPA_RF', @PlotPA_RF); 

%% Conclude setup
% Save all the structures to a .mat file
save('C:\Users\verasonics\Documents\Vantage-4.9.2-2308102000\MatFiles\L7-4PA_PumpLaser.mat');
return





%% UI(1) - Burst acquisition callback function
%BurstAcqCallback
Resource = evalin('base', 'Resource');
PR = evalin('base', 'PR');
Resource.Parameters.startEvent = PR.numImagingEvents+1;
assignin('base', 'Resource', Resource);

% Set Control command to update Event sequence
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Event'};
assignin('base','Control', Control);

assignin('base', 'PR', PR);

return
%BurstAcqCallback

%% UI(2) - Next location callback function
%NextLocCallback
PR = evalin('base', 'PR');
PR.currentLoc = PR.currentLoc + 1;
PR.imageNum = 1; % Reset for new location
disp('Ready to image next location')
assignin('base', 'PR', PR);

return
%NextLocCallback

%% UI(3) - Sensitivity cutoff change
%SensCutCallback
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);

return
%SensCutCallback

%% UI(4) - Color persistence change
%ColPersCallback
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'persistLevel'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);

% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'persistLevel',UIValue};
assignin('base','Control', Control);

return
%ColPersCallback

%% UI(5) - Reject level change
%RejLevelCallback
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'reject'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);

% Set Control.Command to set reject level.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'reject',UIValue};
assignin('base','Control', Control);

return
%RejLevelCallback

%% UI(6) - Toggle saving RF frames
%RFFrameToggle
PR = evalin('base', 'PR');
PR.savingRFToggle = UIState; % Set flag; defaults to 0
assignin('base', 'PR', PR);

return
%RFFrameToggle

%% UI(7) - Switch between B-mode only and B-mode + PA
%CBranchButtonCallback
Control = evalin('base','Control');
Control.Command = 'setBFlag';
assignin('base','Control', Control);
isUsingNewHal = com.verasonics.vantage.hal.VantageHal.isUsingNewHal();
if(isUsingNewHal)
    result = com.verasonics.hal.sequencer.Sequencer.takeSubSequenceBranch();
else
    result = setConditionalBranchFlag;
end

return
%CBranchButtonCallback



%% ================ External Functions ====================
%% EF(1) - Grab B-mode IQ
%SaveIQ2D
saveBModeIQ(IData, QData)
PR = evalin('base', 'PR'); % Get P info
Resource = evalin('base', 'Resource'); % Get resource info

% Create unique folder for each run so users don't overwrite their files
% after running script again
PR.acqSessionDir = [pwd, '/PA-Acqs/PA_AcqSession_', PR.initTime];
if ~exist(PR.acqSessionDir, 'dir')
    mkdir(PR.acqSessionDir);
end

% Create subfolder to hold single location
PR.singleLocationDir = [PR.acqSessionDir, sprintf('/BurstSet%03d', PR.currentLoc)];
if ~exist(PR.singleLocationDir, 'dir')
    mkdir(PR.singleLocationDir);
end

% For initial acq on startup...
if PR.imageNum < 1
    return
end

% Grab data in temp matrix for saving later
PR.tempBModeIData = IData;
PR.tempBModeQData = QData;

assignin('base', 'PR', PR);
%SaveIQ2D

%% EF(2) - Save PA IQ
%SaveIQPA
savePAIQ(IData, QData)
PR = evalin('base', 'PR'); % Get P info
Resource = evalin('base', 'Resource'); % Get resource info

% For initial acq on startup...
if PR.imageNum < 1
    PR.imageNum = PR.imageNum + 1;
    assignin('base', 'PR', PR)
    return
end

% Save IQ frame to .mat file
PAIData = IData;
PAQData = QData;
BModeIData = PR.tempBModeIData;
BModeQData = PR.tempBModeQData;
NumBModeAngles = PR.na;
NumPAAcqs = PR.ne;
    
if PR.savingRFToggle == 0
    save([PR.singleLocationDir, sprintf('/Acq%04d-onlyIQ', PR.imageNum)], 'PAIData', 'PAQData', ...
                                                             'BModeIData', 'BModeQData',...
                                                             'NumBModeAngles', 'NumPAAcqs');
elseif PR.savingRFToggle == 1
    RF = PR.tempRF;
    PAStartSample = PR.PAStart;
    PAAcqLength = PR.PAAcqLength;
    PAEndSample = PR.PAEnd;
    save([PR.singleLocationDir, sprintf('/Acq%04d-plusRF', PR.imageNum)], 'PAIData', 'PAQData', ...
                                                             'BModeIData', 'BModeQData', ...
                                                             'RF', 'PAStartSample', ...
                                                             'PAAcqLength', 'PAEndSample',...
                                                             'NumBModeAngles', 'NumPAAcqs');
end

% Notify user, increment image counter
fprintf('PA and BMode acq %03d saved, resume imaging\n', PR.imageNum)
PR.imageNum = PR.imageNum + 1;

assignin('base', 'PR', PR);
%SaveIQPA

%% EF(3) - Grab RF for temp storage, save with PA data
%SaveRFFrame
saveRFData(RData)
PR = evalin('base', 'PR');

% Do nothing if not saving RF data
if PR.savingRFToggle == 0
    return
end

% If toggled on, grab RF data and place in temp location
PR.tempRF = RData;
assignin('base', 'PR', PR);
%SaveRFFrame

%% EF(4) - Plot RF data 
%PlotPA_RF
function PlotPA_RF(RData) 
        persistent myHandle 
        % We begin by defining a persistent variable, ‘myHandle’, which will keep 
        % the handle to our plot window available over multiple calls to our 
        % function

        tempPR = evalin('base', 'PR'); %keyboard
    %     lims = evalin('base', 'Receive(1).endSample');
        % If 'myPlotChn1' exists read it for the channel to plot 
        if evalin('base', 'exist(''myPlotChn1'', ''var'')')
            channel = evalin('base', 'myPlotChn1'); 
        else 
            channel = 64;
        end 

        % Create the figure if it doesn't exist
        if isempty(myHandle) || ~ishandle(myHandle)
            figure; 
            myHandle = axes('NextPlot', 'replacechildren', 'YLim', [-2000 2000]);
        end 
    %     'XLim', [0 temp(1).endSample], ...
        % To clear and update the plot portion of the figure, we use the ‘NextPlot’,
        % ‘replacechildren’ attribute pair

        % Plot the RF data 
        plot(myHandle, RData([0:tempPR.PAAcqLength]+tempPR.PAStart, channel)); 
    %     plot(myHandle, RData(:, channel)); 
        title(myHandle, sprintf('Channel #%d', channel)); 
        drawnow; 
end    
%PlotPA_RF
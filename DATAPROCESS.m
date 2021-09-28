% Process the data
% This script was written for the CyberSAFE@UALR REU 2016 funded by the NSF\
% The data processing performed here is for a specific experiment
% that we designed and implemented
% Writer: Dennis Frank (frankd@gatech.edu)
% Partner: Jasmine Mabrey (j.n.mabrey@spartans.nsu.edu)
% Mentor: Kenji Yoshigoe (kxyoshigoe@ualr.edu)

SubjectID = input('Subject ID? ','s');
RecordID = input('Record ID? ','s');
Day = input('Day? Format: DD.MM.YY ','s');
StartTime = input('Start Time? Format: Military HH.MM.SS ','s');
filename = [SubjectID '-' RecordID '-' Day '-' StartTime '.csv'];
close all

Mn = {'COUNTER',...
    'INTERPOLATED',...
    'AF3',...
    'T7',...
    'Pz',...
    'T8',...
    'AF4',...
    'RAW_CQ',...
    'GYROX',...
    'GYROY',...
    'MARKER',...
    'MARKER_HARDWARE',...
    'SYNC',...
    'TIME_STAMP_s',...
    'TIME_STAMP_ms',...
    'CQ_AF3',...
    'CQ_T7',...
    'CQ_Pz',...
    'CQ_T8',...
    'CQ_AF4'};
IDmain = fopen(filename);
M = csvread(filename,1,0);

Fs = 128;  % SAMPLING FREQUENCY
% So there will be 128 samples per ERP
[Mr, Mc] = size(M);
tt = 0:1/Fs:Mr/Fs;
tt(end) = [];

% Code for plotting raw data from main file
%{
for i = 1:Mc
    figure
    plot(tt, M(:,i));
    title(Mn(i),'Interpreter','None')
    xlabel('Time (s)')
    ylabel('MAGNITUDE')
end
%}

% Find which channel had the best overall connection quality
CQAF3 = sum(M(:,16));
CQT7 = sum(M(:,17));
CQPz = sum(M(:,18));
CQT8 = sum(M(:,19));
CQAF4 = sum(M(:,20));
MaxVec = [CQAF3 CQT7 CQPz CQT8 CQAF4];
[~, MaxChannel] = max(MaxVec);
if MaxChannel == 1
    display('The channel with the best overall connection quality is AF3');
elseif MaxChannel == 2
    display('The channel with the best overall connection quality is T7');
elseif MaxChannel == 3
    display('The channel with the best overall connection quality is Pz');
elseif MaxChannel == 4
    display('The channel with the best overall connection quality is T8');
elseif MaxChannel == 5
    display('The channel with the best overall connection quality is AF4');
end

% Pluck out channels and markers
AF3 = M(:,3)';
T7 = M(:,4)';
Pz = M(:,5)';
T8 = M(:,6)';
AF4 = M(:,7)';
MARKERS = M(:,11)';

% Make sure FOUR baseline markers are present
% TWO for reference EEG data and TWO for challenge EEG data
% Also, find where the '2' marker is in the EEG data to seperate
    % the reference EEG data from the challenge EEG data
bCounter = 0;
Rbase = [0 0];
Cbase = [0 0];
SWITCHloc = 0;
for i = 1:Mr
    if MARKERS(i)==70 % 70 is the marker for baseline start/stop
        if bCounter == 0
            Rbase(1) = i;
        elseif bCounter == 1
            Rbase(2) = i;
        elseif bCounter == 2
            Cbase(1) = i;
        elseif bCounter == 3
            Cbase(2) = i;
        end
        bCounter = bCounter + 1;
    end
    if MARKERS(i) == 2
       SWITCHloc = i; 
    end
end
if bCounter~=4
    error('There is not the correct number of baseline markers!');
end

% Determine mean of baseline EEG data to subtract from raw EEG data
% For reference 
AF3baseR = sum(AF3(Rbase(1):Rbase(2)))/(Rbase(2)-Rbase(1));
T7baseR = sum(T7(Rbase(1):Rbase(2)))/(Rbase(2)-Rbase(1));
PzbaseR = sum(Pz(Rbase(1):Rbase(2)))/(Rbase(2)-Rbase(1));
T8baseR = sum(T8(Rbase(1):Rbase(2)))/(Rbase(2)-Rbase(1));
AF4baseR = sum(AF4(Rbase(1):Rbase(2)))/(Rbase(2)-Rbase(1));

% For challenge
AF3baseC = sum(AF3(Cbase(1):Cbase(2)))/(Cbase(2)-Cbase(1));
T7baseC = sum(T7(Cbase(1):Cbase(2)))/(Cbase(2)-Cbase(1));
PzbaseC = sum(Pz(Cbase(1):Cbase(2)))/(Cbase(2)-Cbase(1));
T8baseC = sum(T8(Cbase(1):Cbase(2)))/(Cbase(2)-Cbase(1));
AF4baseC = sum(AF4(Cbase(1):Cbase(2)))/(Cbase(2)-Cbase(1));

% Splice the channels into reference EEG and challenge EEG
AF3R = AF3(1:SWITCHloc);
AF3C = AF3(SWITCHloc:end);
T7R = T7(1:SWITCHloc);
T7C = T7(SWITCHloc:end);
PzR = AF3(1:SWITCHloc);
PzC = AF3(SWITCHloc:end);
T8R = T8(1:SWITCHloc);
T8C = T8(SWITCHloc:end);
AF4R = AF4(1:SWITCHloc);
AF4C = AF4(SWITCHloc:end);
% Splice the marker channel
MARKERR = MARKERS(1:SWITCHloc);
MARKERC = MARKERS(SWITCHloc:end);

% Subtract baseline from each channel
AF3R = AF3R - AF3baseR;
T7R = T7R - T7baseR;
PzR = PzR - PzbaseR;
T8R = T8R - T8baseR;
AF4R = AF4R - AF4baseR;

AF3C = AF3C - AF3baseC;
T7C = T7C - T7baseC;
PzC = PzC - PzbaseC;
T8C = T8C - T8baseC;
AF4C = AF4C - AF4baseC;

% Filter channels with butterband filter from 1-55Hz
low = 1;
high = 55;
[b,a]=butter(9,[low/(Fs/2),high/(Fs/2)]); % 9th order
AF3R = filter(b,a,AF3R);
T7R = filter(b,a,T7R);
PzR = filter(b,a,PzR);
T8R = filter(b,a,T8R);
AF4R = filter(b,a,AF4R);
AF3C = filter(b,a,AF3C);
T7C = filter(b,a,T7C);
PzC = filter(b,a,PzC);
T8C = filter(b,a,T8C);
AF4C = filter(b,a,AF4C);


% Initialize ERPs for each of the stimuli
% Rows are: AF3, T7, Pz, T8, AF4
% Visual stimuli have 1000 ms ERP (1 sec * 128 sample/sec)
ERP5 = zeros(5,Fs);
ERP6 = zeros(5,Fs); 
ERP7 = zeros(5,Fs);
ERP8 = zeros(5,Fs);
ERP9 = zeros(5,Fs);
% Audio stimuli have 2000 ms ERP (2 sec * 128 sample/sec)
ERP10 = zeros(5,Fs*2);
ERP11 = zeros(5,Fs*2);
ERP12 = zeros(5,Fs*2);

% Calculate the reference ERPs
[~, length] = size(AF3R);
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
% AF3R
for i = 1:length
   if MARKERR(i) == 5 % SMELL
       if trashcheck(AF3R(i:i+127))
            ERP5(1,:) = ERP5(1,:) + AF3R(i:i+127);
            M5 = M5 + 1;
       end
   elseif MARKERR(i) == 6 % FOOD
       if trashcheck(AF3R(i:i+127))
            ERP6(1,:) = ERP6(1,:) + AF3R(i:i+127);
            M6 = M6 + 1;
       end
   elseif MARKERR(i) == 7 % TOOL
       if trashcheck(AF3R(i:i+127))
            ERP7(1,:) = ERP7(1,:) + AF3R(i:i+127);
            M7 = M7 + 1;
       end
   elseif MARKERR(i) == 8 % GESTURE
       if trashcheck(AF3R(i:i+127))
            ERP8(1,:) = ERP8(1,:) + AF3R(i:i+127);
            M8 = M8 + 1;
       end
   elseif MARKERR(i) == 9 % PERSON
       if trashcheck(AF3R(i:i+127))
            ERP9(1,:) = ERP9(1,:) + AF3R(i:i+127);
            M9 = M9 + 1;
       end
   elseif MARKERR(i) == 10 % SONG 1
       if trashcheck(AF3R(i:i+255))
            ERP10(1,:) = ERP10(1,:) + AF3R(i:i+255);
            M10 = M10 + 1;
       end
   elseif MARKERR(i) == 11 % SONG 2
       if trashcheck(AF3R(i:i+255))
            ERP11(1,:) = ERP11(1,:) + AF3R(i:i+255);
            M11 = M11 + 1;
       end
   elseif MARKERR(i) == 12 % SONG 3
       if trashcheck(AF3R(i:i+255))
            ERP12(1,:) = ERP12(1,:) + AF3R(i:i+255);
            M12 = M12 + 1;
       end
   end
end
ERP5(1,:) = ERP5(1,:)/M5;
ERP6(1,:) = ERP6(1,:)/M6;
ERP7(1,:) = ERP7(1,:)/M7;
ERP8(1,:) = ERP8(1,:)/M8;
ERP9(1,:) = ERP9(1,:)/M9;
ERP10(1,:) = ERP10(1,:)/M10;
ERP11(1,:) = ERP11(1,:)/M11;
ERP12(1,:) = ERP12(1,:)/M12;
% T7R
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERR(i) == 5 % SMELL
       if trashcheck(T7R(i:i+127))
            ERP5(2,:) = ERP5(2,:) + T7R(i:i+127);
            M5 = M5 + 1;
       end
   elseif MARKERR(i) == 6 % FOOD
       if trashcheck(T7R(i:i+127))
            ERP6(2,:) = ERP6(2,:) + T7R(i:i+127);
            M6 = M6 + 1;
       end
   elseif MARKERR(i) == 7 % TOOL
       if trashcheck(T7R(i:i+127))
            ERP7(2,:) = ERP7(2,:) + T7R(i:i+127);
            M7 = M7 + 1;
       end
   elseif MARKERR(i) == 8 % GESTURE
       if trashcheck(T7R(i:i+127))
            ERP8(2,:) = ERP8(2,:) + T7R(i:i+127);
            M8 = M8 + 1;
       end
   elseif MARKERR(i) == 9 % PERSON
       if trashcheck(T7R(i:i+127))
            ERP9(2,:) = ERP9(2,:) + T7R(i:i+127);
            M9 = M9 + 1;
       end
   elseif MARKERR(i) == 10 % SONG 1
       if trashcheck(T7R(i:i+255))
            ERP10(2,:) = ERP10(2,:) + T7R(i:i+255);
            M10 = M10 + 1;
       end
   elseif MARKERR(i) == 11 % SONG 2
       if trashcheck(T7R(i:i+255))
            ERP11(2,:) = ERP11(2,:) + T7R(i:i+255);
            M11 = M11 + 1;
       end
   elseif MARKERR(i) == 12 % SONG 3
       if trashcheck(T7R(i:i+255))
            ERP12(2,:) = ERP12(2,:) + T7R(i:i+255);
            M12 = M12 + 1;
       end
   end
end
ERP5(2,:) = ERP5(2,:)/M5;
ERP6(2,:) = ERP6(2,:)/M6;
ERP7(2,:) = ERP7(2,:)/M7;
ERP8(2,:) = ERP8(2,:)/M8;
ERP9(2,:) = ERP9(2,:)/M9;
ERP10(2,:) = ERP10(2,:)/M10;
ERP11(2,:) = ERP11(2,:)/M11;
ERP12(2,:) = ERP12(2,:)/M12;
%PzR
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERR(i) == 5 % SMELL
       if trashcheck(PzR(i:i+127))
            ERP5(3,:) = ERP5(3,:) + PzR(i:i+127);
            M5 = M5 + 1;
       end
   elseif MARKERR(i) == 6 % FOOD
       if trashcheck(PzR(i:i+127))
            ERP6(3,:) = ERP6(3,:) + PzR(i:i+127);
            M6 = M6 + 1;
       end
   elseif MARKERR(i) == 7 % TOOL
       if trashcheck(PzR(i:i+127))
            ERP7(3,:) = ERP7(3,:) + PzR(i:i+127);
            M7 = M7 + 1;
       end
   elseif MARKERR(i) == 8 % GESTURE
       if trashcheck(PzR(i:i+127))
            ERP8(3,:) = ERP8(3,:) + PzR(i:i+127);
            M8 = M8 + 1;
       end
   elseif MARKERR(i) == 9 % PERSON
       if trashcheck(PzR(i:i+127))
            ERP9(3,:) = ERP9(3,:) + PzR(i:i+127);
            M9 = M9 + 1;
       end
   elseif MARKERR(i) == 10 % SONG 1
       if trashcheck(PzR(i:i+255))
            ERP10(3,:) = ERP10(3,:) + PzR(i:i+255);
            M10 = M10 + 1;
       end
   elseif MARKERR(i) == 11 % SONG 2
       if trashcheck(PzR(i:i+255))
            ERP11(3,:) = ERP11(3,:) + PzR(i:i+255);
            M11 = M11 + 1;
       end
   elseif MARKERR(i) == 12 % SONG 3
       if trashcheck(PzR(i:i+255))
            ERP12(3,:) = ERP12(3,:) + PzR(i:i+255);
            M12 = M12 + 1;
       end
   end
end
ERP5(3,:) = ERP5(3,:)/M5;
ERP6(3,:) = ERP6(3,:)/M6;
ERP7(3,:) = ERP7(3,:)/M7;
ERP8(3,:) = ERP8(3,:)/M8;
ERP9(3,:) = ERP9(3,:)/M9;
ERP10(3,:) = ERP10(3,:)/M10;
ERP11(3,:) = ERP11(3,:)/M11;
ERP12(3,:) = ERP12(3,:)/M12;
%T8R
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERR(i) == 5 % SMELL
       if trashcheck(T8R(i:i+127))
            ERP5(4,:) = ERP5(4,:) + T8R(i:i+127);
            M5 = M5 + 1;
       end
   elseif MARKERR(i) == 6 % FOOD
       if trashcheck(T8R(i:i+127))
            ERP6(4,:) = ERP6(4,:) + T8R(i:i+127);
            M6 = M6 + 1;
       end
   elseif MARKERR(i) == 7 % TOOL
       if trashcheck(T8R(i:i+127))
            ERP7(4,:) = ERP7(4,:) + T8R(i:i+127);
            M7 = M7 + 1;
       end
   elseif MARKERR(i) == 8 % GESTURE
       if trashcheck(T8R(i:i+127))
            ERP8(4,:) = ERP8(4,:) + T8R(i:i+127);
            M8 = M8 + 1;
       end
   elseif MARKERR(i) == 9 % PERSON
       if trashcheck(T8R(i:i+127))
            ERP9(4,:) = ERP9(4,:) + T8R(i:i+127);
            M9 = M9 + 1;
       end
   elseif MARKERR(i) == 10 % SONG 1
       if trashcheck(T8R(i:i+255))
            ERP10(4,:) = ERP10(4,:) + T8R(i:i+255);
            M10 = M10 + 1;
       end
   elseif MARKERR(i) == 11 % SONG 2
       if trashcheck(T8R(i:i+255))
            ERP11(4,:) = ERP11(4,:) + T8R(i:i+255);
            M11 = M11 + 1;
       end
   elseif MARKERR(i) == 12 % SONG 3
       if trashcheck(T8R(i:i+255))
            ERP12(4,:) = ERP12(4,:) + T8R(i:i+255);
            M12 = M12 + 1;
       end
   end
end
ERP5(4,:) = ERP5(4,:)/M5;
ERP6(4,:) = ERP6(4,:)/M6;
ERP7(4,:) = ERP7(4,:)/M7;
ERP8(4,:) = ERP8(4,:)/M8;
ERP9(4,:) = ERP9(4,:)/M9;
ERP10(4,:) = ERP10(4,:)/M10;
ERP11(4,:) = ERP11(4,:)/M11;
ERP12(4,:) = ERP12(4,:)/M12;
%AF4R
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERR(i) == 5 % SMELL
       if trashcheck(AF4R(i:i+127))
            ERP5(5,:) = ERP5(5,:) + AF4R(i:i+127);
            M5 = M5 + 1;
       end
   elseif MARKERR(i) == 6 % FOOD
       if trashcheck(AF4R(i:i+127))
            ERP6(5,:) = ERP6(5,:) + AF4R(i:i+127);
            M6 = M6 + 1;
       end
   elseif MARKERR(i) == 7 % TOOL
       if trashcheck(AF4R(i:i+127))
            ERP7(5,:) = ERP7(5,:) + AF4R(i:i+127);
            M7 = M7 + 1;
       end
   elseif MARKERR(i) == 8 % GESTURE
       if trashcheck(AF4R(i:i+127))
            ERP8(5,:) = ERP8(5,:) + AF4R(i:i+127);
            M8 = M8 + 1;
       end
   elseif MARKERR(i) == 9 % PERSON
       if trashcheck(AF4R(i:i+127))
            ERP9(5,:) = ERP9(5,:) + AF4R(i:i+127);
            M9 = M9 + 1;
       end
   elseif MARKERR(i) == 10 % SONG 1
       if trashcheck(AF4R(i:i+255))
            ERP10(5,:) = ERP10(5,:) + AF4R(i:i+255);
            M10 = M10 + 1;
       end
   elseif MARKERR(i) == 11 % SONG 2
       if trashcheck(AF4R(i:i+255))
            ERP11(5,:) = ERP11(5,:) + AF4R(i:i+255);
            M11 = M11 + 1;
       end
   elseif MARKERR(i) == 12 % SONG 3
       if trashcheck(AF4R(i:i+255))
            ERP12(5,:) = ERP12(5,:) + AF4R(i:i+255);
            M12 = M12 + 1;
       end
   end
end
ERP5(5,:) = ERP5(5,:)/M5;
ERP6(5,:) = ERP6(5,:)/M6;
ERP7(5,:) = ERP7(5,:)/M7;
ERP8(5,:) = ERP8(5,:)/M8;
ERP9(5,:) = ERP9(5,:)/M9;
ERP10(5,:) = ERP10(5,:)/M10;
ERP11(5,:) = ERP11(5,:)/M11;
ERP12(5,:) = ERP12(5,:)/M12;

ERP5 = scaler(ERP5);
ERP6 = scaler(ERP6);
ERP7 = scaler(ERP7);
ERP8 = scaler(ERP8);
ERP9 = scaler(ERP9);
ERP10 = scaler(ERP10);
ERP11 = scaler(ERP11);
ERP12 = scaler(ERP12);

% Find autocorrelation maximum for each ERP to be used for normalizing the 
% range to [-1, 1] when cross-correlating with challenge ERPs
ERP5ACM = dot(ERP5,ERP5,2);
ERP6ACM = dot(ERP6,ERP6,2);
ERP7ACM = dot(ERP7,ERP7,2);
ERP8ACM = dot(ERP8,ERP8,2);
ERP9ACM = dot(ERP9,ERP9,2);
ERP10ACM = dot(ERP10,ERP10,2);
ERP11ACM = dot(ERP11,ERP11,2);
ERP12ACM = dot(ERP12,ERP12,2);

%--------------------------------------------------------------------------
% Start analysis of challenge ERPs
% These vectors keep track of dot products of the challenge ERP to the
% reference ERP as the number of trials used to calculate the challenge
% ERP increases 
ERP5COMP = zeros(5,10);
ERP6COMP = zeros(5,10);
ERP7COMP = zeros(5,10);
ERP8COMP = zeros(5,10);
ERP9COMP = zeros(5,10);
ERP10COMP = zeros(5,15);
ERP11COMP = zeros(5,15);
ERP12COMP = zeros(5,15);

% Initialize challenge ERPs for each of the stimuli
% Rows are: AF3, T7, Pz, T8, AF4
% Visual stimuli have 1000 ms ERP
ERP5C = zeros(5,Fs);
ERP6C = zeros(5,Fs); 
ERP7C = zeros(5,Fs);
ERP8C = zeros(5,Fs);
ERP9C = zeros(5,Fs);
% Audio stimuli have 2000 ms ERP
ERP10C = zeros(5,Fs*2);
ERP11C = zeros(5,Fs*2);
ERP12C = zeros(5,Fs*2);

% Challenge ERP calculations
[~, length] = size(AF3C);
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
% AF3C
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       if trashcheck(AF3C(i:i+127))
            ERP5C(1,:) = ERP5C(1,:) + AF3C(i:i+127);
            M5 = M5 + 1;
            tempERP = ERP5C(1,:)/M5;
            tempdot = dot(ERP5(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP5ACM(1);
            ERP5COMP(1,M5) = tempnorm;
       end
   elseif MARKERC(i) == 6 % FOOD
       if trashcheck(AF3C(i:i+127))
            ERP6C(1,:) = ERP6C(1,:) + AF3C(i:i+127);
            M6 = M6 + 1;
            tempERP = ERP6C(1,:)/M6;
            tempdot = dot(ERP6(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP6ACM(1);
            ERP6COMP(1,M6) = tempnorm;
       end
   elseif MARKERC(i) == 7 % TOOL
       if trashcheck(AF3C(i:i+127))
            ERP7C(1,:) = ERP7C(1,:) + AF3C(i:i+127);
            M7 = M7 + 1;
            tempERP = ERP7C(1,:)/M7;
            tempdot = dot(ERP7(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP7ACM(1);
            ERP7COMP(1,M7) = tempnorm;
       end
   elseif MARKERC(i) == 8 % GESTURE
       if trashcheck(AF3C(i:i+127))
            ERP8C(1,:) = ERP8C(1,:) + AF3C(i:i+127);
            M8 = M8 + 1;
            tempERP = ERP8C(1,:)/M8;
            tempdot = dot(ERP8(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP8ACM(1);
            ERP8COMP(1,M8) = tempnorm;
       end
   elseif MARKERC(i) == 9 % PERSON
       if trashcheck(AF3C(i:i+127))
            ERP9C(1,:) = ERP9C(1,:) + AF3C(i:i+127);
            M9 = M9 + 1;
            tempERP = ERP9C(1,:)/M9;
            tempdot = dot(ERP9(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP9ACM(1);
            ERP9COMP(1,M9) = tempnorm;
       end
   elseif MARKERC(i) == 10 % SONG 1
       if trashcheck(AF3C(i:i+255))
            ERP10C(1,:) = ERP10C(1,:) + AF3C(i:i+255);
            M10 = M10 + 1;
            tempERP = ERP10C(1,:)/M10;
            tempdot = dot(ERP10(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP10ACM(1);
            ERP10COMP(1,M10) = tempnorm;
       end
   elseif MARKERC(i) == 11 % SONG 2
       if trashcheck(AF3C(i:i+255))
            ERP11C(1,:) = ERP11C(1,:) + AF3C(i:i+255);
            M11 = M11 + 1;
            tempERP = ERP11C(1,:)/M11;
            tempdot = dot(ERP11(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP11ACM(1);
            ERP11COMP(1,M11) = tempnorm;
       end
   elseif MARKERC(i) == 12 % SONG 3
       if trashcheck(AF3C(i:i+255))
            ERP12C(1,:) = ERP12C(1,:) + AF3C(i:i+255);
            M12 = M12 + 1;
            tempERP = ERP12C(1,:)/M12;
            tempdot = dot(ERP12(1,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP12ACM(1);
            ERP12COMP(1,M12) = tempnorm;
       end
   end
end
ERP5C(1,:) = ERP5C(1,:)/M5;
ERP6C(1,:) = ERP6C(1,:)/M6;
ERP7C(1,:) = ERP7C(1,:)/M7;
ERP8C(1,:) = ERP8C(1,:)/M8;
ERP9C(1,:) = ERP9C(1,:)/M9;
ERP10C(1,:) = ERP10C(1,:)/M10;
ERP11C(1,:) = ERP11C(1,:)/M11;
ERP12C(1,:) = ERP12C(1,:)/M12;
% T7C
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       if trashcheck(T7C(i:i+127))
            ERP5C(2,:) = ERP5C(2,:) + T7C(i:i+127);
            M5 = M5 + 1;
            tempERP = ERP5C(2,:)/M5;
            tempdot = dot(ERP5(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP5ACM(2);
            ERP5COMP(2,M5) = tempnorm;
       end
   elseif MARKERC(i) == 6 % FOOD
       if trashcheck(T7C(i:i+127))
            ERP6C(2,:) = ERP6C(2,:) + T7C(i:i+127);
            M6 = M6 + 1;
            tempERP = ERP6C(2,:)/M6;
            tempdot = dot(ERP6(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP6ACM(2);
            ERP6COMP(2,M6) = tempnorm;
       end
   elseif MARKERC(i) == 7 % TOOL
       if trashcheck(T7C(i:i+127))
            ERP7C(2,:) = ERP7C(2,:) + T7C(i:i+127);
            M7 = M7 + 1;
            tempERP = ERP7C(2,:)/M7;
            tempdot = dot(ERP7(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP7ACM(2);
            ERP7COMP(2,M7) = tempnorm;
       end
   elseif MARKERC(i) == 8 % GESTURE
       if trashcheck(T7C(i:i+127))
            ERP8C(2,:) = ERP8C(2,:) + T7C(i:i+127);
            M8 = M8 + 1;
            tempERP = ERP8C(2,:)/M8;
            tempdot = dot(ERP8(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP8ACM(2);
            ERP8COMP(2,M8) = tempnorm;
       end
   elseif MARKERC(i) == 9 % PERSON
       if trashcheck(T7C(i:i+127))
            ERP9C(2,:) = ERP9C(2,:) + T7C(i:i+127);
            M9 = M9 + 1;
            tempERP = ERP9C(2,:)/M9;
            tempdot = dot(ERP9(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP9ACM(2);
            ERP9COMP(2,M9) = tempnorm;
       end
   elseif MARKERC(i) == 10 % SONG 1
       if trashcheck(T7C(i:i+255))
            ERP10C(2,:) = ERP10C(2,:) + T7C(i:i+255);
            M10 = M10 + 1;
            tempERP = ERP10C(2,:)/M10;
            tempdot = dot(ERP10(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP10ACM(2);
            ERP10COMP(2,M10) = tempnorm;
       end
   elseif MARKERC(i) == 11 % SONG 2
       if trashcheck(T7C(i:i+255))
            ERP11C(2,:) = ERP11C(2,:) + T7C(i:i+255);
            M11 = M11 + 1;
            tempERP = ERP11C(2,:)/M11;
            tempdot = dot(ERP11(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP11ACM(2);
            ERP11COMP(2,M11) = tempnorm;
       end
   elseif MARKERC(i) == 12 % SONG 3
       if trashcheck(T7C(i:i+255))
            ERP12C(2,:) = ERP12C(2,:) + T7C(i:i+255);
            M12 = M12 + 1;
            tempERP = ERP12C(2,:)/M12;
            tempdot = dot(ERP12(2,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP12ACM(2);
            ERP12COMP(2,M12) = tempnorm;
       end
   end
end
ERP5C(2,:) = ERP5C(2,:)/M5;
ERP6C(2,:) = ERP6C(2,:)/M6;
ERP7C(2,:) = ERP7C(2,:)/M7;
ERP8C(2,:) = ERP8C(2,:)/M8;
ERP9C(2,:) = ERP9C(2,:)/M9;
ERP10C(2,:) = ERP10C(2,:)/M10;
ERP11C(2,:) = ERP11C(2,:)/M11;
ERP12C(2,:) = ERP12C(2,:)/M12;
%PzC
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       if trashcheck(PzC(i:i+127))
            ERP5C(3,:) = ERP5C(3,:) + PzC(i:i+127);
            M5 = M5 + 1;
            tempERP = ERP5C(3,:)/M5;
            tempdot = dot(ERP5(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP5ACM(3);
            ERP5COMP(3,M5) = tempnorm;
       end
   elseif MARKERC(i) == 6 % FOOD
       if trashcheck(PzC(i:i+127))
            ERP6C(3,:) = ERP6C(3,:) + PzC(i:i+127);
            M6 = M6 + 1;
            tempERP = ERP6C(3,:)/M6;
            tempdot = dot(ERP6(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP6ACM(3);
            ERP6COMP(3,M6) = tempnorm;
       end
   elseif MARKERC(i) == 7 % TOOL
       if trashcheck(PzC(i:i+127))
            ERP7C(3,:) = ERP7C(3,:) + PzC(i:i+127);
            M7 = M7 + 1;
            tempERP = ERP7C(3,:)/M7;
            tempdot = dot(ERP7(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP7ACM(3);
            ERP7COMP(3,M7) = tempnorm;
       end
   elseif MARKERC(i) == 8 % GESTURE
       if trashcheck(PzC(i:i+127))
            ERP8C(3,:) = ERP8C(3,:) + PzC(i:i+127);
            M8 = M8 + 1;
            tempERP = ERP8C(3,:)/M8;
            tempdot = dot(ERP8(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP8ACM(3);
            ERP8COMP(3,M8) = tempnorm;
       end
   elseif MARKERC(i) == 9 % PERSON
       if trashcheck(PzC(i:i+127))
            ERP9C(3,:) = ERP9C(3,:) + PzC(i:i+127);
            M9 = M9 + 1;
            tempERP = ERP9C(3,:)/M9;
            tempdot = dot(ERP9(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP9ACM(3);
            ERP9COMP(3,M9) = tempnorm;
       end
   elseif MARKERC(i) == 10 % SONG 1
       if trashcheck(PzC(i:i+255))
            ERP10C(3,:) = ERP10C(3,:) + PzC(i:i+255);
            M10 = M10 + 1;
            tempERP = ERP10C(3,:)/M10;
            tempdot = dot(ERP10(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP10ACM(3);
            ERP10COMP(3,M10) = tempnorm;
       end
   elseif MARKERC(i) == 11 % SONG 2
       if trashcheck(PzC(i:i+255))
            ERP11C(3,:) = ERP11C(3,:) + PzC(i:i+255);
            M11 = M11 + 1;
            tempERP = ERP11C(3,:)/M11;
            tempdot = dot(ERP11(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP11ACM(3);
            ERP11COMP(3,M11) = tempnorm;
       end
   elseif MARKERC(i) == 12 % SONG 3
       if trashcheck(PzC(i:i+255))
            ERP12C(3,:) = ERP12C(3,:) + PzC(i:i+255);
            M12 = M12 + 1;
            tempERP = ERP12C(3,:)/M12;
            tempdot = dot(ERP12(3,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP12ACM(3);
            ERP12COMP(3,M12) = tempnorm;
       end
   end
end
ERP5C(3,:) = ERP5C(3,:)/M5;
ERP6C(3,:) = ERP6C(3,:)/M6;
ERP7C(3,:) = ERP7C(3,:)/M7;
ERP8C(3,:) = ERP8C(3,:)/M8;
ERP9C(3,:) = ERP9C(3,:)/M9;
ERP10C(3,:) = ERP10C(3,:)/M10;
ERP11C(3,:) = ERP11C(3,:)/M11;
ERP12C(3,:) = ERP12C(3,:)/M12;
%T8C
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       if trashcheck(T8C(i:i+127))
            ERP5C(4,:) = ERP5C(4,:) + T8C(i:i+127);
            M5 = M5 + 1;
            tempERP = ERP5C(4,:)/M5;
            tempdot = dot(ERP5(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP5ACM(4);
            ERP5COMP(4,M5) = tempnorm;
       end
   elseif MARKERC(i) == 6 % FOOD
       if trashcheck(T8C(i:i+127))
            ERP6C(4,:) = ERP6C(4,:) + T8C(i:i+127);
            M6 = M6 + 1;
            tempERP = ERP6C(4,:)/M6;
            tempdot = dot(ERP6(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP6ACM(4);
            ERP6COMP(4,M6) = tempnorm;
       end
   elseif MARKERC(i) == 7 % TOOL
       if trashcheck(T8C(i:i+127))
            ERP7C(4,:) = ERP7C(4,:) + T8C(i:i+127);
            M7 = M7 + 1;
            tempERP = ERP7C(4,:)/M7;
            tempdot = dot(ERP7(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP7ACM(4);
            ERP7COMP(4,M7) = tempnorm;
       end
   elseif MARKERC(i) == 8 % GESTURE
       if trashcheck(T8C(i:i+127))
            ERP8C(4,:) = ERP8C(4,:) + T8C(i:i+127);
            M8 = M8 + 1;
            tempERP = ERP8C(4,:)/M8;
            tempdot = dot(ERP8(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP8ACM(4);
            ERP8COMP(4,M8) = tempnorm;
       end
   elseif MARKERC(i) == 9 % PERSON
       if trashcheck(T8C(i:i+127))
            ERP9C(4,:) = ERP9C(4,:) + T8C(i:i+127);
            M9 = M9 + 1;
            tempERP = ERP9C(4,:)/M9;
            tempdot = dot(ERP9(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP9ACM(4);
            ERP9COMP(4,M9) = tempnorm;
       end
   elseif MARKERC(i) == 10 % SONG 1
       if trashcheck(T8C(i:i+255))
            ERP10C(4,:) = ERP10C(4,:) + T8C(i:i+255);
            M10 = M10 + 1;
            tempERP = ERP10C(4,:)/M10;
            tempdot = dot(ERP10(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP10ACM(4);
            ERP10COMP(4,M10) = tempnorm;
       end
   elseif MARKERC(i) == 11 % SONG 2
       if trashcheck(T8C(i:i+255))
            ERP11C(4,:) = ERP11C(4,:) + T8C(i:i+255);
            M11 = M11 + 1;
            tempERP = ERP11C(4,:)/M11;
            tempdot = dot(ERP11(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP11ACM(4);
            ERP11COMP(4,M11) = tempnorm;
       end
   elseif MARKERC(i) == 12 % SONG 3
       if trashcheck(T8C(i:i+255))
            ERP12C(4,:) = ERP12C(4,:) + T8C(i:i+255);
            M12 = M12 + 1;
            tempERP = ERP12C(4,:)/M12;
            tempdot = dot(ERP12(4,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP12ACM(4);
            ERP12COMP(4,M12) = tempnorm;
       end
   end
end
ERP5C(4,:) = ERP5C(4,:)/M5;
ERP6C(4,:) = ERP6C(4,:)/M6;
ERP7C(4,:) = ERP7C(4,:)/M7;
ERP8C(4,:) = ERP8C(4,:)/M8;
ERP9C(4,:) = ERP9C(4,:)/M9;
ERP10C(4,:) = ERP10C(4,:)/M10;
ERP11C(4,:) = ERP11C(4,:)/M11;
ERP12C(4,:) = ERP12C(4,:)/M12;
%AF4C
M5 = 0;
M6 = 0;
M7 = 0;
M8 = 0;
M9 = 0;
M10 = 0;
M11 = 0;
M12 = 0;
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       if trashcheck(AF4C(i:i+127))
            ERP5C(5,:) = ERP5C(5,:) + AF4C(i:i+127);
            M5 = M5 + 1;
            tempERP = ERP5C(5,:)/M5;
            tempdot = dot(ERP5(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP5ACM(5);
            ERP5COMP(5,M5) = tempnorm;
       end
   elseif MARKERC(i) == 6 % FOOD
       if trashcheck(AF4C(i:i+127))
            ERP6C(5,:) = ERP6C(5,:) + AF4C(i:i+127);
            M6 = M6 + 1;
            tempERP = ERP6C(5,:)/M6;
            tempdot = dot(ERP6(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP6ACM(5);
            ERP6COMP(5,M6) = tempnorm;
       end
   elseif MARKERC(i) == 7 % TOOL
       if trashcheck(AF4C(i:i+127))
            ERP7C(5,:) = ERP7C(5,:) + AF4C(i:i+127);
            M7 = M7 + 1;
            tempERP = ERP7C(5,:)/M7;
            tempdot = dot(ERP7(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP7ACM(5);
            ERP7COMP(5,M7) = tempnorm;
       end
   elseif MARKERC(i) == 8 % GESTURE
       if trashcheck(AF4C(i:i+127))
            ERP8C(5,:) = ERP8C(5,:) + AF4C(i:i+127);
            M8 = M8 + 1;
            tempERP = ERP8C(5,:)/M8;
            tempdot = dot(ERP8(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP8ACM(5);
            ERP8COMP(5,M8) = tempnorm;
       end
   elseif MARKERC(i) == 9 % PERSON
       if trashcheck(AF4C(i:i+127))
            ERP9C(5,:) = ERP9C(5,:) + AF4C(i:i+127);
            M9 = M9 + 1;
            tempERP = ERP9C(5,:)/M9;
            tempdot = dot(ERP9(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP9ACM(5);
            ERP9COMP(5,M9) = tempnorm;
       end
   elseif MARKERC(i) == 10 % SONG 1
       if trashcheck(AF4C(i:i+255))
            ERP10C(5,:) = ERP10C(5,:) + AF4C(i:i+255);
            M10 = M10 + 1;
            tempERP = ERP10C(5,:)/M10;
            tempdot = dot(ERP10(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP10ACM(5);
            ERP10COMP(5,M10) = tempnorm;
       end
   elseif MARKERC(i) == 11 % SONG 2
       if trashcheck(AF4C(i:i+255))
            ERP11C(5,:) = ERP11C(5,:) + AF4C(i:i+255);
            M11 = M11 + 1;
            tempERP = ERP11C(5,:)/M11;
            tempdot = dot(ERP11(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP11ACM(5);
            ERP11COMP(5,M11) = tempnorm;
       end
   elseif MARKERC(i) == 12 % SONG 3
       if trashcheck(AF4C(i:i+255))
            ERP12C(5,:) = ERP12C(5,:) + AF4C(i:i+255);
            M12 = M12 + 1;
            tempERP = ERP12C(5,:)/M12;
            tempdot = dot(ERP12(5,:),sscaler(tempERP),2);
            tempnorm = tempdot./ERP12ACM(5);
            ERP12COMP(5,M12) = tempnorm;
       end
   end
end
ERP5C(5,:) = ERP5C(5,:)/M5;
ERP6C(5,:) = ERP6C(5,:)/M6;
ERP7C(5,:) = ERP7C(5,:)/M7;
ERP8C(5,:) = ERP8C(5,:)/M8;
ERP9C(5,:) = ERP9C(5,:)/M9;
ERP10C(5,:) = ERP10C(5,:)/M10;
ERP11C(5,:) = ERP11C(5,:)/M11;
ERP12C(5,:) = ERP12C(5,:)/M12;
%{
% Challenge ERP derivation
[~, length] = size(AF3C);
N5 = 0;
N6 = 0;
N7 = 0;
N8 = 0;
N9 = 0;
N10 = 0;
N11 = 0;
N12 = 0;
for i = 1:length
   if MARKERC(i) == 5 % SMELL
       ERP5C(1,:) = ERP5C(1,:) + AF3C(i:i+127);
       ERP5C(2,:) = ERP5C(2,:) + T7C(i:i+127);
       ERP5C(3,:) = ERP5C(3,:) + PzC(i:i+127);
       ERP5C(4,:) = ERP5C(4,:) + T8C(i:i+127);
       ERP5C(5,:) = ERP5C(5,:) + AF4C(i:i+127);
       N5 = N5 + 1;
       tempERP = [ERP5C(1,:);
           ERP5C(2,:);
           ERP5C(3,:);
           ERP5C(4,:);
           ERP5C(5,:)]/N5;
       tempdot = dot(ERP5,scaler(tempERP),2);
       tempnorm = tempdot./ERP5ACM;
       ERP5COMP(:,N5) = tempnorm;
   elseif MARKERC(i) == 6 % FOOD
       ERP6C(1,:) = ERP6C(1,:) + AF3C(i:i+127);
       ERP6C(2,:) = ERP6C(2,:) + T7C(i:i+127);
       ERP6C(3,:) = ERP6C(3,:) + PzC(i:i+127);
       ERP6C(4,:) = ERP6C(4,:) + T8C(i:i+127);
       ERP6C(5,:) = ERP6C(5,:) + AF4C(i:i+127);
       N6 = N6 + 1;
       tempERP = [ERP6C(1,:);
           ERP6C(2,:);
           ERP6C(3,:);
           ERP6C(4,:);
           ERP6C(5,:)]/N6;
       tempdot = dot(ERP6,scaler(tempERP),2);
       tempnorm = tempdot./ERP6ACM;
       ERP6COMP(:,N6) = tempnorm;
   elseif MARKERC(i) == 7 % TOOL
       ERP7C(1,:) = ERP7C(1,:) + AF3C(i:i+127);
       ERP7C(2,:) = ERP7C(2,:) + T7C(i:i+127);
       ERP7C(3,:) = ERP7C(3,:) + PzC(i:i+127);
       ERP7C(4,:) = ERP7C(4,:) + T8C(i:i+127);
       ERP7C(5,:) = ERP7C(5,:) + AF4C(i:i+127);
       N7 = N7 + 1;
       tempERP = [ERP7C(1,:);
           ERP7C(2,:);
           ERP7C(3,:);
           ERP7C(4,:);
           ERP7C(5,:)]/N7;
       tempdot = dot(ERP7,scaler(tempERP),2);
       tempnorm = tempdot./ERP7ACM;
       ERP7COMP(:,N7) = tempnorm;
   elseif MARKERC(i) == 8 % GESTURE
       ERP8C(1,:) = ERP8C(1,:) + AF3C(i:i+127);
       ERP8C(2,:) = ERP8C(2,:) + T7C(i:i+127);
       ERP8C(3,:) = ERP8C(3,:) + PzC(i:i+127);
       ERP8C(4,:) = ERP8C(4,:) + T8C(i:i+127);
       ERP8C(5,:) = ERP8C(5,:) + AF4C(i:i+127);
       N8 = N8 + 1;
       tempERP = [ERP8C(1,:);
           ERP8C(2,:);
           ERP8C(3,:);
           ERP8C(4,:);
           ERP8C(5,:)]/N8;
       tempdot = dot(ERP8,scaler(tempERP),2);
       tempnorm = tempdot./ERP8ACM;
       ERP8COMP(:,N8) = tempnorm;
   elseif MARKERC(i) == 9 % PERSON
       ERP9C(1,:) = ERP9C(1,:) + AF3C(i:i+127);
       ERP9C(2,:) = ERP9C(2,:) + T7C(i:i+127);
       ERP9C(3,:) = ERP9C(3,:) + PzC(i:i+127);
       ERP9C(4,:) = ERP9C(4,:) + T8C(i:i+127);
       ERP9C(5,:) = ERP9C(5,:) + AF4C(i:i+127);
       N9 = N9 + 1;
       tempERP = [ERP9C(1,:);
           ERP9C(2,:);
           ERP9C(3,:);
           ERP9C(4,:);
           ERP9C(5,:)]/N9;
       tempdot = dot(ERP9,scaler(tempERP),2);
       tempnorm = tempdot./ERP9ACM;
       ERP9COMP(:,N9) = tempnorm;
   elseif MARKERC(i) == 10 % SONG 1
       ERP10C(1,:) = ERP10C(1,:) + AF3C(i:i+255);
       ERP10C(2,:) = ERP10C(2,:) + T7C(i:i+255);
       ERP10C(3,:) = ERP10C(3,:) + PzC(i:i+255);
       ERP10C(4,:) = ERP10C(4,:) + T8C(i:i+255);
       ERP10C(5,:) = ERP10C(5,:) + AF4C(i:i+255);
       N10 = N10 + 1;
       tempERP = [ERP10C(1,:);
           ERP10C(2,:);
           ERP10C(3,:);
           ERP10C(4,:);
           ERP10C(5,:)]/N10;
       tempdot = dot(ERP10,scaler(tempERP),2);
       tempnorm = tempdot./ERP10ACM;
       ERP10COMP(:,N10) = tempnorm;
   elseif MARKERC(i) == 11 % SONG 2
       ERP11C(1,:) = ERP11C(1,:) + AF3C(i:i+255);
       ERP11C(2,:) = ERP11C(2,:) + T7C(i:i+255);
       ERP11C(3,:) = ERP11C(3,:) + PzC(i:i+255);
       ERP11C(4,:) = ERP11C(4,:) + T8C(i:i+255);
       ERP11C(5,:) = ERP11C(5,:) + AF4C(i:i+255);
       N11 = N11 + 1;
       tempERP = [ERP11C(1,:);
           ERP11C(2,:);
           ERP11C(3,:);
           ERP11C(4,:);
           ERP11C(5,:)]/N11;
       %A = corrcoef(ERP11,tempERP);
       %testcore = [testcore A(1,2)];
       tempdot = dot(ERP11,scaler(tempERP),2);
       tempnorm = tempdot./ERP11ACM;
       ERP11COMP(:,N11) = tempnorm;
   elseif MARKERC(i) == 12 % SONG 3
       ERP12C(1,:) = ERP12C(1,:) + AF3C(i:i+255);
       ERP12C(2,:) = ERP12C(2,:) + T7C(i:i+255);
       ERP12C(3,:) = ERP12C(3,:) + PzC(i:i+255);
       ERP12C(4,:) = ERP12C(4,:) + T8C(i:i+255);
       ERP12C(5,:) = ERP12C(5,:) + AF4C(i:i+255);
       N12 = N12 + 1;
       tempERP = [ERP12C(1,:);
           ERP12C(2,:);
           ERP12C(3,:);
           ERP12C(4,:);
           ERP12C(5,:)]/N12;
       x = scaler(tempERP);
       tempdot = dot(ERP12,scaler(tempERP),2);
       tempnorm = tempdot./ERP12ACM;
       ERP12COMP(:,N12) = tempnorm;
   end
end
%}
%Plot
CHANZ = {'AF3','T7','Pz','T8','AF4'};
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP5COMP(i,:))
    title(sprintf('%s - SMELL ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP6COMP(i,:))
    title(sprintf('%s - FOOD ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP7COMP(i,:))
    title(sprintf('%s - MOVE ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP8COMP(i,:))
    title(sprintf('%s - GESTURE ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP9COMP(i,:))
    title(sprintf('%s - PERSON ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP10COMP(i,:))
    title(sprintf('%s - SONG 1 ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP11COMP(i,:))
    title(sprintf('%s - SONG 2 ERP',CHANZ{i}))
end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:5
    subplot(5,1,i)
    plot(ERP12COMP(i,:))
    title(sprintf('%s - SONG 3 ERP',CHANZ{i}))
end

% Close
fclose(IDmain);
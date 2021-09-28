clear sound
close all
clear
filename = input('What''s the filename? ','s');
%songnum = input('Which song is this? ','s');
songnum = filename(1);
pause('on');
n = 60; % number of samples to make

[audio, Fs] = audioread(filename);
length = length(audio);
maxlength = length - Fs - 1;
for i = 1:n
   a = rand;
   audioloc = maxlength*a;
   audioloc = round(audioloc);
   tempaudio = audio(audioloc:audioloc+Fs); %1 sec sample
   
   rng = 0.5*rand()+2.5;   % Between 2.5-3 seconds of silence
   zeromat = zeros(1,round(Fs*rng));

   tempaudio = [tempaudio zeromat];
   tempname = sprintf('%s_%d.wav',songnum,i);
   audiowrite(tempname,tempaudio,Fs)
   %soundsc(tempaudio,Fs) % COMMENT OUT WHEN ACTUALLY RUNNING
   
   %pause(4); % COMMENT OUT WHEN ACTUALLY RUNNING
end
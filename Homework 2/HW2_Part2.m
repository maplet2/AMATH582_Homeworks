clear all; close all; clc

% In this code we will look at two different music files playing the song
% 'Mary had a Little Lamb' played on a piano and on a recorder. The goal is
% then to use the Gabor Transformation to create a spectrogram that would
% replicate the music sheet for the song and also highlight the difference
% between the song played on two different instruments.

% Due to memory limitations, we will only look at the first verse of the song
tr_piano=16 / 2; % total record time in seconds divided in half
py=audioread('music1.wav');
py = py';
% Take seconds 0-8 of the song
pstart = 1;
pstop = (length(py)/2);
py = py(1:pstop);
%%
Fs=length(py)/tr_piano;
plot((pstart:pstop)/Fs,py);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow
%p8 = audioplayer(yp,Fs); playblocking(p8);
pn = length(py);
%kp = (2*pi)*[0:np/2 -np/2:-1]; % if n is odd
pk = (2*pi)*[0:pn/2-1 -pn/2:-1]; % if n is even
%%
timeslide=0:0.05:tr_piano;
pt = linspace(0,4,pn);
p_spc = zeros(length(timeslide), length(py));
%%
width = 10000;
for n=1:length(timeslide)
    g = exp(-width*(pt-timeslide(n)).^10);
    pyf = py.*g;
    pyft = fft(pyf);
    p_spc(n,:) = abs(fftshift(pyft));
end

%%
figure(2)
clf;
pcolor(timeslide, fftshift(pk), p_spc.'), shading interp, colormap(hot)
ybounds = 5*10^4;
ylim([0 ybounds])
title('Piano Music Score')
xlabel('Time in Seconds')
ylabel('Frequency')

%%
figure(3)
clf;
tr_rec=4; % record time in seconds
yr=audioread('music2.wav'); 
yr = yr';
yr_end = length(yr)-8;
yr = yr(1:yr_end);
rstart = (length(yr)/14)*7;
rstop = (length(yr)/14)*11;
yr=yr(rstart:rstop);

Fs=length(yr)/tr_rec;
tr=(rstart:rstop)/Fs;
plot(tr,yr);
xlabel('Time [sec]'); 
ylabel('Amplitude');
xlim([min(tr) max(tr)])
title('Mary had a little lamb (recorder)');
%p8 = audioplayer(yr,Fs); playblocking(p8);

nr = length(yr);
kr = (2*pi)*[0:nr/2 -nr/2:-1];
%%
timeslide=0:0.05:tr_rec;
tr = linspace(0,4,nr);
spc_r = zeros(length(timeslide), length(tr));

%%
for n=1:length(timeslide)
    g = exp(-width*(tr-timeslide(n)).^10);
    yrf = yr.*g;
    yrft = fft(yrf);
    spc_r(n,:) = abs(fftshift(yrft));
end

%%
figure(4)
clf;
timeslide_act = 7:0.05:11;
pcolor(timeslide_act, fftshift(kr), spc_r.'), shading interp, colormap(hot)
ybounds = 5*10^4;
ylim([0 ybounds])
title('Recorder Music Score')
xlabel('Time in Seconds')
ylabel('Frequency')

clear all; close all; clc

% Goal of this code is to apply different types of Gabor Transformation to
% a song by Handel to see how different types of filters would affect the
% resulting spectrogram. 

load handel

% Setting up the Data
v = y'/2;
n = length(v); 
t = (1:length(v))/Fs; 
L = max(t);
k = (2*pi)*[0:n/2 -n/2:-1];

% %% Try Different types of filters:
% 
% % Guassian Filter
% midpoint = max(t) / 2;
% width = 1;
% g=exp(-width(1)*((t-midpoint).^2));
% 
% % Shannon Filter
% Fs = zeros(1,n);
% midind = round(n / 2);
% width_sh = 5000;
% Fs(midind-width_sh:1:midind+width_sh)=ones(1, 2*width_sh+1);
% 
% % Mexican Hat Filter
% sigma = 1;
% m = (2/(sqrt(3*sigma)*pi^(1/4))).*((1-((t-midpoint)/sigma).^2).*exp(-((t-midpoint).^2)/(2*sigma^2)));
% %%
% tslide = 0:0.1:10;
% indjump = ones(1, 101);
% for j=1:length(tslide)-1
%     indjump(j+1) = indjump(j) + 725; 
% end
% %%
% spc_guass = [];
% spc_shan = [];
% spc_mehat = [];
% for j=1:length(tslide)
%     ga=exp(-width*(t-tslide(j)).^10);
%     mehat = (2/(sqrt(3*sigma)*pi^(1/4))).*((1-((t-tslide(j))/sigma).^2).*exp(-((t-tslide(j)).^2)/(2*sigma^2)));
%     sh = zeros(1,n);
%     sh_start = indjump(j)-width_sh; 
%     if sh_start < 1
%         sh_start = 1;
%     end   
%     sh_end = indjump(j)+width_sh;
%     if sh_end > n 
%         sh_end = n;
%     end
%     current_shannon_width = sh_end-sh_start;
%     sh(sh_start:1:sh_end)=ones(1, current_shannon_width+1);
%      vfg = g.*v;
%      vftg = fft(vfg);
%      spc_guass = [spc_guass; abs(fftshift(vftg))];
%      
%      
%      vfsh = sh.*v;
%      vftsh = fft(vfsh);
%      spc_shan = [spc_shan; abs(fftshift(vftsh))];
%      
%      
%      vfm = mehat.*v;
%      vftm = fft(vfm);
%      spc_mehat = [spc_mehat; abs(fftshift(vftm))];
% end
% 
% 
% %% Plot Top Part of Figure
% set(gcf, 'Position',  [100, 100, 1200, 800])
% clf
% subplot(4,3,1)
% plot(t,v, 'b')
% hold on
% plot(t, g, 'r','linewidth',2)
% xlim([0, t(end)])
% title('Signal of Interest with Guassian Filter');
% 
% subplot(4,3,2)
% plot(t,v, 'b')
% hold on
% plot(t,v, 'b');
% plot(t, Fs, 'r', 'linewidth', 2);
% xlim([0, t(end)])
% title('Signal of Interest with Shannon Filter');
% 
% subplot(4,3,3)
% plot(t,v, 'b')
% hold on
% plot(t, m, 'r','linewidth',2)
% title('Signal of Interest with Mexican Hat Filter');
% xlim([0, t(end)])
% %% Plot Bottom Part of Figure
% 
% subplot(4,3, [4 7 10])
% pcolor(tslide, fftshift(k), spc_guass.'), shading interp, colormap(hot)
% title('Spectrogram of Signal')
% xlim([0 t(end)]);
% xlabel('Time')
% ylabel('Frequency')
% 
% subplot(4,3, [5 8 11])
% pcolor(tslide, fftshift(k), spc_shan.'), shading interp, colormap(hot)
% title('Spectrogram of Signal')
% xlim([0 t(end)]);
% xlabel('Time')
% ylabel('Frequency')
% 
% subplot(4,3, [6 9 12])
% pcolor(tslide, fftshift(k), spc_mehat.'), shading interp, colormap(hot)
% title('Spectrogram of Signal')
% xlim([0 t(end)]);
% xlabel('Time')
% ylabel('Frequency')



%% Different Time Slides
midpoint = max(t) / 2;
width = 1;
g=exp(-width(1)*((t-midpoint).^10));

%% Step = 0.1
step1 = 0.1;
tslide1 = 0:step1:10;
spc1 = zeros(length(tslide1), n);
width = 1;
for j=1:length(tslide1)
    g1=exp(-width*(t-tslide1(j)).^10);
    vf = g1.*v;
    vft = fft(vf);
    spc1(j,:) = abs(fftshift(vft));
end
%% Step2 = 1
step2 = 1;
tslide2 = 0:step2:10;
spc2 = zeros(length(tslide2), n);
width = 1;
for j=1:length(tslide2)
    g2=exp(-width*(t-tslide2(j)).^10);
    vf = g2.*v;
    vft = fft(vf);
    spc2(j,:) = abs(fftshift(vft));
end
%% Step3 = 0.01
step3 = 0.01;
tslide3 = 0:step3:10;
spc3 = zeros(length(tslide3),n);
for j=1:length(tslide3)
    g3=exp(-width*(t-tslide3(j)).^10);
    vf = g3.*v;
    vft = fft(vf);
    spc3(j,:) = abs(fftshift(vft));
end

%% Ploting Top Part
g1=exp(-width*((t-midpoint-step1).^10));
g2=exp(-width*((t-midpoint-step2).^10));
g3=exp(-width*((t-midpoint-step3).^10));
clf;
set(gcf, 'Position',  [100, 100, 1200, 800])
subplot(4,3,1)
plot(t,v, 'b')
hold on
plot(t, g, 'r','linewidth',2)
plot(t, g1, 'g','linewidth',2)
xlim([0, t(end)])
title('Signal of Interest with Base Time Slide; slide = 0.1');

subplot(4,3,2)
plot(t,v, 'b')
hold on
plot(t, g, 'r','linewidth',2)
plot(t, g3, 'g','linewidth',2)
xlim([0, t(end)])
title('Signal of Interest with Small Time Slide; slide = 0.01');

subplot(4,3,3)
plot(t,v, 'b')
hold on
plot(t, g, 'r','linewidth',2)
plot(t, g2, 'g','linewidth',2)
title('Signal of Interest with Large Time Slide; slide = 1');
xlim([0, t(end)])
%% Plotting Bottom Part
subplot(4,3, [4 7 10])
pcolor(tslide1, fftshift(k), spc1.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0, t(end)])
xlabel('Time')
ylabel('Frequency')


subplot(4,3, [5 8 11])
pcolor(tslide3, fftshift(k), spc3.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0, t(end)])
xlabel('Time')
ylabel('Frequency')


subplot(4,3, [6 9 12])
pcolor(tslide2, fftshift(k), spc2.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0, t(end)])
xlabel('Time')
ylabel('Frequency')



%% Different Window Widths
% Filtering with Different Widths
tslide = 0:0.1:10;
spc1 = [];
spc2 = [];
spc3 = [];
width = [1, 10000, 0.0001];
for j=1:length(tslide)
    g1=exp(-width(1)*(t-tslide(j)).^10);
    g2=exp(-width(2)*(t-tslide(j)).^10);
    g3=exp(-width(3)*(t-tslide(j)).^10);
    vf1 = g1.*v;
    vft1 = fft(vf1);
    vf2 = g2.*v;
    vft2 = fft(vf2);
    vf3 = g3.*v;
    vft3 = fft(vf3);
    spc1 = [spc1; abs(fftshift(vft1))];
    spc2 = [spc2; abs(fftshift(vft2))];
    spc3 = [spc3; abs(fftshift(vft3))];
end

midpoint = max(t) / 2;
g1=exp(-width(1)*((t-midpoint).^10));
g2=exp(-width(2)*((t-midpoint).^10));
g3=exp(-width(3)*((t-midpoint).^10));

%% Plotting Different Window Widths Top Part
clf
set(gcf, 'Position',  [100, 100, 1200, 800])
subplot(4,3,1)
hold on
plot(t,v, 'b')
plot(t, g1, 'r','linewidth',2)
title('Signal of Interest with Base Width Filter');
xlabel('Time [sec]');
xlim([0 t(end)]);
ylabel('Amplitude');

subplot(4,3,2)
hold on
plot(t,v, 'b')
plot(t, g2, 'r','linewidth',2)
title('Signal of Interest with Small Width Filter');
xlabel('Time [sec]');
ylabel('Amplitude');
xlim([0 t(end)]);

subplot(4,3,3)
hold on
plot(t,v, 'b')
plot(t, g3, 'r','linewidth',2)
title('Signal of Interest with Large Width Filter');
xlabel('Time [sec]');
ylabel('Amplitude');
xlim([0 t(end)]);

%% Bottom Part

subplot(4,3, [4 7 10])
pcolor(tslide, fftshift(k), spc1.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0 t(end)]);
xlabel('Time')
ylabel('Frequency')

subplot(4,3, [5 8 11])
pcolor(tslide, fftshift(k), spc2.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0 t(end)]);
xlabel('Time')
ylabel('Frequency')

subplot(4,3, [6 9 12])
pcolor(tslide, fftshift(k), spc3.'), shading interp, colormap(hot)
title('Spectrogram of Signal')
xlim([0 t(end)]);
xlabel('Time')
ylabel('Frequency')
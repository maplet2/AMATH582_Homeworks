clear all; close all; clc
% This code will take data that contains a an ultrasound of a dog's stomach
% that has swallowed a marble and track/predict the trajectory of the
% marble through the dog.

% Load Initial Data
load('Testdata.mat')

L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); 
x=x2(1:n); 
y=x; 
z=x;

% Create Frequency Domain
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k);

% Create Meshgrid for both Domains 
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Plot one snapshot of the data
figure(1)
Un(:,:,:)=reshape(Undata(1,:),n,n,n);
close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Raw Ultrasound Data at First Timeshot')
xlabel('x')
ylabel('y')
zlabel('z')

% Fourier Transform the Data and Take the Sum over the 20 Timeshots
Untavg = zeros(n,n,n);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Unf = fftn(Un);
    Untavg = Untavg + Unf;
end

% Take the Average over the 20 Timeshots, Normalize the Data with
%   the Maximum and then Graph the Averaged Signal in Frequency Domain
Utavg = Untavg / 20;
Utavgshift = fftshift(Untavg) / 20;
[maxUtavg, max_ind] = max(Utavgshift(:));

figure(2)
isosurface(Kx,Ky,Kz,abs(Utavgshift)./abs(maxUtavg),0.4)
axis([-10 10 -10 10 -10 10]), grid on, drawnow 
title('Averaged Frequency over Different Snapshots')
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')


% Pinpoint Location of Maximum signal Strength which is our Marble's Signal
%	Note that ind2sub returns in [y,x,z] form 
[Ymax, Xmax, Zmax] = ind2sub([n,n,n], max_ind);
kxmax = ks(Xmax);
kymax = ks(Ymax);
kzmax = ks(Zmax);


%%
% Use Marble Location to Create Filter
filter_eq = @(x,y,z)(exp(-0.2*((x-kxmax).^2 + (y-kymax).^2 + (z-kzmax).^2)));
filter = filter_eq(Kx,Ky,Kz);

%%
% Apply Filter to Fourier Transformed Data in all 20 different timeshots
%   before inverse transforming the data back to the time domain. Within
%   the time domain, we will locate the marble's location and track its
%   movement through the different timeshots. 
X_trav = [];
Y_trav = [];
Z_trav = [];
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Unf = fftn(Un);
    Unfs = fftshift(Unf);
    Unfs_fil = Unfs.*filter;
    Unfs_fil_us = ifftshift(Unfs_fil);
    Un_filtered = ifftn(Unfs_fil_us);
    [time_max, ind] = max(Un_filtered(:));
    [Ytmax,Xtmax,Ztmax] = ind2sub([n,n,n], ind);
    X_trav(j) = x(Xtmax);
    Y_trav(j) = y(Ytmax);
    Z_trav(j) = z(Ztmax);
end

% Plot Trajectory of Marble
clf;
figure(3)
hold on 
view(3)
plot3(X_trav, Y_trav, Z_trav)
plot3(X_trav(end), Y_trav(end), Z_trav(end), 'r.','MarkerSize',20)
title('Trajectory of Marble Through Dog''s Inside')
xlabel('x')
ylabel('y')
zlabel('z')
axis([-15 15 -15 15 -15 15])
grid on

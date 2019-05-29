clear; clc; close all; 

% Load Data
test1 = load_dat('camera1_1.mat','camera2_1.mat','camera3_1.mat');
test2 = load_dat('camera1_2.mat','camera2_2.mat','camera3_2.mat');
test3 = load_dat('camera1_3.mat','camera2_3.mat','camera3_3.mat');
test4 = load_dat('camera1_4.mat','camera2_4.mat','camera3_4.mat');


%%
X = test3;
[~,n] = size(X);

%%
mn = mean(X, 2);
X=X-repmat(mn,1,n);
[u,s,v]=svd(X'/sqrt(n-1),0); % perform the SVD
lambda=diag(s).^2; % produce diagonal variances

figure(1)
clf;
% Plot Variance in Principal Components
subplot(1,2,1)
plot(lambda/sum(lambda), 'ro', 'linewidth', 2)
title('Principal Components of Data')
xlabel('PCA Component number')
ylabel('Variance in Percentage')
ylim([0 1])

% Plot Displacement in New Coordinate System
subplot(1,2,2)
project = v'*X;
hold on
plot(1:length(X), project(1:2,:), 'linewidth', 1.5)
plot(1:length(X), project(3,:), 'linewidth', 0.5)
title('Graph of Displacement in New Coordinate System')
xlabel('Time')
ylabel('Position')
xlim([0, length(X)])
legend('PC 1','PC 2', 'PC 3', 'location', 'southeast')
print(gcf, '-dpng', 'test3_PCA.png')


% % figure(2)   
% % clf;
% % hold on
% % plot(X(1,:), X(2,:), 'r.', 'markersize', 15)
% % plot(X(3,:), X(4,:), 'g.', 'markersize', 15)
% % plot(X(5,:), X(6,:), 'b.', 'markersize', 15)
% % xlim([0 640])
% % ylim([0 480])
% 
% %%
% X = test2;
% [m,n] = size(X);
% 
% figure(2)
% clf;
% hold on
% plot(X(1,:), X(2,:), 'r.', 'markersize', 15)
% plot(X(3,:), X(4,:), 'g.', 'markersize', 15)
% plot(X(5,:), X(6,:), 'b.', 'markersize', 15)
% xlim([0 640])
% ylim([0 480])
% 
% mn = mean(X, 2);
% X=X-repmat(mn,1,n);
% [u,s,v]=svd(X'/sqrt(n-1),'econ'); % perform the SVD
% lambda=diag(s).^2; % produce diagonal variances
% figure(1)
% clf;
% project = u'*X; 
% plot(lambda/sum(lambda), 'ro', 'linewidth', 2)
% ylim([0 1])
% 
% %%
% X = test3;
% [m,n] = size(X);
% 
% figure(2)
% clf;
% hold on
% plot(X(1,:), X(2,:), 'r.', 'markersize', 15)
% plot(X(3,:), X(4,:), 'g.', 'markersize', 15)
% plot(X(5,:), X(6,:), 'b.', 'markersize', 15)
% xlim([0 640])
% ylim([0 480])
% 
% mn = mean(X, 2);
% X=X-repmat(mn,1,n);
% [u,s,v]=svd(X/sqrt(n-1),0); % perform the SVD
% lambda=diag(s).^2; % produce diagonal variances
% figure(1)
% clf;
% project = u'*X; 
% plot(lambda/sum(lambda), 'ro', 'linewidth', 2)
% ylim([0 1])
% %%
% X = test4;
% [m,n] = size(X);
% 
% figure(2)
% clf;
% hold on
% plot(X(1,:), X(2,:), 'r.', 'markersize', 15)
% plot(X(3,:), X(4,:), 'g.', 'markersize', 15)
% plot(X(5,:), X(6,:), 'b.', 'markersize', 15)
% xlim([0 640])
% ylim([0 480])
% 
% mn = mean(X, 2);
% X=X-repmat(mn,1,n);
% [u,s,v]=svd(X/sqrt(n-1),0); % perform the SVD
% lambda=diag(s).^2; % produce diagonal variances
% figure(1)
% clf;
% plot(lambda/sum(lambda), 'ro', 'linewidth', 2)
% ylim([0 1])
% 
% figure(3)
% clf;
% project = u'*X;
% plot(1:length(X), project(1:3,:))

%%

% This Function will take in the three different file names for each
%   of the tests and will create a 6 row matrix that contains the data
%   for each camera. It will also align the phases of the data. 
function test = load_dat(mat1, mat2, mat3)
    cam1 = cell2mat(struct2cell(load(mat1, 'min*')));
    cam2 = cell2mat(struct2cell(load(mat2, 'min*')));
    cam3 = cell2mat(struct2cell(load(mat3, 'min*')));
    
    [~, ind1] = min(cam1(2,1:40));
    [~, ind2] = min(cam2(2,1:40));
    [~, ind3] = max(cam3(1,1:40));
    
    cam1_al = cam1(:,ind1:length(cam1));
    cam2_al = cam2(:, ind2:length(cam2));
    cam3_al = cam3(:, ind3:length(cam3));
    
    len = min([length(cam1_al),length(cam2_al),length(cam3_al)]);
    
    test = [cam1_al(:, 1:len); cam2_al(:,1:len); cam3_al(:,1:len)]; 
end



clear; close all; clc

%% Load datasets
load('D:\JinYang\MATLAB\3D_FE_Global_DVC\images_StantonGodshall\vol_1.mat');
Img{1} = vol{1};
load('D:\JinYang\MATLAB\3D_FE_Global_DVC\images_StantonGodshall\vol_4.mat');
Img{2} = vol{1};

%%%%% Visualization: before registration %%%%%
figure, imagesc3(double(Img{1})-double(Img{2})); 
colorbar; title('Before registration: f - g')


%%%%% Assign img registration option %%%%%
[optimizer, metric] = imregconfig('monomodal')

%%%%% Image size %%%%%
imageSize = [1379,517,35];
Rmoving = imref3d(imageSize);
Rfixed = imref3d(imageSize);


%%%%% Image registration %%%%%
tform = imregtform(Img{2},Rmoving,Img{1},Rfixed, 'translation', optimizer, metric);
%%%%% Solved displacement vector is:   tform.T(4,1:3);
tform.T

%%%%% Warp original image using the registrated translation %%%%%
movingRegistered = imwarp(Img{2},tform,'OutputView',imref3d(size(Img{1})));


%%%%% Visualization: after registration %%%%%
figure, imagesc3(double(Img{1})-double(movingRegistered)); 
colorbar; title('After registration: f - g')
 
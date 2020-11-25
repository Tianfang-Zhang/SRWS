%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is an implement of the Self-Regualrized Weighted Sparse (SRWS)
% model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code in your publications, please cite:
% Zhang T, Peng Z, Wu H, et al. Infrared small target detection via self-regularized weighted sparse model[J]. Neurocomputing, 420: 124-148.
% @article{zhang420infrared,
%   title={Infrared small target detection via self-regularized weighted sparse model},
%   author={Zhang, Tianfang and Peng, Zhenming and Wu, Hao and He, Yanmin and Li, Chaohai and Yang, Chunping},
%   journal={Neurocomputing},
%   volume={420},
%   pages={124--148},
%   publisher={Elsevier}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Tianfang Zhang
% Email: sparkcarleton@gmail.com
% Copyright:  University of Electronic Science and Technology of China
% Date: 2019/7/3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* License: Our code is only available for non-commercial research use.

clc;    clear;  close all;

Img = imread('images\1.bmp');

if ndims(Img) == 3
    Img = rgb2gray(Img);
end
Img = im2double(Img);
figure,subplot(121),imshow(Img),title('Original Image');

alg = SRWS;
alg = alg.process(Img);
subplot(122),imshow(alg.result, []),title('result');

figure,
subplot(131),imshow(alg.rstB, []),title('Background');
subplot(132),imshow(alg.rstE, []),title('Noise');
subplot(133),imshow(alg.Z, []),title('Coefficient Matrix');
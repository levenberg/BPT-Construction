function calculate()
clc
clear all;
I=imread('airplane.tif');
YCbCr=rgb2hsv(I);
Y =YCbCr(:,:,1);
Cb=YCbCr(:,:,2);
Cr=YCbCr(:,:,3);
[bh, bw,dim] = size(YCbCr);
norder=5;
amin=15;
amax=8000;
Y_data = double(reshape(Y', 1, bh*bw));
Cb_data = double(reshape(Cb', 1, bh*bw));
Cr_data = double(reshape(Cr', 1, bh*bw));
tic
scopstmp = scop_mex(Y_data,Cb_data,Cr_data, bh, bw, norder,amin, amax);
toc
scops = prepross1(scopstmp);
save('scop.mat','scops');
end
function scops = prepross1(scopstmp)
    nfeat = scopstmp.m_attri;
    scops.shapes = adfeature(single(scopstmp.shapes), nfeat);
    scops.pattern2 = transpose(uint32(reshape(scopstmp.pattern2, 2, length(scopstmp.pattern2)/2)));
    scops.pattern3 = transpose(uint32(reshape(scopstmp.pattern3, 3, length(scopstmp.pattern3)/3)));
    scops.pattern4 = transpose(uint32(reshape(scopstmp.pattern4, 3, length(scopstmp.pattern4)/3)));
end

function features = adfeature(features, nfeat)
if( isempty(features))
    features = zeros(1, nfeat);
end
features = reshape(single(features), nfeat, length(features)/nfeat);
features = features';
features(isnan(features)) = 0;
end


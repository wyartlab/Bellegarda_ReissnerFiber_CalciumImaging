
function [dff, dff2, raw] = calc_dff(Mr, Mask,n);

% This function takes the input Mask and rois and calculate deltaF/F by
% first fitting a polynomial to the raw data to remove photobleaching then
% allowing input by user of timepoints A and B from which the baseline will
% be selected. 

% n number of ROIs, 
%Mask cells for ROI defined with x and y based on ginput function see txt file with commands

% Mr is the registered M file
% Mask is the mask for all ROIs
% n is the number of ROIs selected

for i=1:n;
    
    ind=[];bgf=[];obgf=[];cbgf=[];fresult=[];gof=[];output=[];reply=[];
    
ind = repmat(find(Mask{i}),[1 size(Mr,3)]);
ind = ind + repmat((0:(size(Mr,3)-1))*(size(Mr,1)*size(Mr,2)),[nnz(Mask{i}) 1]);
ind = ind(:);

bgf = double(Mr(ind));
bgf = reshape(bgf,nnz(Mask{i}),size(Mr,3));
bgf = mean(bgf).';

obgf = bgf;

% t is time in frames
t = (0:(size(Mr,3)-1)).';
[fresult,gof,output] = fit(t,bgf,'exp2','Normalize','on');

% cbgf est le vecteur pour chaque frame du fit a 2 exp
cbgf = feval(fresult,t);

% normalise a 1
cbgf = cbgf/cbgf(1);
tr = bgf./cbgf;

figure;plot(tr);title(['CSFcN number ' num2str(i)],'FontSize', 18)
baseline = ginput;
time1(i) = round(baseline(1,1));
time2(i) = round(baseline(2,1));
clear baseline;
close all

raw (:,i)= bgf;

% baseline is taken from chosen timepoints

dff(:,i) = (tr/mean(tr(time1(i):time2(i)))-1)*100;
dff2(:,i) = tr;

end;





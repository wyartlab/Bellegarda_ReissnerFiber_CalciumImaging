function output = dffintegral_CB_ventral_MM5(framesel,fps, rois_v, trial, condition);
 
output=struc('M', [], 'framesel', [], 'fps', [], 'base_KAs',[], 'dff_KAs', [], 'dff_MM5_KAs', [],'raw_KAs', [], 'raw_MM5_KAs', [], 'int_KAs', [], 'ind_int', [], 'rois_KAs', [], 'val_int', []);

%clear all
close all
disp ('Choose image acquisition file')
[image folder] = uigetfile('*.*');
imagefile = strcat(folder, image);


% Select the image with high contrast to choose ROIs from

disp ('Choose image to select ROIs')
[file2, folder2]= uigetfile('*.*');
StanDev = strcat(folder2, file2);

% StanDev = strcat(folder, 'AVG_', image);
%nframes = input('How many frames is the image file? ');
%freq = input('What is the frequency of acquisition in Hz? ');
nframes = framesel(end);
freq = fps;

M=multitiff2M(imagefile,1:nframes);

output.M=M;
output.framesel=framesel;
output.fps=fps;

% if you need to register due to motion artefact
% Mr = registerfilm_ROI(M(:,:,1:nframes),1);

x_dim = length(M(1,:,1));
y_dim = length(M(:,1,1));

if isempty(trial)
    trial =1;
end

trial=num2str(trial);
condition=num2str(condition);

if isempty(rois_v)
    
figure;I=imread(StanDev);imshow(I);
disp('Choose ROIs for ventral CSF-contacting neurons from rostral to caudal')
rois_v = getROIcell;
close all

end;
    
save 

% Ventral: generate mask from rois, calculate dff, find and threshold
% minima, fit poly to minima, subtract dff trace from poly function,
% calculate integral normalized per roi by integral of the signal per second (used to be *60 in the code)


for i=1:size(rois_v,2); 
    Mask_v{i} = roipoly(imread(imagefile),rois_v{i}(:,1),rois_v{i}(:,2));
    i
end;

output.rois_KAs=rois_v;

% Save masks for rois_v variable:
    savefile = ['rois ' trial condition '.mat'];
    save(savefile, 'rois_v')

%save figure ROIs:
figure; imshow(StanDev);
for i=1:length(rois_v);
    patch(rois_v{1,i}(:,1),rois_v{1,i}(:,2),'r','FaceAlpha',0.55);
    text(rois_v{1,i}(1,1),rois_v{1,i}(1,2),num2str(i),'Color','g','FontSize', 8);
    hold on;
end;
title('ROIs','FontSize', 18);
rois = strcat(folder, 'ROIs');
saveas(gcf, rois,'fig'); 
saveas(gcf, rois,'png');

[dff_v, raw_v] = calc_dff(M, Mask_v, size(rois_v,2));

% Save dff and raw variables:
    savefile = ['dff and raw ' trial condition '.mat'];
    save(savefile, 'dff_v', 'raw_v')

% keep the raw and running average with factor 5

dff_MM5_v(:,i)= movmean(dff_v(:,i),5);
raw_MM5_v(:,i)= movmean(raw_v(:,i),5);

dff_KAs=dff_v;
raw_KAs=raw_v;

output.dff_KAs=dff_v;
output.dff_MM5_KAs=dff_MM5_v;

output.raw_KAs=raw_v;
output.raw_MM5_v=raw_MM5_v;


% save workspace here
    savefile = ['rois dff and raw ' trial condition '.mat'];
    save(savefile, 'rois_v','dff_v', 'raw_v', 'output')


for i=1:size(rois_v,2); 
    i
    min_v=[];
    ind_v=[];
    
    [min_v, ind_v] = lmin(dff_v(:,i),100)
    
    if length(ind_v)>2
        
        [fit_v, gof_v, out_v] = fit(ind_v', min_v', 'poly2', 'Normalize', 'on');
        base_v (:,i) = feval(fit_v, [1:nframes])';
%    adj_v(:,i) = dff_MM5_v(:,i)-base_v(:,i);   
        adj_v(:,i) = dff_v(:,i)-base_v(:,i);
        int_v(:,i) = trapz(adj_v(:,i));
    %int_v(:,i) = int_v(:,i)/nframes*freq*60;
        int_v(:,i) = int_v(:,i)/nframes*freq;
        
    else if length(ind_v)<=2
            
         adj_v(:,i) = dff_v(:,i)-mean(dff_v(:,i));
         int_v(:,i) = trapz(adj_v(:,i));
         int_v(:,i) = int_v(:,i)/nframes*freq;
    end;
    end;
end;


output.base_KAs=base_v;
output.int_KAs=int_v;

% sort the ROIs based on the integral of the calcium signals

val_int=[];
ind_int=[];
[val_int,ind_int]=sort(int_v,2,'descend')

output.val_int=val_int;
output.ind_int=ind_int;

% save each DFF per ROI and the corresponding baseline in red

figure; hold on;title(['DFF of CSF-contacting neurons from 1 to 10  and baseline in red ' trial],'FontSize', 18);

for i=1:min(10,size(rois_v,2)); 
    
    plot(dff_v(:,i)+100*i,'o-','LineWidth',2); hold on; plot (base_v(:,i)+100*i,'r');

end;

dff_baseline = strcat(folder, ['DFF of CSF-contacting neurons from 1 to 10  and baseline in red ' trial]);
saveas(gcf, dff_baseline,'fig'); 
saveas(gcf, dff_baseline,'png');

% save each DFF per ROI and the corresponding baseline in red

figure; hold on;title(['DFF of CSF-contacting neurons from last 10 and baseline in red ' trial],'FontSize', 18);

for i=1:min(10,size(rois_v,2)); 
    
    plot(dff_v(:,end+1-i)+100*i,'o-','LineWidth',2); hold on; plot (base_v(:,end+1-i)+100*i,'r');

end;

dff_baseline = strcat(folder, ['DFF of CSF-contacting neurons from last 10 and baseline in red ' trial]);
saveas(gcf, dff_baseline,'fig'); 
saveas(gcf, dff_baseline,'png');


% save values of normalized integrals

figure; title(['Normalized integrals per ROI selected from rostral to caudal ' trial]);hold on; 
for i=1:size(int_v); 
    plot(int_v,'ko','LineWidth',5); 
end; 
integral = strcat(folder, ['Normalized integrals per ROI selected from rostral to caudal ' trial]);
saveas (gcf, integral,'fig');
saveas (gcf, integral,'png');

figure; title(['10 CSF-contacting neurons dff sorted descending by the integral value: largest integral ' trial],'FontSize', 18);

for i=1:min(10,size(rois_v,2)); 
    subplot(5,2,i);
    plot(raw_v(:,ind_int(i)),'k'); hold on; 
    plot(dff_v(:,ind_int(i)),'b'); 
    plot(adj_v(:,ind_int(i)),'r'); 
    title(['ROI', num2str(ind_int(i))]);
    hold off; 
end;
vccd = strcat(folder, ['10 CSF-cNs DFF sorted descending by the integral value: largest integral ' trial]);
saveas(gcf,vccd,'fig');
saveas(gcf,vccd,'png');

figure; title(['CSF-contacting neurons dff sorted ascending by the integral value: lowest integral' trial],'FontSize', 18);

for i=1:min(10,size(rois_v,2)); 
    subplot(5,2,i);
    plot(raw_v(:,ind_int(end+1-i)),'k'); hold on; 
    plot(dff_v(:,ind_int(end+1-i)),'b'); 
    plot(adj_v(:,ind_int(end+1-i)),'r'); 
    title(['CSFcN ', num2str(ind_int(end+1-i))]);
    hold off; 
end;
vccd = strcat(folder, ['CSF-cNs DFF sorted descending by the integral value' trial]);
saveas(gcf,vccd,'fig');
saveas(gcf,vccd,'png');

save analysis_CSFcNs_integral
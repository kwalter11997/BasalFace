%% Blur Threshold Analysis

%faces
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Blur\FindFacesBlur2')
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Blur\FindFacesBlur2')
addpath('C:\Users\Kerri\Dropbox\Kerri_Walter\FInD functions')

fnames = d;
fnames = {fnames.name};
fnames = fnames([3:end]); %remove . and ..

for subj = 1:23
    load(fnames{subj})
    faceThresh(subj,:) = thresholdEst(1:5)
end

%Remove outliers
upper = mean(faceThresh(:)) + (std(faceThresh(:)))*3; lower = mean(faceThresh(:)) - (std(faceThresh(:)))*3; 
outliers = upper<faceThresh | faceThresh<lower;
faceThresh(outliers)=NaN; %none

%VA
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Blur\AIM_VA2')
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Blur\AIM_VA2')

fnames = d;
fnames = {fnames.name};
fnames = fnames([3:end]); %remove . and ..

for subj = 1:23
    load(fnames{subj})
    aimThresh(subj,:) = threshAcuity(1:5)
end

%Remove outliers
upper = mean(aimThresh(:)) + (std(aimThresh(:)))*3; lower = mean(aimThresh(:)) - (std(aimThresh(:)))*3; 
outliers = upper<aimThresh | aimThresh<lower;
aimThresh(outliers)=NaN; %none

% cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Blur')
% load('faceThresh')
% load('aimThresh')
%% Fit
% % Get coefficients of a line fit through the data.
% coefficients = polyfit([1:5], nanmean(aimThresh),4);
% % Create a new x axis with exactly 1000 points (or whatever you want).
% xFit = linspace(min([1:5]), max([1:5]), 1000);
% % Get the estimated yFit value for each of those 1000 new x locations.
% yFit = polyval(coefficients , xFit);
% % Plot everything.
% figure
% plot([1:5], nanmean(aimThresh), 'b.', 'MarkerSize', 15); % Plot training data.
% hold on; % Set hold on so the next plot does not blow away the one we just drew.
% plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
% grid on;

bFace = polyfit([1:5], nanmean(faceThresh), 2);
slopeFace = bFace(1)

bAcuity = polyfit([1:5], nanmean(aimThresh), 2);
slopeAcuity = bAcuity(1)

%normalize
normFace = (faceThresh-min(faceThresh(:))) / (max(faceThresh(:)) - min(faceThresh(:)))
normAim = (aimThresh-min(aimThresh(:))) / (max(aimThresh(:)) - min(aimThresh(:)))

figure
subplot(1,2,1)
boxplot(normFace,'PlotStyle','compact','Colors','b')
subplot(1,2,2)
boxplot(normAim,'PlotStyle','compact','Colors','r')

bNormFace = polyfit([1:5], nanmean(normFace), 1);
slopeNormFace = bNormFace(1)

bNormAcuity = polyfit([1:5], nanmean(normAim), 1);
slopeNormAcuity = bNormAcuity(1)

%exponential fit
f = fit([1:5]',nanmean(normFace)','exp2')
figure
plot(f,[1:5]',nanmean(normFace)')

f = fit([1:5]',nanmean(normAim)','exp2')
figure
plot(f,[1:5]',nanmean(normAim)')

%% Figs
% 
% pdiff = 380/571;
% x = [0 .1*pdiff .2*pdiff .4*pdiff .8*pdiff]

figure 
subplot(1,2,2)
boxplot(normFace,'PlotStyle','compact','Colors','k')
hold on
xticks([1:5])
xticklabels({'0','.1','.2','.4','.8'})
xlabel('Blur Width (deg)')
ylabel('Normalized Threshold Estimate')
plot([1:5],nanmean(normFace), 'k')
title('Face Thresholds')

subplot(1,2,1)
boxplot(normAim,'PlotStyle','compact','Colors','k')
hold on
xticks([1:5])
xticklabels({'0','.1','.2','.4','.8'})
xlabel('Blur Width (deg)')
ylabel('Normalized Threshold Estimate')
plot([1:5],nanmean(normAim), 'k')
title('Acuity Thresholds')

% figure
% subplot(1,2,2)
% boxplot(normFace,'PlotStyle','compact','Colors','b')
% hold on
% xticks([1:5])
% xticklabels({'0','.0665','.1331','.2662','.5324'})
% xlabel('Blur Width (deg)')
% ylabel('Normalized Threshold Estimate')
% plot([1:5],nanmean(normFace), 'b')
% title('Face Thresholds')
% 
% subplot(1,2,1)
% boxplot(normAim,'PlotStyle','compact','Colors','r')
% hold on
% xticks([1:5])
% xticklabels({'0','.1','.2','.4','.8'})
% xlabel('Blur Width (deg)')
% ylabel('Normalized Threshold Estimate')
% plot([1:5],nanmean(normAim), 'r')
% title('Acuity Thresholds')

%fit lines
%linear
figure
scatter([1:5], nanmean(normFace),'ob','filled')
hold on
scatter([1:5], nanmean(normAim), 'or','filled')
h = lsline;
h(1).Color = 'r';
h(2).Color = 'b';
xticks([1:5])
xticklabels({'0','.1','.2','.4','.8'})
xlabel('Blur Width (deg)')
ylabel('Normalized Threshold Estimate')
title('Linear Least-Squares')

%exponential
figure
f1 = fit([1:5]',nanmean(normFace)','exp2') %f1(x) = 0.12*exp(-0.01*x) + 1.57e-06*exp(2.35*x)
f2 = fit([1:5]',nanmean(normAim)','exp2') %f2(x) = 1348*exp(-0.16*x) + -1348*exp(-0.16*x)

scatter([1:5], nanmean(normFace),'ob','filled')
hold on
scatter([1:5], nanmean(normAim), 'or','filled')
plot(f1,'b')
plot(f2,'r')

legend('Face', 'Acuity')

%power
figure
f1 = fit([1:5]',nanmean(normFace)','power1') 
f2 = fit([1:5]',nanmean(normAim)','power1') 
scatter([1:5], nanmean(normFace),'ob','filled')
hold on
scatter([1:5], nanmean(normAim), 'or','filled')
plot(f1,'b')
plot(f2,'r')
%% Stats
%find individual slopes
figure
hold on
for t = 1:23
    data = normFace(t,:)'; 
%     x = [0 .1*pdiff .2*pdiff .4*pdiff .8*pdiff]';
%     if isnan(data(5))
%         data = normFace(t,1:4)';
% %         x = [0 .1*pdiff .2*pdiff .4*pdiff]';
%     end
    %ff = fit(x,data','exp2');
%     ff=polyfit(x, data, 1)
    ff=polyfit([0,.1,.2,.4,.8]', data, 1)
    %faceSlope(t) = ff.b
    faceSlope(t) = ff(1)
%     scatter(x,data,'ob')
    scatter([0,.1,.2,.4,.8]',data,'ob')
    
    %fa = fit([1:5]',normAim(t,:)','exp2');
    fa=polyfit([0,.1,.2,.4,.8]', normAim(t,:)', 1)
    %aimSlope(t) = fa.b
    aimSlope(t) = fa(1)
    scatter([0,.1,.2,.4,.8], normAim(t,:),'or')
end
ylim([0,1])
lsline

[h,p,ci,stats] = ttest(faceSlope',aimSlope') %t(22)=-11.291; p<.001)
computeCohen_d(aimSlope',faceSlope') %d=4.152

%% Threshold analysis - FInD Faces and FInD Faces Combos (N=8)
%remote study

d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Singles');
fnamesSing = d;
fnamesSing = {fnamesSing.name};
fnamesSing = fnamesSing(3:end); %remove . and ..
    
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Combos');
fnamesComb = d;
fnamesComb = {fnamesComb.name};
fnamesComb = fnamesComb(3:end); %remove . and ..

%build single matrix
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Singles')
for subj = 1:8
    load(fnamesSing{subj})
    singThresh(subj,:) = thresholdEst
end
%build combo matrix
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Combos')
for subj = 1:8
    load(fnamesComb{subj})
    combThresh(subj,:) = thresholdEst
end

mean(singThresh')
mean(combThresh')

%ratio of difference
ratio = mean(singThresh(:))/mean(combThresh(:)) %total
r = mean(singThresh,2) ./ mean(combThresh,2) %by subj

[h,p,ci,stats] = ttest(r,1.4) %test against probability summation supperadditivty (1.4)
d = computeCohen_d(1.4, r, 'independent')

%means of individual parameters when combined
for p=1:10
    indAvg = mean(singThresh);
    if p==10
       indComb(p) = (indAvg(p)+indAvg(1))/2
    else
       indComb(p) = (indAvg(p)+indAvg(p+1))/2
    end
end
indComb

%box
figure
boxplot(singThresh,'PlotStyle','compact','Colors','b')
hold on
boxplot(combThresh,'positions', [1.5:10.5],'PlotStyle','compact','Colors','r')
xticks(1:11)
xticklabels({'1','2','3','4','5','6','7','8','9','10','1'})
xlim([0.5,11])
ylim([0,13])
xlabel('Component')
ylabel('Threshold Estimate')
%mean lines
% plot([1:10],mean(singThresh), 'b')
% hold on
% plot([1.5:10.5],mean(combThresh), 'r')
plot([1.5:10.5],mean(combThresh), 'r')
plot([1.5:10.5],indComb, 'b')
%add manual legend 
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s','MarkerFaceColor','b');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','r');
legend(h, 'Single','Combined','Location','northwest');

%across subjs
figure
for s=1:8
    [~,p] = ttest(singThresh(s,:)',combThresh(s,:)')
    subplot(2,4,s)
    boxplot([singThresh(s,:)',combThresh(s,:)'],'Labels',{'Single','Combined'})
    if p<.05
        sigstar([1,2],p)
    end
    title(sprintf('Subj%d', s))
    ylabel('Threshold Estimate')
end

[h,p,ci,stats] = ttest(mean(singThresh'),mean(combThresh')) %paired ttest
d = computeCohen_d(mean(singThresh'), mean(combThresh'), 'paired')

singcomb = [singThresh;combThresh]
std(singcomb(:))

%% Look at all data for singles (includes pilot data. N=16)

d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\pilotData_1-50');
fnamesPilot = d;
fnamesPilot = {fnamesPilot.name};
fnamesPilot = fnamesPilot(3:end); %remove . and ..

%fnamesAll = [fnamesPilot,fnamesSing];

%build matrix
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\pilotData_1-50');
for subj = 1:8
    load(fnamesPilot{subj})
    pilotThresh(subj,:) = thresholdEst(1:10) %comparing only thresholds 1-10
end

allThresh = [pilotThresh;singThresh]

%box
figure
boxplot(allThresh,'PlotStyle','compact','Colors','k')
hold on
xlabel('Parameter Number')
ylabel('Threshold Estimate')
%mean lines
plot([1:10],mean(allThresh), 'k')

[p,tbl,stats] = anova1(allThresh)
results = multcompare(stats)

%% compare in-lab and remote studies

boxplot(pilotThresh); hold on
boxplot(singThresh)

[h,p,ci,stats] = ttest2(mean(pilotThresh,2),mean(singThresh,2)) %2sample ttest
d = computeCohen_d(mean(pilotThresh,2), mean(singThresh,2), 'independent')

%% find average chart time
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Singles')
for subj = 1:8
    load(fnamesSing{subj})
    singTime(subj,:) = expDuration
end
singChartTimes = [singTime/40]'
singThreshTimes = [singTime/10]'
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\Combos')
for subj = 1:8
    load(fnamesComb{subj})
    combTime(subj,:) = expDuration
end
singCombTimes = [singTime;combTime]
expTime = mean(singCombTimes)
chartTimesec = expTime/40
threshTimesec = expTime/10
threshTimemin = threshTimesec/60

%charts
chartTimesec = mean(singCombTimes/40)
chartTimestd = std(singCombTimes/40)

%% Single and Combo Yaw

d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\SingleYaw');
fnamesYawSing = d;
fnamesYawSing = {fnamesYawSing.name};
fnamesYawSing = fnamesYawSing(3:end); %remove . and ..
    
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\ComboYaw');
fnamesYawComb = d;
fnamesYawComb = {fnamesYawComb.name};
fnamesYawComb = fnamesYawComb(3:end); %remove . and ..

%build single matrix
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\SingleYaw')
for subj = 1:8
    load(fnamesYawSing{subj})
    yawSingThresh(subj,:) = thresholdEst
end
%build combo matrix
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\ComboYaw')
for subj = 1:8
    load(fnamesYawComb{subj})
    yawCombThresh(subj,:) = thresholdEst
end

mean(yawSingThresh')
mean(yawCombThresh')

%Ratio
ratio = mean(yawSingThresh(:)) / mean(yawCombThresh(:))
r = mean(yawSingThresh,2) ./ mean(yawCombThresh,2)

[h,p,ci,stats] = ttest(r,1.4)
d = computeCohen_d(1.4,r, 'independent')

%means of individual parameters when combined
for p=1:10
    indYawAvg = mean(yawSingThresh);
    if p==10
       indYawComb(p) = (indYawAvg(p)+indYawAvg(1))/2
    else
       indYawComb(p) = (indYawAvg(p)+indYawAvg(p+1))/2
    end
end
indYawComb
mean(yawCombThresh)

singcombyaw = [yawSingThresh;yawCombThresh]
std(singcombyaw(:))

%box
figure
boxplot(yawSingThresh,'PlotStyle','compact','Colors','b')
hold on
boxplot(yawCombThresh,'positions', [1.5:10.5],'PlotStyle','compact','Colors','r')
xticks(1:11)
xticklabels({'1','2','3','4','5','6','7','8','9','10','1'})
xlim([0.5,11])
ylim([0,30])
xlabel('Component')
ylabel('Threshold Estimate')
%mean lines
% plot([1:10],mean(singThresh), 'b')
% hold on
% plot([1.5:10.5],mean(combThresh), 'r')
plot([1.5:10.5],mean(yawCombThresh), 'r')
plot([1.5:10.5],indYawComb, 'b')
%add manual legend 
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s','MarkerFaceColor','b');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','r');
legend(h, 'Single','Combined','Location','northwest');

[h,p,ci,stats] = ttest(mean(yawSingThresh'),mean(yawCombThresh')) %paired ttest sing vs comb
d = computeCohen_d(mean(yawSingThresh,2), mean(yawCombThresh,2), 'paired')

indYawComb-mean(yawCombThresh) %difference per parameter

%across subjs
figure
for s=1:8
    [~,p] = ttest(yawSingThresh(s,:)',yawCombThresh(s,:)')
    subplot(2,4,s)
    boxplot([yawSingThresh(s,:)',yawCombThresh(s,:)'],'Labels',{'Single','Combined'})
    if p<.05
        sigstar([1,2],p)
    end
    title(sprintf('Subj%d', s))
    ylabel('Threshold Estimate')
end

[h,p,ci,stats] = ttest2([mean(yawSingThresh'),mean(yawCombThresh')],[mean(singThresh'),mean(combThresh')]) %two sample ttest foward vs yaw
d = computeCohen_d([mean(yawSingThresh'),mean(yawCombThresh')], [mean(singThresh'),mean(combThresh')], 'independent')

%chart time
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\SingleYaw')
for subj = 1:8
    load(fnamesYawSing{subj})
    yawSingTime(subj,:) = expDuration
end
cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\ComboYaw')
for subj = 1:8
    load(fnamesYawComb{subj})
    yawCombTime(subj,:) = expDuration
end
yawSingCombTimes = [yawSingTime;yawCombTime]
expTime = mean(yawSingCombTimes)
chartTime = expTime/40
threshTime = expTime/10

%charts
chartTimesec = mean(yawSingCombTimes/40)
chartTimestd = std(yawSingCombTimes/40)

[h,p,ci,stats] = ttest2(singCombTimes,yawSingCombTimes) %2sample ttest for times
d = computeCohen_d(singCombTimes, yawSingCombTimes, 'independent')

% [p,tbl,stats] = anova1(yawThresh)
% results = multcompare(stats)

%% sing/comb + yaw sing/comb ANOVA

% anovaMat = [singThresh;combThresh;yawSingThresh;yawCombThresh]
% [p,tbl,stats] = anova2(anovaMat,8) %two way anova
% multcompare(stats)

%% Confidence intervals
x = yawCombThresh(:);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals

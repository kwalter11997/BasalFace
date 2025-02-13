%% CVI Threshold Analysis
%foward facing and roatated

cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_Data')
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_Data')
addpath('C:\Users\Kerri\Dropbox\Kerri_Walter\FInD functions')

fnamesCVI = d;
fnamesCVI = {fnamesCVI.name};
fnamesCVI = fnamesCVI([5:end-1]); %remove . and .. and folders
% 
% %cvi thresholds from raw datafiles were collected with old version of FitDPrime, recalculate
for subj = 1:8
    load(fnamesCVI{subj})
    nPCs = 10; %10 parameters being tested
    thresholdEst=zeros(1, nPCs);
    for pcNum=1:nPCs % process data foreach Spatial Frequency
        testLevel=[]; % set blank list of all levels for this SF
        sYes=[]; % blank list of all stimuli seen for this SF
        for trialNo=1:nTrials % work through all trials
            testLevel=[testLevel trialRecord(pcNum,trialNo).targLevelPerLoc]; % concatenate all levels
            sYes=[sYes trialRecord(pcNum,trialNo).stimYes]; % concatenate all Yes responses
        end
        myData=arrangeData(testLevel, sYes);
        tLevel=myData(:,1);
        pYes=myData(:,4);
        errEst=myData(:,8);
        tLevel(1)=tLevel(2)/2;
        [fitObj, ~, thresholdEst(pcNum), err95ci]=fitDPrimeSensFunc( tLevel, pYes, errEst, 1 );
    end
    cviThresh(subj,:) = thresholdEst(1:10)
end

load('cviThresh_updated-FitDPrime_8subj')

%Remove outliers
upper = mean(cviThresh(:)) + (std(cviThresh(:)))*3; lower = mean(cviThresh(:)) - (std(cviThresh(:)))*3; 
outliers = upper<cviThresh | cviThresh<lower;
cviThresh(outliers)=NaN;

nanmean(cviThresh,2)

d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_Data\Yaw')
fnamesCVI_yaw = d;
fnamesCVI_yaw = {fnamesCVI_yaw.name};
fnamesCVI_yaw = fnamesCVI_yaw(3:end); %remove . and .. and folder

% cvi thresholds from raw datafiles were collected with old version of FitDPrime, recalculate
for subj = 1:8
    load(fnamesCVI_yaw{subj})
    nPCs = 10; %10 parameters being tested
    thresholdEst=zeros(1, nPCs);
    
    for pcNum=1:nPCs % process data foreach Spatial Frequency
%         figure(pcNum); % create new figure
        testLevel=[]; % set blank list of all levels for this SF
        sYes=[]; % blank list of all stimuli seen for this SF

        for trialNo=1:nTrials % work through all trials            
            %chart times 
            if isempty(trialRecord(pcNum,trialNo).chartTime) %if subj didn't get to the next block of 10
               chartTimes(pcNum,trialNo) = NaN;
            else
               chartTimes(pcNum,trialNo) = trialRecord(pcNum,trialNo).chartTime;
            end            
%             if trialRecord(pcNum,trialNo).chartTime < 2 %if subj moved through this trial in less than 2 seconds
%                trialRecord(pcNum,trialNo).chartTime = []; %remove trial
%                trialRecord(pcNum,trialNo).stimYes = [];
%                trialRecord(pcNum,trialNo).testLevel = [];
%                trialRecord(pcNum,trialNo).targLevelPerLoc = [];
%             end
            
            testLevel=[testLevel trialRecord(pcNum,trialNo).targLevelPerLoc]; % concatenate all levels
            sYes=[sYes trialRecord(pcNum,trialNo).stimYes]; % concatenate all Yes responses 
        end
    
        myData=arrangeData(testLevel, sYes);
        if isempty(myData)
            thresholdEst(pcNum) = NaN;
        else
            tLevel=myData(:,1);
            pYes=myData(:,4);
            errEst=myData(:,8);
            tLevel(1)=tLevel(2)/2;
            [fitObj, ~, thresholdEst(pcNum), err95ci]=fitDPrimeSensFunc( tLevel, pYes, errEst, 1 );
        end
    end
    
    cvi_yawThresh(subj,:) = thresholdEst(1:10)
    cviyawChartTimes(subj) = nanmean(chartTimes(:)); %avg chart time per subj
    cviyawNumTrials(subj) = unique(sum(~isnan(chartTimes),2)); %trial#
         
end

load('cvi_yawThresh_8subj')

%Remove outliers
upper = mean(cvi_yawThresh(:)) + (std(cvi_yawThresh(:)))*3; lower = mean(cvi_yawThresh(:)) - (std(cvi_yawThresh(:)))*3; 
outliers = upper<cvi_yawThresh | cvi_yawThresh<lower;
cvi_yawThresh(outliers)=NaN;

nanmean(cvi_yawThresh,2)

% Controls
d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_ControlData\Foward')
fnamesCon = d;
fnamesCon = {fnamesCon.name};
fnamesCon = fnamesCon(3:end); %remove . and .. and folder

% cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_ControlData\Foward')
% for subj = 1:8
%     load(fnamesCon{subj})
%     conThresh(subj,:) = thresholdEst(1:10)
% end

load('conThresh')

%Remove outliers
upper = mean(conThresh(:)) + (std(conThresh(:)))*3; lower = mean(conThresh(:)) - (std(conThresh(:)))*3; 
outliers = upper<conThresh | conThresh<lower;
conThresh(outliers)=NaN;

mean(conThresh,2)

d = dir('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_ControlData\Tilt')
fnamesCon_yaw = d;
fnamesCon_yaw = {fnamesCon_yaw.name};
fnamesCon_yaw = fnamesCon_yaw(3:end); %remove . and .. and folder
% 
% cd('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace\CVI_ControlData\Tilt')
% for subj = 1:8
%     load(fnamesCon_yaw{subj})
%     con_yawThresh(subj,:) = thresholdEst(1:10)
% end

load('con_yawThresh')

%Remove outliers
upper = mean(con_yawThresh(:)) + (std(con_yawThresh(:)))*3; lower = mean(con_yawThresh(:)) - (std(con_yawThresh(:)))*3; 
outliers = upper<con_yawThresh | con_yawThresh<lower;
con_yawThresh(outliers)=NaN;

mean(con_yawThresh,2)

%% Figs
figure
boxplot(conThresh,'PlotStyle','compact','Colors','b','Jitter',0)
%boxplot(conThresh,'Colors','b')
hold on
boxplot(cviThresh,'positions', [1.2:10.2],'PlotStyle','compact','Colors','r','Jitter',0)
%boxplot(cviThresh,'positions', [1.2:10.2],'Colors','r')
ylim([0,23])
xticks(1.1:11)
xticklabels({'1','2','3','4','5','6','7','8','9','10','1'})
xlim([0.5,11])
ylim([0,35])
xlabel('Parameter Number')
ylabel('Threshold Estimate')
%plot([1:10],mean(conThresh), 'b')
%plot([1.2:10.2],mean(cviThresh), 'r')
%mean lines
yline(mean(conThresh(:)),'--b')
yline(mean(cviThresh(:)),'--r')

%add manual legend 
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s','MarkerFaceColor','b');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','r');
legend(h, 'Control','CVI','Location','northwest');

%Yaw
figure
boxplot(con_yawThresh,'PlotStyle','compact','Colors','b','Jitter',0)
hold on
boxplot(cvi_yawThresh,'positions', [1.2:10.2],'PlotStyle','compact','Colors','r','Jitter',0)
ylim([0,23])
xticks(1.1:11)
xticklabels({'1','2','3','4','5','6','7','8','9','10','1'})
xlim([0.5,11])
ylim([0,35])
xlabel('Parameter Number')
ylabel('Threshold Estimate')
%plot([1:10],mean(con_yawThresh), 'b')
%plot([1.2:10.2],nanmean(cvi_yawThresh), 'r')

yline(mean(con_yawThresh(:)),'--b')
yline(nanmean(cvi_yawThresh(:)),'--r')

%add manual legend 
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'s','MarkerFaceColor','b');
h(2) = plot(NaN,NaN,'s','MarkerFaceColor','r');
legend(h, 'Control','CVI','Location','northwest');

%% Stats

%control vs cvi foward 2sample ttest
[h,p,ci,stats] = ttest2(mean(conThresh'),mean(cviThresh')) %cvi is sig higher than control (t(13) = -3.4387, p = .004)
computeCohen_d(mean(cviThresh'),mean(conThresh')) %1.7194

%control vs cvi yaw (temporary permutation)
% [p, observeddifference, effectsize] = permutationTest(mean(con_yawThresh,2), mean(cvi_yawThresh,2), 0, 'exact', 1, 'plotresult', 1) %cvi is sig higher than control (p=.018)

%control vs cvi yaw 2sample ttest
[h,p,ci,stats] = ttest2(mean(con_yawThresh'),nanmean(cvi_yawThresh')) %sig dif (t(13) = -2.1629, p = 0.0483)
computeCohen_d(mean(cvi_yawThresh'),mean(con_yawThresh')) %0.964

%control foward anova
[p,tbl,stats] = anova1(conThresh)
results = multcompare(stats) %4 lower than 7 8 and 10

%cvi foward anova
[p,tbl,stats] = anova1(cviThresh)
results = multcompare(stats) %no differences

%control yaw anova
[p,tbl,stats] = anova1(con_yawThresh)
results = multcompare(stats) %bunch of stuff

%cvi yaw anova
[p,tbl,stats] = anova1(cvi_yawThresh)
results = multcompare(stats) %no differences

%control foward vs yaw
[h,p,ci,stats] = ttest(mean(conThresh'),mean(con_yawThresh')) %paired ttest
computeCohen_d(mean(con_yawThresh'),mean(conThresh')) %1.574
%control yaw is sig higher than foward (t(7) = -3.8527; p=.006)

%cvi foward vs yaw
[h,p,ci,stats] = ttest(mean(cviThresh'),nanmean(cvi_yawThresh')) %paired ttest
computeCohen_d(mean(cvi_yawThresh'),mean(cviThresh'))
% %cvi yaw is not sig higher than foward (t(7) = -1.3551; p=0.2175)

%median
[h,p,ci,stats] = ttest(nanmedian(cviThresh'),nanmedian(cvi_yawThresh')) %paired ttest
computeCohen_d(nanmedian(cvi_yawThresh'),nanmedian(cviThresh'))
% (t(7) = -0.0851; p=0.9346)

% [p, observeddifference, effectsize] = permutationTest(mean(cviThresh,2), mean(cvi_yawThresh,2), 0, 'exact', 1, 'plotresult', 1) 
% %cvi yaw is not sig higher than foward (p=.777)



%control vs cvi ANOVA
data = [mean(conThresh,2);mean(cviThresh,2);mean(con_yawThresh,2);nanmean(cvi_yawThresh,2)]
concvi = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]'
fowyaw = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]'

[p,tbl,stats] = anovan(data,{concvi,fowyaw})
multcompare(stats)

%% Bar chart comparison

y = [mean(conThresh(:)) mean(cviThresh(:));mean(con_yawThresh(:)) nanmean(cvi_yawThresh(:))]
% x = categorical({'Foward', 'Tilt'})
b = bar(y)

% b(1).FaceColor = [1 1 1]; b(2).FaceColor = [.5 .5 .5];
hold on

% errlow = [mean(conThresh(:)) - std(conThresh(:)) mean(cviThresh(:)) - std(cviThresh(:)); mean(con_yawThresh(:)) - std(con_yawThresh(:)) nanmean(cvi_yawThresh(:)) - nanstd(cvi_yawThresh(:))]
% errhigh = [mean(conThresh(:)) + std(conThresh(:)) mean(cviThresh(:)) + std(cviThresh(:)); mean(con_yawThresh(:)) + std(con_yawThresh(:)) nanmean(cvi_yawThresh(:)) + nanstd(cvi_yawThresh(:))]

error = [std(conThresh(:)) std(cviThresh(:)); std(con_yawThresh(:)) nanstd(cvi_yawThresh(:))]

[ngroups,nbars] = size(y);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    t = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(t, y(:,i), error(:,i), 'k', 'linestyle', 'none');
end

set(gca,'xticklabel',{'Foward', 'Tilt'})
ylabel('Average Threshold')
ylim([0,17])

sigstar({[.85,1.15] [1.85 2.15]},[0.004,0.048])

legend(b, 'Control','CVI','Location','northwest');

%% chart time

%control foward
for subj=1:8
    load(fnamesCon{subj})
    for pc = 1:10
        for trial = 1:4
            if isempty(trialRecord(pc,trial).chartTime) %if subj didn't get to the next block of 10
               chartTimes(pc,trial) = NaN;
            else
               chartTimes(pc,trial) = trialRecord(pc,trial).chartTime;
            end
        end
    end
    conChartTimes(subj) = nanmean(chartTimes(:)); %avg chart time per subj
    conNumTrials(subj) = unique(sum(~isnan(chartTimes),2)); %trial#
    conExpDur(subj) = expDuration; %total exp duration
end

mean(conChartTimes)
std(conChartTimes)

%cvi foward
for subj=1:8
    load(fnamesCVI{subj})
    for pc = 1:10
        for trial = 1:4
            if isempty(trialRecord(pc,trial).chartTime) %if subj didn't get to the next block of 10
               chartTimes(pc,trial) = NaN;
            else
               chartTimes(pc,trial) = trialRecord(pc,trial).chartTime;
            end
        end
    end
    cviChartTimes(subj) = nanmean(chartTimes(:)); %avg chart time per subj
    cviNumTrials(subj) = unique(sum(~isnan(chartTimes),2)); %trial#
    cviExpDur(subj) = expDuration; %total exp duration
end

mean(cviChartTimes)
std(cviChartTimes)

%control tilt
for subj=1:8
    load(fnamesCon_yaw{subj})
    for pc = 1:10
        for trial = 1:4
            if isempty(trialRecord(pc,trial).chartTime) %if subj didn't get to the next block of 10
               chartTimes(pc,trial) = NaN;
            else
               chartTimes(pc,trial) = trialRecord(pc,trial).chartTime;
            end
        end
    end
    conyawChartTimes(subj) = nanmean(chartTimes(:)); %avg chart time per subj
    conyawNumTrials(subj) = unique(sum(~isnan(chartTimes),2)); %trial#
    conyawExpDur(subj) = expDuration; %total exp duration
end

mean(conyawChartTimes)
std(conyawChartTimes)

%cvi tilt
for subj=1:8
    load(fnamesCVI_yaw{subj})
    for pc = 1:10
        for trial = 1:4
            if isempty(trialRecord(pc,trial).chartTime) %if subj didn't get to the next block of 10
               chartTimes(pc,trial) = NaN;
            else
               chartTimes(pc,trial) = trialRecord(pc,trial).chartTime;
            end
        end
    end
    cviyawChartTimes(subj) = nanmean(chartTimes(:)); %avg chart time per subj
    cviyawNumTrials(subj) = unique(sum(~isnan(chartTimes),2)); %trial#
    cviyawExpDur(subj) = expDuration; %total exp duration
end

mean(cviyawChartTimes)
std(cviyawChartTimes)

timedata = [conChartTimes,cviChartTimes,conyawChartTimes,cviyawChartTimes]
g1 = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]; %control vs cvi
g2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %foward vs tilt
anovan(timedata,{g1,g2})

%% plot average random faces

con = mean(conThresh(:))
cvi = mean(cviThresh(:))
cony = mean(mean(con_yawThresh,2))
cviy = nanmean(cvi_yawThresh(:))

std(conThresh(:))
std(cviThresh(:))
std(con_yawThresh(:))
nanstd(cvi_yawThresh(:))

rng('default');
bfm_dir = 'C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace';
cd(bfm_dir); % switch to Basle Face Model Directory
[model msz] = load_model(); % load BFM
rp = defrp; % rp = render parameters for BFM
rp.width=400; % rp = render parameters for BFM
rp.height=400;

%Parameter 5 changing
for conditions = 1:4   
    shape_coords = normrnd(0, 1, [199,1]); % create a random vector for the 199-D structural PCA coordinates
    tex_coords = normrnd(0, 1, [199,1]); % ...and the same for the 199-D texture PCA coordinates
    
    shape_coords1 = shape_coords;
    shape_coords2 = shape_coords;
    
    % render and draw first face
    if conditions == 1
        shape_coords1(5) = shape_coords(5)-con/2 %control foward
        rp.phi = 0; % render front-facing view
    elseif conditions == 2
        shape_coords1(5) = shape_coords(5)-cvi/2 %cvi foward
        rp.phi = 0; % render front-facing view
    elseif conditions == 3
        shape_coords1(5) = shape_coords(5)-cony/2 %control tilt
        rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between -5 and +5 degrees
    elseif conditions == 4
        shape_coords1(5) = shape_coords(5)-cviy/2 %cvi tilt
        rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between -5 and +5 degrees
    end

    mesh1  = coef2object(shape_coords1, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
    tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
    h = figure(conditions*2-1);
    display_face(mesh1, tex1, model.tl, rp);
    set(gcf, 'Color', [ 0.5 0.5 0.5 ]);

    % render and draw second face
    if conditions == 1
        shape_coords2(5) = shape_coords(5)+con/2 %control foward
        rp.phi = 0; % render front-facing view
    elseif conditions == 2
        shape_coords2(5) = shape_coords(5)+cvi/2 %cvi foward
        rp.phi = 0; % render front-facing view
    elseif conditions == 3
        shape_coords2(5) = shape_coords(5)+cony/2 %control tilt
        rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between -5 and +5 degrees
    elseif conditions == 4
        shape_coords2(5) = shape_coords(5)+cviy/2 %cvi tilt
        rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between -5 and +5 degrees
    end
    mesh2  = coef2object(shape_coords2, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
    h = figure(conditions*2);
    display_face(mesh2, tex1, model.tl, rp);
    set(gcf, 'Color', [ 0.5 0.5 0.5 ]);
end

%% CVI demogaphics

% cFKE0796 = 27 20/20 F
% cFNW1194 = 27 20/20 F
% cFSD1001 = 22 20/20 F
% cFSM1198 = 25 20/15 F
% cFTC0798 = 25 20/70 F
% cMCF0999 = 24 20/40 M
% cMJB0304 = 19 20/20 M
% cMKS0804 = 19 20/25 M

x = [20,20,20,15,70,40,20,25]';
y = mean(cviThresh,2);
c = [0.8500 0.3250 0.0980];
scatter(x,y,'filled','MarkerFaceColor',c)
hold on 

corr(x,y) %pearsons r

% Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot
plot(xFit, yFit, '-', 'LineWidth', 1,'Color',c); % Plot fitted line.

% hl = lsline
% B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
% Slope = B(2)
% Intercept = B(1)
% %'y = .007x + 6.661'
% mdl = fitlm(x,y)
%%'R^2 = 0.0016, p = .925'

x = [20,20,20,15,70,40,20,25]';
y = nanmean(cvi_yawThresh,2);
c = [0.4660 0.6740 0.1880];
scatter(x,y,'filled','MarkerFaceColor',c)
hold on 

corr(x,y) %pearsons r

% Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot
plot(xFit, yFit, '-', 'LineWidth', 1, 'Color', c); % Plot fitted line.

% hl = lsline
% B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
% Slope = B(2)
% Intercept = B(1)
% %'y = -0.022x + 8.516'
% mdl = fitlm(x,y)
%%'R^2 = 0.010, p = .813'

xlabel('Acuity')
ylabel('Average Threshold')
xlim([10,75])

legend('Foward','R^2 = 0.002, p = .925','Tilt','R^2 = 0.010, p = .813')
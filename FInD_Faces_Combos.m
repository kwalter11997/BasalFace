% FInD Faces Combos
% Program to measure discrimination thresholds for Faces, based on Basel face database
% A sequence of grids containing pairs of faces is presented
% Target face pairs shares 198 random PCs, and differs along 2 PC of interest (PCs 1-10 are measured, stacked pairs (ie. 1+2, 2+3, 2+4))
% Distracter face pairs shares 199 random PCs
% Observer clicks on all cells contaiing different faces
% Staircase adjusts the PC difference on successive trials, based on d' sensitivity analysis
%
% 2022  PJB

clear all;

% addpath(genpath('C:\Users\bexlab\Dropbox\Kerri_Walter\BasalFace'))
% addpath(genpath('C:\Users\bexlab\Dropbox\Kerri_Walter\FInD functions'))
addpath(genpath('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace'))
addpath(genpath('C:\Users\Kerri\Dropbox\Kerri_Walter\FInD functions'))

% whichScreen=max(Screen('Screens')); % use highest display # for stimuli
whichScreen=2;
meanLUT=127;
% read in testing parameters:
prompt = {'subject initials','Targ Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)', 'Tilt? (0=no 1=yes)'};
dlg_title = 'FInD Depth';
num_lines = 1;
def = {'XX', '6', '3', '3', '4', '50', '50', '0'};
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1));
%targPC=str2num(char(answer(2,1)));
% targPC1=[1:10];
% targPC2=[2:10,1]; %shifted
targPC1=[1:10];
targPC2=[2:10,1]; %shifted
%threshEsts=str2num(char(answer(3,1)));
threshEsts=repmat(4,1,length(targPC1));
faceSize=str2num(char(answer(2,1)));
nRows=str2num(char(answer(3,1)));
nCols=str2num(char(answer(4,1)));
nTrials=str2num(char(answer(5,1)));
scrnWidthCm=str2num(char(answer(6,1))); % screen width (cm)
viewDistCm=str2num(char(answer(7,1))); % observer's distance from screen
tilt = str2num(char(answer(8,1)));

rng('default');
% bfm_dir = 'C:\Users\bexlab\Dropbox\Kerri_Walter\BasalFace';
bfm_dir = 'C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace';
cd(bfm_dir); % switch to Basle Face Model Directory
[model msz] = load_model(); % load BFM
rp = defrp; % rp = render parameters for BFM
rp.phi = 0; % render front-facing view

centerPCValMu=0; % TODO - this controls the mid point of the discrimination vector
centerPCValSD=0.5; % TODO - this controls the standard deviation of the mid point of the discrimination vector

saveImage=0;

try
    Screen('Preference','SkipSyncTests', 1);
    [windowPtr, winRect]=Screen('OpenWindow', whichScreen, meanLUT);%,[0 0 3200 1800]);
    
    frameRate=Screen('FrameRate', windowPtr);          % screen timing parameters
    scrnWidthDeg=2*atand((0.5*scrnWidthCm)/viewDistCm); % convert screen width to degrees visual angle
    scrnWidthPix=winRect(3)-winRect(1); % width of screen in pixels
    scrnHeightPix=winRect(4)-winRect(2); % width of screen in pixels
    pixPerDeg=scrnWidthPix/scrnWidthDeg; % # pixels per degree
    screenCenter(1)=winRect(1)+scrnWidthPix/2; % center of screen
    screenCenter(2)=winRect(2)+scrnHeightPix/2;
    imSize=pixPerDeg*faceSize; % standard deviation of Gaussian depth blob
    
    Screen('TextSize',windowPtr,48); % put some instructions on screen, specify font parameters
    Screen('TextFont',windowPtr,'Arial');
    textToObserver=sprintf('Screen Parameters %d by %d pixels at %3.2f Hz. Click on cells containing different faces. Loading...', winRect(3)-winRect(1),winRect(4)-winRect(2), frameRate);
    Screen('DrawText', windowPtr, textToObserver, 100, 100, 0, meanLUT);
    Screen('TextSize',windowPtr,48);                               % put some instructions on screen
    
    nPCs=length(targPC1); % how many SFs to test
    slopeEsts=2*ones(size(targPC1)); % start with a rough estimate of slope
    threshStd=2*ones(size(targPC1)); % start with a rough estimate of threshold range
    slopeStd=2*ones(size(targPC1)); % start with a rough estimate of slope range
    
    srcRectIm=[0 0 imSize imSize];
    [srcX,srcY]=meshgrid(1:imSize,1:imSize); % pixel co-ordinates of veridical image before depth offsets
    xOffset=imSize/2;
    cellRect=[0 0 imSize*2 imSize];
    nextSizeDeg=3.5; % size of next chart button
    nextRect=CenterRectOnPoint([0 0 nextSizeDeg*pixPerDeg nextSizeDeg*pixPerDeg],winRect(3)-nextSizeDeg*pixPerDeg,screenCenter(2)); % halfway down the screen, away from the left

    nCells=nRows*nCols;
    [xFactor, yFactor]=meshgrid((1:nCols)-nCols/2-0.5,(1:nRows)-nRows/2-0.5); % relative distance from center of screen in multiples of stim size
    trialSeed=zeros(nTrials, nPCs);
    
    for trialNo=1:nTrials
        for pcNum=1:nPCs
            trialRecord(pcNum,trialNo) = struct('trialSeed',0,'targLevelPerLoc',[],'stimYes',[], 'nTargs', [],'targLocs',[],'targLevel', [],'targSFreq', [], 'centerVal1', [], 'centerVal2', []); % start empty array
        end
    end
    
    %% run experiment
    minTargs=round(nRows*nCols*0.66); % min 2/3 full cells,1/3 empty
    maxTargs=nRows*nCols-1; % at least ones empty cell
    tcount = 1; %trial counter
    
    tic;
    for trialNo=1:nTrials % work through all trials
        
        for pcNum=1:nPCs % work though the list of SFs
            
            if trialNo>1 % there has been at least 1 trial
                testLevel=[]; % set blank list of all levels for this SF
                sYes=[]; % blank list of all stimuli seen for this SF
                for trialSoFar=1:nTrials % work through all trials
                    testLevel=[testLevel trialRecord(pcNum, trialSoFar).targLevelPerLoc]; % concatenate all levels
                    sYes=[sYes trialRecord(pcNum,trialSoFar).stimYes]; % concatenate all Yes responses
                end
                myData=arrangeData(testLevel, sYes);
                tLevel=myData(:,1);
                pYes=myData(:,4);
                errEst=myData(:,8);
                tLevel(1)=tLevel(2)/2;
                fitObj=fitDPrimeSensFunc( tLevel, pYes, errEst, 0 );
                threshEsts(pcNum)=fitObj.thresh;
                slopeEsts(pcNum)=fitObj.slope;
                %                 [threshEsts(sfNum), threshStd(sfNum)]=posteriorDistribution(threshEsts(sfNum), threshStd(sfNum), fitobject.thresh, (ci(1,1)-ci(1,2))/2, length(testLevel) ); % update current threshold estimate
                %                 [slopeEsts(sfNum), threshStd(sfNum)]=posteriorDistribution(slopeEsts(sfNum), threshStd(sfNum), fitobject.slope, (ci(2,1)-ci(2,2))/2, length(testLevel) ); % update current slope estimate
            end
            
            minTestLevel=findDPrimeLevel(threshEsts(pcNum),slopeEsts(pcNum),0.1,[0.1 5],'Log'); % find level with d' closest to near zero - low bound minimum disparity (0.5 pixels), high bound maximum disparity (0.5 sigma)
            maxTestLevel=min([100 findDPrimeLevel(threshEsts(pcNum),slopeEsts(pcNum),4.5,[0.5 5],'Log')]); % find level with d' closest to max(5), cap at 100% contrast
            
            Screen('Flip', windowPtr); % clear screen
            
            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(pcNum,trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(pcNum,trialNo).trialSeed); % seed the random number generator
            
            trialRecord(pcNum,trialNo).nTargs=minTargs-1+Randi(maxTargs-minTargs); % random # targets between min and max
            trialRecord(pcNum,trialNo).targLocs=Shuffle(1:nCells); % pick random target locations - first nTargs locations, rest are blank
%             trialRecord(pcNum,trialNo).targLevel=logspace(log10(minTestLevel), log10(maxTestLevel), trialRecord(pcNum,trialNo).nTargs); % log spaced range of contrasts between min and max
            trialRecord(pcNum,trialNo).targLevel=linspace(minTestLevel, maxTestLevel, trialRecord(pcNum,trialNo).nTargs); % log spaced range of contrasts between min and max
            trialRecord(pcNum,trialNo).targLevelPerLoc=zeros(1,nCells); % fill all cell locations with zero contrast
            trialRecord(pcNum,trialNo).targLevelPerLoc(trialRecord(pcNum,trialNo).targLocs(1:trialRecord(pcNum,trialNo).nTargs))=trialRecord(pcNum,trialNo).targLevel; % fill target locations with corresponding contrast
            
            centerPCVal1=centerPCValMu+centerPCValSD*randn(1, nCells); % generate new rndom level of center structur values
            centerPCVal2=centerPCValMu+centerPCValSD*randn(1, nCells); % generate new rndom level of center structur values

            trialRecord(pcNum,trialNo).centerVal1=centerPCVal1;
            trialRecord(pcNum,trialNo).centerVal2=centerPCVal2;
            
            stimSeen=-1*ones(1,nCells); % set up response grid - default all no target seen
            
            xStimCenter=screenCenter(1)+xFactor*imSize*2; % x center of each cell on threen for this size target
            yStimCenter=screenCenter(2)+yFactor*imSize; % y center of each cell on threen for this size target
            xStimCenter(:,1) = xStimCenter(:,1)-imSize/2; xStimCenter(:,3) = xStimCenter(:,3)+imSize/2; %leave space between cells
            yStimCenter(1,:) = yStimCenter(1,:)-imSize/2; yStimCenter(3,:) = yStimCenter(3,:)+imSize/2; %leave space between cells
                       
            for cellNum=1:nCells % work through each cell
                shape_coords = normrnd(0, 1, [199,1]); % create a random vector for the 199-D structural PCA coordinates
                tex_coords = normrnd(0, 1, [199,1]); % ...and the same for the 199-D texture PCA coordinates

                % render and draw first face
                shape_coords(targPC1(pcNum))=centerPCVal1(cellNum)+trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC1 to mid level + threshold deviation
                shape_coords(targPC2(pcNum))=centerPCVal2(cellNum)+trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC2 to mid level + threshold deviation
                
                mesh1  = coef2object(shape_coords, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
                tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
                h = figure(1);
                rp.width=400; % rp = render parameters for BFM
                rp.height=400;
                if tilt == 1
                    rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between +5 and -5 degrees
                end
                display_face(mesh1, tex1, model.tl, rp);
                set(gcf, 'Color', [ 0.5 0.5 0.5 ]);
                f1 = getframe;
                img1 = f1.cdata;
                destRect=CenterRectOnPoint(srcRectIm, xStimCenter(cellNum)+xOffset, yStimCenter(cellNum)); % add this translation to dest rect
                stimTex=Screen('MakeTexture', windowPtr,img1);
                Screen('DrawTexture', windowPtr, stimTex, [],destRect);
                Screen('Close', stimTex);
                
                % render and draw second face
                shape_coords(targPC1(pcNum))=centerPCVal1(cellNum)-trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC1 to mid level + threshold deviation
                shape_coords(targPC2(pcNum))=centerPCVal2(cellNum)-trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC2 to mid level + threshold deviation
                
                mesh1  = coef2object(shape_coords, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
                tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
                h = figure(1);
                rp.width=400; % rp = render parameters for BFM
                rp.height=400;
                if tilt == 1
                    rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between +5 and -5 degrees
                end
                display_face(mesh1, tex1, model.tl, rp)
                set(gcf, 'Color', [ 0.5 0.5 0.5 ])
                f1 = getframe;
                img1 = f1.cdata;
                destRect=CenterRectOnPoint(srcRectIm, xStimCenter(cellNum)-xOffset, yStimCenter(cellNum)); % add this translation to dest rect
                stimTex=Screen('MakeTexture', windowPtr,img1);
                Screen('DrawTexture', windowPtr, stimTex, [],destRect);
                Screen('Close', stimTex);
                
                destRect=CenterRectOnPoint(cellRect, xStimCenter(cellNum), yStimCenter(cellNum));
                Screen('FrameRect', windowPtr,1,destRect); % draw grey line around cells
            end
            
            DrawFormattedText(windowPtr, sprintf('Trial %d', tcount)); %trialnum
            Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw next button
            Screen('TextBackgroundColor', windowPtr, [32 255 64]); %green background for button
            DrawFormattedText(windowPtr, sprintf('%s', 'next'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
            Screen('TextBackgroundColor', windowPtr, [0,0,0,0]); %transparent background for trial num
            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
            
            if saveImage
                myImfileName=sprintf('C:\\Users\\bexlab\\Dropbox\\Kerri_Walter\\BasalFace\\%dTrial%d.jpg', pcNum,trialNo);
                myImage=Screen('GetImage', windowPtr);
                imwrite(myImage,myImfileName);
            end
            
            ShowCursor('Arrow');
            clickedExit=0;
            while clickedExit==0
                [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
                while any(buttons) % if already down, wait for release
                    [mx,my,buttons] = GetMouse(windowPtr);
                end
                while ~any(buttons) % wait for new press
                    [mx,my,buttons] = GetMouse(windowPtr);
                end
                
                    if mx > nextRect(1) % observer clicked finished word
                        clickedExit=1;
                    else
                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                    [~,respNum]=min(mouseDistFromEachBox(:));
                    stimSeen(respNum)=stimSeen(respNum)*-1;
                    if stimSeen(respNum)==1 % seen
                        Screen('FrameRect', windowPtr, [0 255 0],CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2); %
                    else % reversed from seen back to unseen
                         Screen('FrameRect', windowPtr, 0,CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2);
                    end
                    Screen('Flip', windowPtr, [], 1); % show seen ring and don't erase buffer
                end
            end
            stimSeen(stimSeen==-1)=0; % convert -1 to 0
            trialRecord(pcNum,trialNo).stimYes=stimSeen; % mark this cell as target seen
            
            if mod(tcount,10) == 0 %if trial is a multiple of 10 (every 10 trials, break)
                Screen('Flip', windowPtr); % clear screen
                textToObserver=sprintf('Break! Press SPACE to continue');
                Screen('DrawText', windowPtr, textToObserver, 100, 100, 0, meanLUT);
                Screen('TextSize',windowPtr,48); %break screen
                Screen('Flip', windowPtr); %show
                while 1 % wait indefinitely (until loop is exited)
                    [keyIsDown, secs, keyCode] = KbCheck;
                    if keyCode(KbName('Space')) %break out of while loop, continue script to next trial
                        break;
                    end
                end
            elseif tcount == nTrials*nPCs %unless it's the last trial, just quit
                %do nothing
            end
            
            tcount= tcount+1; %increase trial counter
            
            Screen('Flip', windowPtr); % clear screen
            textToObserver=sprintf('Loading...');
            Screen('DrawText', windowPtr, textToObserver, 100, 100, 0, meanLUT);
            Screen('TextSize',windowPtr,48); %loading screen
            
        end % end SFs loop
    end % end trials loop
    expDuration=toc;
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'findFaces.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sfindFaces_combos.mat', testSName);       % file path to Mat files
    save(dataFile,'trialRecord','expDuration','targPC1','targPC2','threshEsts','faceSize','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');
    
    thresholdEst=zeros(1, nPCs);
    fitLB=zeros(1, nPCs);
    fitUB=zeros(1, nPCs);
    for pcNum=1:nPCs % process data foreach Spatial Frequency
        figure(pcNum); % create new figure
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
        fitLB(pcNum)=err95ci(1,2);
        fitUB(pcNum)=err95ci(2,2);
    end
    fitLB(isnan(fitLB))=max(fitLB); % fix any nans
    fitUB(isnan(fitUB))=max(fitUB);
    figure();
    errorbar(targPC1,thresholdEst, fitLB, fitUB );
    xlabel('Principal Component #'); % label the x axis
    ylabel('Threshold (std)'); % label y axis
    title(sprintf('%s Face Discrimination Threshold', sName));

    save(dataFile,'trialRecord','thresholdEst','fitLB','fitUB','expDuration','targPC1','targPC2','threshEsts','faceSize','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');
catch
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'findFacesFailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sfindFaces_combosFailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'trialRecord','expDuration','targPC1','targPC2','threshEsts','faceSize','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');   
    ShowCursor;
end




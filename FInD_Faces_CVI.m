% FInD Faces
% Program to measure discrimination thresholds for Faces, based on Basel face database
% A sequence of grids containing pairs of faces is presented
% Target face pairs shares 198 random PCs, and differs along 1 PC of interest (PCs 1-10 are measured here)
% Distracter face pairs shares 199 random PCs
% Observer clicks on all cells contaiing different faces
% Staircase adjusts the PC difference on successive trials, based on d' sensitivity analysis
%
% HIDDEN-CELL VERSION
%
% 2022  PJB

clear all;

addpath('C:\Users\bexlab\Dropbox\Kerri_Walter\BasalFace')
addpath('C:\Users\bexlab\Dropbox\Kerri_Walter\FInD functions')
% addpath(genpath('C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace'))
% addpath(genpath('C:\Users\Kerri\Dropbox\Kerri_Walter\FInD functions'))
% addpath(genpath('C:\KerriWork\BasalFace'))

whichScreen=max(Screen('Screens')); % use highest display # for stimuli
% whichScreen=2;
meanLUT=127;
% read in testing parameters:
prompt = {'subject initials','Targ Size (deg)','# rows','# columns','# Trials','Screen Width (cm)','Viewing Distance (cm)'};
dlg_title = 'FInD Depth';
num_lines = 1;
def = {'XX', '6', '3', '3', '4', '50', '50'};
answer = inputdlg(prompt,dlg_title,num_lines,def); % read in parameters from GUI
sName=char(answer(1,1));
%targPC=str2num(char(answer(2,1)));
targPC=[1:10];
%threshEsts=str2num(char(answer(3,1)));
threshEsts=repmat(4,1,length(targPC));
faceSize=str2num(char(answer(2,1)));
nRows=str2num(char(answer(3,1)));
nCols=str2num(char(answer(4,1)));
nTrials=str2num(char(answer(5,1)));
scrnWidthCm=str2num(char(answer(6,1))); % screen width (cm)
viewDistCm=str2num(char(answer(7,1))); % observer's distance from screen

rng('default');
% bfm_dir = 'C:\Users\bexlab\Dropbox\Kerri_Walter\BasalFace';
% bfm_dir = 'C:\Users\Kerri\Dropbox\Kerri_Walter\BasalFace';
% cd(bfm_dir); % switch to Basle Face Model Directory
[model msz] = load_model(); % load BFM
rp = defrp; % rp = render parameters for BFM
rp.phi = 0; % render front-facing view

centerPCValMu=0; % TODO - this controls the mid point of the discrimination vector
centerPCValSD=0.5; % TODO - this controls the standard deviation of the mid point of the discrimination vector

saveImage=0;
q=0; %quit variable

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
    Screen('TextSize',windowPtr,48); % put some instructions on screen
    
    nPCs=length(targPC); % how many SFs to test
    slopeEsts=2*ones(size(targPC)); % start with a rough estimate of slope
    threshStd=2*ones(size(targPC)); % start with a rough estimate of threshold range
    slopeStd=2*ones(size(targPC)); % start with a rough estimate of slope range
    
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
            trialRecord(pcNum,trialNo) = struct('trialSeed',0,'targLevelPerLoc',[],'stimYes',[], 'nTargs', [],'targLocs',[],'targLevel', [],'targSFreq', [], 'centerVal', []); % start empty array
        end
    end
    
    %% run experiment
    minTargs=round(nRows*nCols*0.66); % min 2/3 full cells,1/3 empty
    maxTargs=nRows*nCols-1; % at least ones empty cell
    tcount = 1; %trial counter
    
    tic;
    for trialNo=1:nTrials % work through all trials
        
        for pcNum=1:nPCs % work though the list of SFs
            my_coordinate=[]; mx_coordinate=[]; buttons_coordinate=[]; %create new mouse position variables each trial
            
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
            %             vbl=Screen('Flip', windowPtr); % clear screen
            %             buttonPressTime=vbl;
            
            SetMouse(screenCenter(1), screenCenter(2), windowPtr); % move mouse to screen center
            trialRecord(pcNum,trialNo).trialSeed=round(sum(100*clock)); % use current time to generate new seed number for random number generation so that all trial parameters can be reproduced if necessary
            rng(trialRecord(pcNum,trialNo).trialSeed); % seed the random number generator
            
            trialRecord(pcNum,trialNo).nTargs=minTargs-1+Randi(maxTargs-minTargs); % random # targets between min and max
            trialRecord(pcNum,trialNo).targLocs=Shuffle(1:nCells); % pick random target locations - first nTargs locations, rest are blank
            %             trialRecord(pcNum,trialNo).targLevel=logspace(log10(minTestLevel), log10(maxTestLevel), trialRecord(pcNum,trialNo).nTargs); % log spaced range of contrasts between min and max
            trialRecord(pcNum,trialNo).targLevel=linspace(minTestLevel, maxTestLevel, trialRecord(pcNum,trialNo).nTargs); % log spaced range of contrasts between min and max
            trialRecord(pcNum,trialNo).targLevelPerLoc=zeros(1,nCells); % fill all cell locations with zero contrast
            trialRecord(pcNum,trialNo).targLevelPerLoc(trialRecord(pcNum,trialNo).targLocs(1:trialRecord(pcNum,trialNo).nTargs))=trialRecord(pcNum,trialNo).targLevel; % fill target locations with corresponding contrast
            
            centerPCVal=centerPCValMu+centerPCValSD*randn(1, nCells); % generate new rndom level of center structur values
            trialRecord(pcNum,trialNo).centerVal=centerPCVal;
            
            stimSeen=-1*ones(1,nCells); % set up response grid - default all no target seen
            
            xStimCenter=screenCenter(1)+xFactor*imSize*2; % x center of each cell on threen for this size target
            yStimCenter=screenCenter(2)+yFactor*imSize; % y center of each cell on threen for this size target
            xStimCenter(:,1) = xStimCenter(:,1)-imSize/2; xStimCenter(:,3) = xStimCenter(:,3)+imSize/2; %leave space between cells
            yStimCenter(1,:) = yStimCenter(1,:)-imSize/2; yStimCenter(3,:) = yStimCenter(3,:)+imSize/2; %leave space between cells
            
            for cellNum=1:nCells % work through each cell
                shape_coords = normrnd(0, 1, [199,1]); % create a random vector for the 199-D structural PCA coordinates
                tex_coords = normrnd(0, 1, [199,1]); % ...and the same for the 199-D texture PCA coordinates
                
                % render and draw first face
                shape_coords(targPC(pcNum))=centerPCVal(cellNum)+trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC to mid level + threshold deviation
                
                mesh1  = coef2object(shape_coords, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
                tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
                h = figure(1);
                rp.width=400; % rp = render parameters for BFM
                rp.height=400;
                rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between +5 and -5 degrees
                display_face(mesh1, tex1, model.tl, rp);
                set(gcf, 'Color', [ 0.5 0.5 0.5 ]);
                f1 = getframe;
                img1(cellNum,1) = {f1.cdata};
                destRect(cellNum,1)={CenterRectOnPoint(srcRectIm, xStimCenter(cellNum)+xOffset, yStimCenter(cellNum))}; % add this translation to dest rect
                %                 stimTex=Screen('MakeTexture', windowPtr,img1);
                %                 Screen('DrawTexture', windowPtr, stimTex, [],destRect);
                %                 Screen('Close', stimTex);
                
                % render and draw second face
                shape_coords(targPC(pcNum))=centerPCVal(cellNum)-trialRecord(pcNum,trialNo).targLevelPerLoc(cellNum)/2; % set target PC to mid level + threshold deviation
                mesh1  = coef2object(shape_coords, model.shapeMU, model.shapePC, model.shapeEV ); % Convert into vertex space
                tex1 = coef2object(tex_coords, model.texMU,  model.texPC,   model.texEV); % Convert into texture RGB space
                h = figure(1);
                rp.width=400; % rp = render parameters for BFM
                rp.height=400;
                rp.phi=deg2rad(-5)+deg2rad(10)*rand(); %random pitch between +5 and -5 degrees
                display_face(mesh1, tex1, model.tl, rp)
                set(gcf, 'Color', [ 0.5 0.5 0.5 ])
                f1 = getframe;
                img1(cellNum,2) = {f1.cdata};
                destRect(cellNum,2)={CenterRectOnPoint(srcRectIm, xStimCenter(cellNum)-xOffset, yStimCenter(cellNum))}; % add this translation to dest rect
                %                 stimTex=Screen('MakeTexture', windowPtr,img1);
                %                 Screen('DrawTexture', windowPtr, stimTex, [],destRect);
                %                 Screen('Close', stimTex);
                
                destRectFull=CenterRectOnPoint(cellRect, xStimCenter(cellNum), yStimCenter(cellNum));
                %                 rects = Screen('MakeTexture', windowPtr, [127.5 127.5 127.5]);
                %                 Screen('DrawTexture', windowPtr, rects, [],destRect);
                Screen('FrameRect', windowPtr,0,destRectFull); % draw black line around cells
                
            end
            
            %             DrawFormattedText(windowPtr,  'Click on cells where faces are different:', 'center', 100, 0);
            Screen('FillOval',windowPtr,[32 255 64], nextRect);  % draw next button
            DrawFormattedText(windowPtr, sprintf('%s', 'next'), 'center', 'center', 255, [],[],[],[],[],nextRect); %  char(26)
            DrawFormattedText(windowPtr, sprintf('Trial %d', tcount)); %trialnum
            Screen('Flip', windowPtr, [], 1); % show stimulus and don't erase buffer
            chartStart=tic; % time for each chart
            
            if saveImage
                myImfileName=sprintf('C:\\Users\\bexlab\\Dropbox\\Kerri_Walter\\BasalFace\\%dTrial%d.jpg', pcNum,trialNo);
                myImage=Screen('GetImage', windowPtr);
                imwrite(myImage,myImfileName);
            end
            
            ShowCursor('Arrow');
            clickedExit=0;
            pressCheck=0;
            
            while clickedExit==0
                
                [mx,my,buttons] = GetMouse(windowPtr); % wait for mouse button release before processing response
                mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                [~,respClick]=min(mouseDistFromEachBox(:));
                
                if any(buttons) && pressCheck == 0% if already down, wait for release
                    stimSeen(respClick)=stimSeen(respClick)*-1;
                    pressCheck = 1; %only do this once
                end
                
                if any(buttons) % Check if any buttons are being pressed
                    my_coordinate=[my_coordinate;my]; mx_coordinate=[mx_coordinate;mx]; buttons_coordinate=[buttons_coordinate;buttons]; %save continuous mouse position data during button press
                end
                
                while ~any(buttons) % wait for new press
                    [mx,my,buttons] = GetMouse(windowPtr);
                    my_coordinate=[my_coordinate;my]; mx_coordinate=[mx_coordinate;mx]; buttons_coordinate=[buttons_coordinate;buttons]; %save continuous mouse position data while button is not pressed
                    
                    pressCheck=0;
                    
                    mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                    [~,respNum]=min(mouseDistFromEachBox(:));
                    
                    for cellNum=1:nCells % work through each cell
                        if cellNum==respNum % only draw faces where mouse currently is
                            stimTex1=Screen('MakeTexture', windowPtr,img1{respNum,1}); %pull the faces for this cell
                            Screen('DrawTexture', windowPtr, stimTex1, [],destRect{respNum,1}); %draw faces for this cell
                            stimTex2=Screen('MakeTexture', windowPtr,img1{respNum,2}); %pull the faces for this cell
                            Screen('DrawTexture', windowPtr, stimTex2, [],destRect{respNum,2}); %draw faces for this cell
                            
                            if stimSeen(respNum)==1 % seen
                                Screen('FrameRect', windowPtr, [0 255 0],CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2); %
                            else % reversed from seen back to unseen
                                Screen('FrameRect', windowPtr, 0,CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2);
                            end
                            
                        end
                    end
                    Screen('Flip', windowPtr, [], 1); % show faces
                    Screen('FillRect', windowPtr, [127.5 127.5 127.5],CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2); %redraw box
                    if stimSeen(respNum)==1 % seen
                        Screen('FrameRect', windowPtr, [0 255 0],CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2); %
                    else % reversed from seen back to unseen
                        Screen('FrameRect', windowPtr, 0,CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2);
                    end
                end
                
                if mx > nextRect(1) % observer clicked finished word
                    clickedExit=1;
                    %                 else
                    %                     mouseDistFromEachBox=sqrt((mx-xStimCenter).^2+(my-yStimCenter).^2); % calculate distance from each box
                    %                     [~,respNum]=min(mouseDistFromEachBox(:));
                    %                     stimSeen(respNum)=stimSeen(respNum)*-1;
                    %                     if stimSeen(respNum)==1 % seen
                    %                         Screen('FrameRect', windowPtr, [0 255 0],CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2); %
                    %                     else % reversed from seen back to unseen
                    %                          Screen('FrameRect', windowPtr, 0,CenterRectOnPoint(cellRect, xStimCenter(respNum), yStimCenter(respNum)),2);
                    %                     end
                    %                     Screen('Flip', windowPtr, [], 1); % show seen ring and don't erase buffer
                end
                
            end
            stimSeen(stimSeen==-1)=0; % convert -1 to 0
            trialRecord(pcNum,trialNo).stimYes=stimSeen; % mark this cell as target seen
            
            trialRecord(pcNum,trialNo).chartTime=toc(chartStart); %save chart time
            
            dataMouse(pcNum,trialNo).my=my_coordinate;
            dataMouse(pcNum,trialNo).mx=mx_coordinate;
            dataMouse(pcNum,trialNo).buttons=buttons_coordinate;
            
            if tcount == nTrials*nPCs %if it's the last trial, just quit
                %do nothing
            elseif mod(tcount,10) == 0 %if trial is a multiple of 10 (every 10 trials, break)
                Screen('Flip', windowPtr); % clear screen
                textToObserver=sprintf('Break! Press SPACE to continue');
                Screen('DrawText', windowPtr, textToObserver, 100, 100, 0, meanLUT);
                Screen('TextSize',windowPtr,48); %break screen
                Screen('Flip', windowPtr); %show
                while 1 % wait indefinitely (until loop is exited)
                    [keyIsDown, secs, keyCode] = KbCheck;
                    if keyCode(KbName('Space')) %break out of while loop, continue script to next trial
                        break;
                    elseif keyCode(KbName('Q'))
                        q=1; %quit
                        break
                    end
                end
            end
            
            if q==1
                break; %quit
            end
            
            tcount= tcount+1; %increase trial counter
            
            Screen('Flip', windowPtr); % clear screen
            textToObserver=sprintf('Loading...');
            Screen('DrawText', windowPtr, textToObserver, 100, 100, 0, meanLUT);
            Screen('TextSize',windowPtr,48); %loading screen
            
        end % end SFs loop
        
        if q==1
            break; %quit
        end
        
    end % end trials loop
    expDuration=toc;
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'findFaces.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sfindFaces.mat', testSName);       % file path to Mat files
    save(dataFile,'trialRecord','expDuration','targPC','threshEsts','faceSize','dataMouse','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');
    
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
    errorbar(targPC,thresholdEst, fitLB, fitUB );
    xlabel('Principal Component #'); % label the x axis
    ylabel('Threshold (std)'); % label y axis
    title(sprintf('%s Face Discrimination Threshold', sName));
    
    save(dataFile,'trialRecord','thresholdEst','fitLB','fitUB','expDuration','targPC','threshEsts','faceSize','dataMouse','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');
catch
    Screen('CloseAll');
    testSName=sName;    % modify sName if subject already exists
    n=1;
    while exist([testSName,'findFacesFailedRun.mat'],'file') ~= 0
        n=n+1;
        testSName=[sName,num2str(n)];
    end
    dataFile=sprintf('%sfindFacesFailedRun.mat', testSName);       % file path to Mat files
    save(dataFile,'trialRecord','thresholdEst','fitLB','fitUB','expDuration','targPC','threshEsts','faceSize','dataMouse','nRows','nCols','nTrials','scrnWidthCm','viewDistCm','frameRate','winRect');
    ShowCursor;
end




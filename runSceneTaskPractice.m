%% Scene Task Practice Task
%  Training & Practice task for Scene Task experiment. Designed to be run
%  in fMRI with Eye Tracking.
%
%  This script and all related files are created
%  by and for the Dilks Lab and may not be used, reproduced, or replicated
%  without express permission. Please contact Sam Weiller at
%  sam.weiller@gmail.com with any questions.

%% History
%  3.10     11/18/14    Added comments & Eye Tracking

function [SCPRAC] = runSceneTaskPractice(sub, cbl)
%% Start me up
clc
SCPRAC.curDir = cd;
if ~exist('sub', 'var'); SCPRAC.subID = input('\nPlease Enter Your Participant Code #: ', 's'); else SCPRAC.subID = sub; end;
if ~exist('cbl', 'var'); SCPRAC.run =  input('\nPlease Enter The Run #: ', 's'); else  SCPRAC.run = cbl; end;
PATH = fullfile(SCPRAC.curDir, sprintf('SCPRAC_S%d_C%d.mat', SCPRAC.subID, SCPRAC.run));
save(PATH);
if ~exist(PATH, 'file');
    [Path, File] = uigetfile('*.mat', 'Select .MAT with data');
    PATH = fullfile(Path, File);
end
load(PATH);

pause(.2);
fprintf('Scene Task\n');
pause(.1);
fprintf('  Version 3.10\n');
fprintf('  Nov. 18, 2014\n');
pause(.1);
fprintf('Sam Weiller & Greg Adler\n');
pause(.3);
clc;

%% Control Panel
STIMS = [];

stimuliMatFileName = 'sceneTaskPracticeStims1.mat';
load sceneTaskParamsPRACTICE.mat
fprintf('Looking for stimuli...\n')
if exist(stimuliMatFileName, 'file')
    load(stimuliMatFileName);
    fprintf('Stimuli loaded!\n');
else
    fprintf('Please run makeStims first.\n');
    return;
end;

% Determines condition order
designs = [...
    0 1 0 2 0 3 0 1 0 2 0 3;
    0 1 0 2 0 3 0 1 0 2 0 3;
    0 1 0 2 0 3 0 1 0 2 0 3;
    0 1 0 2 0 3 0 1 0 2 0 3;
    0 1 0 2 0 3 0 1 0 2 0 3;
    ];

numStimSets = size(STIMS,2);
imgsPerSet = size(STIMS{1},2);
numBlocks = max(max(size(designs)));
imagesPerBlock = 12;
totalImages = 6;
fixationTime = 10;
KbName('UnifyKeyNames');
% FOR BEHAVORIAL, USE THESE
% choice1key = KbName('7&');
% choice2key = KbName('8*');
% choice3key = KbName('9(');
% % FOR SCANNER, USE THESE
choice1key = KbName('b');
choice2key = KbName('y');
choice3key = KbName('g');
escapeKey = KbName('q');
triggerKey = KbName('t');
screens = Screen('Screens');
screenNumber = min(screens);

screenWidth = 520; 
viewingDistance = 660; 
res = Screen('Resolution', screenNumber);
resWidth = res.width;
visualAngle = 18;
PPD = tand(.5).*2.*viewingDistance.*(resWidth./screenWidth);
visualAngle = PPD*visualAngle;

originalStimSize.horizontal = 1040;
originalStimSize.vertical   = 784;

stimSize.horizontal = visualAngle;
stimSize.vertical = (stimSize.horizontal*originalStimSize.vertical)/originalStimSize.horizontal;

stimPresentTime = .5;
ISITime = 1.3;
trialLength = stimPresentTime + ISITime;
taskPresentationTime = 4;



UserAns = 0;
conditionOrder = designs(cbl, :);

taskNames = {'NAVIGATE\n\n\nLeft     Center     Right', 'CATEGORIZE\n\n\nBedroom     Kitchen     Living Room', 'COLOR\n\n\n 1     2     3'};

%% PTB Setup
Screen('Preference', 'SkipSyncTests', 2);
[w, ~, xMid, yMid] = startPTB(screenNumber, 1, [50 50 50]);
HideCursor;

fixBoxSize = 2*PPD; % 2 degrees
fixBox = [round(xMid - fixBoxSize/2), round(yMid - fixBoxSize/2), round(xMid + fixBoxSize/2), round(yMid + fixBoxSize/2)];

%% Eyelink Setup
el = EyelinkInitDefaults(w);

if ~EyelinkInit(dummymode)
    fprintf('Eyelink Init Aborted.\n');
    Eyelink('Shutdown');
    return;
end;

[v, vs] = Eyelink('GetTrackerVersion');
fprintf('Running Experiment on a "%s" tracker.\n', vs);


i = Eyelink('Openfile', edfFile);

if i ~=0
    fprintf('Cannot create EDF file "%s"', edffilename);
    Eyelink('Shutdown');
    return;
end;

Eyelink('command', 'add_file_preamble_text "Recorded by EyelinkToolbox. Script by SKW"');

[width, height] = Screen('WindowSize', w);

Eyelink('command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);

Eyelink('command', 'calibration_type = HV9');
Eyelink('command', 'saccade_velocity_threshold = 35');
Eyelink('command', 'saccade_acceleration_threshold = 9500');

Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');

if ( Eyelink('IsConnected') ~= 1 && ~dummymode )
    Eyelink('Shutdown');
    return;
end;

moveOn = 0;

el.backgroundcolour = [128 128 128];
el.foregroundcolour = [0 0 0];

EyelinkDoTrackerSetup(el);


%% Create Stimuli & Preallocate
tex = cell(numStimSets, 1);
ANSMAT = cell(numBlocks, 1);

for set = 1:numStimSets
    for img = 1:imgsPerSet
        % Making the cell array
        tex{set}{img} = Screen('MakeTexture', w, STIMS{set}{img});
    end;
end;

DrawFormattedText(w, 'Waiting for trigger...', 'center', 'center');
Screen('Flip', w);
trigger(triggerKey);
Screen('Flip', w);

%% Main Loop
Eyelink('Command', 'set_idle_mode');
Eyelink('Command', 'clear_screen 0')

WaitSecs(0.05);  
Eyelink('StartRecording');    
WaitSecs(0.1);

eyeLinkTrial = 1;

expStart = GetSecs;

cpuTimeExpected = expStart;
realTimeExpected = 0;
SCPRAC.trialCount = zeros(1,numBlocks);

for block = 1:numBlocks
    blockStart = GetSecs;
    timeLogger.block(block).blockStart = GetSecs - expStart;
    timeLogger.block(block).conditionN = conditionOrder(block);
    
    Eyelink('message', 'TRIALID %d', eyeLinkTrial); % Start Eyelink Trial
    
    if conditionOrder(block) == 0
        realTimeExpected = realTimeExpected + fixationTime;
        cpuTimeExpected = cpuTimeExpected + fixationTime;
        
        fixate(w);
        while GetSecs <= cpuTimeExpected - taskPresentationTime
            % Waiting for fixation
        end;
        
        if block ~= numBlocks
            DrawFormattedText(w, taskNames{conditionOrder(block + 1)}, 'center', 'center');
            Screen('Flip', w);
            % Displays button mappings for all but final block.
        end;
        
        while GetSecs <= cpuTimeExpected
            % finish waiting
        end;
    else
        tLstart = GetSecs;
        
        imageMatrix = randsample(totalImages, imagesPerBlock, 1);
        
        correctCount = 0; 
        incorrectCount = 0;
        trial = 1; 
        
        while correctCount < 10
            SCPRAC.trialCount(block) = SCPRAC.trialCount(block) + 1;
            
            if trial == imagesPerBlock
                trial = 1;
            end;
            
            touch = 0;
            
            timeLogger.block(block).trial(trial).start = GetSecs-tLstart;
            loggingIsDone = 0;            
            
            trialStart = GetSecs;
            stimEnd = trialStart + stimPresentTime;
            cpuTimeExpected = cpuTimeExpected + trialLength;
            realTimeExpected = realTimeExpected + trialLength;            
            
            %Draw cell array
            Screen('DrawTexture', w, tex{TRIMAT(imageMatrix(trial), 3)}{TRIMAT(imageMatrix(trial), 2)}, [], [xMid-(stimSize.horizontal/2) yMid-(stimSize.vertical/2) xMid+(stimSize.horizontal/2) yMid+(stimSize.vertical/2)]);
            fixation(w);
            Screen('Flip', w);
            
            while GetSecs <= stimEnd   %checks for keypress during stim presentation
                [touch, ~, keyCode] = KbCheck(-1);
                if touch && ~keyCode(triggerKey)
                    UserAns = find(keyCode);
                    break;
                end;
            end;
            
            while GetSecs <= stimEnd
                % Wait remaining time
            end;
            
            fixation(w);
            Screen('Flip', w);
            timeLogger.block(block).trial(trial).imageEnd = GetSecs-tLstart;
            
            if touch == 0
                while GetSecs < cpuTimeExpected-.1   %checks for keypress in fixation immediately following stim pres up until next stim pres.
                    [touch, ~, keyCode] = KbCheck(-1);
                    if touch && ~keyCode(triggerKey)
                        UserAns = find(keyCode);
                        break;
                    end;
                end;
            end;
            
            while GetSecs <= cpuTimeExpected
                if ~loggingIsDone
                    ANSMAT{block}(trial, 1) = conditionOrder(block); % Task number
                    ANSMAT{block}(trial, 2) = imageMatrix(trial); %trial number
                    ANSMAT{block}(trial, 3) = TRIMAT(imageMatrix(trial), 4); % direction
                    ANSMAT{block}(trial, 4) = TRIMAT(imageMatrix(trial), 3); % category
                    ANSMAT{block}(trial, 5) = TRIMAT(imageMatrix(trial), 5); % color
                    
                    fprintf('User Response: %d\n', UserAns);

                    if length(UserAns) > 1
                        ANSMAT{block}(trial, 6) = 0;
                    else
                        switch UserAns % User Response
                            case choice1key
                                ANSMAT{block}(trial, 6) = 1;
                            case choice2key
                                ANSMAT{block}(trial, 6) = 2;
                            case choice3key
                                ANSMAT{block}(trial, 6) = 3;
                            case escapeKey
                                SCPRAC.ANSMAT = ANSMAT;
                                save(PATH, 'SCPRAC', 'timeLogger');
                                Screen('CloseAll');
                                return;
                            otherwise
                                ANSMAT{block}(trial, 6) = 0;
                        end;
                    end;
                    
                    fprintf('ANSMAT Log: %d\n', ANSMAT{block}(trial, 6));
                    
                    switch conditionOrder(block) % correct/incorrect
                        case 1 % Direction
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 3) == ANSMAT{block}(trial, 6);
                        case 2 % Category
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 4) == ANSMAT{block}(trial, 6);
                        case 3 % Color
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 5) == ANSMAT{block}(trial, 6);
                    end;
                    
                    switch ANSMAT{block}(trial, 7) == 1
                        case 1
                            correctCount = correctCount + 1;
                        case 0
                            incorrectCount = incorrectCount + 1;

                            if incorrectCount >= 3
                                correctCount = 0;
                                incorrectCount = 0;
                            end;
                    end;
                    
                    UserAns = 0;
                    loggingIsDone = 1;
                end;
            end;
            
            timeLogger.block(block).trial(trial).end = GetSecs-tLstart;
            timeLogger.block(block).trial(trial).imageLength = timeLogger.block(block).trial(trial).imageEnd - timeLogger.block(block).trial(trial).start;
            timeLogger.block(block).trial(trial).blankLength = timeLogger.block(block).trial(trial).end - timeLogger.block(block).trial(trial).imageEnd;
            timeLogger.block(block).trial(trial).trialLength = timeLogger.block(block).trial(trial).end - timeLogger.block(block).trial(trial).start;
            
            trial= trial+1;
            SCPRAC.ANSMAT = ANSMAT;
            save(PATH, 'SCPRAC', 'timeLogger');
        end;

    end;
    
    SCPRAC.ANSMAT = ANSMAT;
    timeLogger.block(block).blockLength = GetSecs - blockStart;
    timeLogger.block(block).blockEnd = realTimeExpected;
    fprintf('Block Time: %1.4f\n', timeLogger.block(block).blockLength);    
    save(PATH, 'SCPRAC', 'timeLogger');
end;

%% Shutdown Procedures
ShowCursor;
Screen('CloseAll');

function [w, rect, xc, yc] = startPTB(screenNumber, oGl, color)

if nargin == 0
    oGl = 0;
    color = [0 0 0];
elseif nargin == 1;
    color = [0 0 0];
end;

[w, rect] = Screen('OpenWindow', screenNumber, color);
xc = rect(3)/2;
yc = rect(4)/2;

if oGl == 1
    AssertOpenGL;
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
end;

function fixate(w)
Screen('TextSize', w, 40);
DrawFormattedText(w, '+', 'center', 'center', [200 200 200]);
Screen('TextSize', w, 50);
Screen('Flip', w);

function fixation(w)
Screen('TextSize', w, 40);
DrawFormattedText(w, '+', 'center', 'center', [200 200 200]);
Screen('TextSize', w, 50);

function trigger(triggerKey)
KbName('UnifyKeyNames');

go = 0;
while go == 0
    [touch, ~, keyCode] = KbCheck(-1);
    WaitSecs(.0001);
    if touch && keyCode(triggerKey)
        go = 1;
    end;
end;
%% PnP Localizer
%  TO DO
%  -Add options flag. (mat file name, screen max or min, etc.)
%  -Add COV maker
%  -print number of volumes?
%  -cross color
%  -cross during whitespace
%  -button explanation
%  -Task description higher on screen
function [SCTASK] = runSceneTaskET(sub, cbl, acq)
%% Start me up
clc
SCTASK.curDir = cd;
if ~exist('sub', 'var'); SCTASK.subID = input('\nPlease Enter Your Participant Code #: ', 's'); else SCTASK.subID = sub; end;
if ~exist('cbl', 'var'); SCTASK.run =  input('\nPlease Enter The Run #: ', 's'); else  SCTASK.run = cbl; end;
if ~exist('acq', 'var'); SCTASK.acq =  input('\nPlease Enter The Aquisition #: ', 's'); else  SCTASK.acq = acq; end;
PATH = fullfile(SCTASK.curDir, sprintf('SCTASK_S%d_C%d_A%d.mat', SCTASK.subID, SCTASK.run, SCTASK.acq));
save(PATH);
if ~exist(PATH, 'file');
    [Path, File] = uigetfile('*.mat', 'Select .MAT with data');
    PATH = fullfile(Path, File);
end
load(PATH);

pause(.8);
fprintf('Scene Task\n');
pause(.7);
fprintf('  Version 0.10\n');
fprintf('  Jun. 6, 2014\n');
pause(.4);
fprintf('Sam Weiller & Greg Adler\n');
pause(.5);
clc

%% Control Panel
STIMS = [];

stimuliMatFileName = 'sceneTaskStims1.mat';
load sceneTaskParams.mat
fprintf('Looking for stimuli...\n')
if exist(stimuliMatFileName, 'file')
    load(stimuliMatFileName);
    fprintf('Stimuli loaded!\n');
else
    fprintf('Please run makeStims first.\n');
    return;
end;

% CREATE NEW DESIGN VECTORS
designs = [...
    0 1 0 2 0 3 0 1 0 3 0 2 0 2 0 3 0 1 0 3 0 2 0 1 0;
    0 3 0 1 0 2 0 2 0 1 0 3 0 3 0 1 0 2 0 2 0 1 0 3 0;
    0 2 0 3 0 1 0 1 0 2 0 3 0 3 0 2 0 1 0 1 0 3 0 2 0;
    0 3 0 2 0 1 0 3 0 1 0 2 0 2 0 1 0 3 0 1 0 2 0 3 0;
    0 1 0 3 0 2 0 2 0 3 0 1 0 1 0 3 0 2 0 2 0 3 0 1 0;
    ];

fixCovariate = max(max(designs)) + 1;          % Defines a non-zero covariate number for fixation covariate files.
numStimSets = size(STIMS,2);
imgsPerSet = size(STIMS{1},2);
numBlocks = max(max(size(designs)));
imagesPerBlock = 16;
totalImages = 36;
% numberOfTargets = 2;
fixationTime = 10;
KbName('UnifyKeyNames');
% FOR BEHAVORIAL, USE THESE
choice1key = KbName('7&');
choice2key = KbName('8*');
choice3key = KbName('9(');
% % FOR SCANNER, USE THESE
% choice1key = KbName('r');
% choice2key = KbName('g');
% choice3key = KbName('b');
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

dummymode = 0;

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

taskNames = {'NAVIGATE\n\n\nLeft     Center     Right', 'CATEGORIZE\n\n\nBedroom     Kitchen     Living Room', 'COLOR\n\n\nOne     Two     Three'};

%% Eyelink Params
ET = 1;
pref_eye = 1; % 0 is left, 1 is right, 2 is both
dummymode = 0;

prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
dlg_title = 'Create EDF File';
num_lines = 1;
def = {'DEMO'};
answer = inputdlg(prompt, dlg_title, num_lines, def);
edfFile = answer{1};
fprintf('EDFFile: %s\n', edfFile);

%% PTB Setup & EyeLink Setup
Screen('Preference', 'SkipSyncTests', 2);
[w, ~, xMid, yMid] = startPTB(screenNumber, 1, [50 50 50]);
HideCursor;

% stimSize.rect = [xc-(stimSize.horizontal/2) yc-(stimSize.vertical/2) xc+(stimSize.horizontal/2) yc+(stimSize.vertical/2)];

fixBoxSize = 5*PPD; % 2 degrees
fixBox = [round(xMid - fixBoxSize/2), round(yMid - fixBoxSize/2), round(xMid + fixBoxSize/2), round(yMid + fixBoxSize/2)];

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

el.backgroundcolour = [128 128 128];
el.foregroundcolour = [0 0 0];

EyelinkDoTrackerSetup(el);

%% Create Stimuli & Preallocate
tex = cell(numStimSets, 1);
ANSMAT = cell(numBlocks, 1);

for set = 1:numStimSets
    for img = 1:imgsPerSet
        %Making the cell array
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

for block = 1:numBlocks
    blockStart = GetSecs;
    timeLogger.block(block).blockStart = GetSecs - expStart;
    timeLogger.block(block).conditionN = conditionOrder(block);

    eye_used = pref_eye;
    
    if conditionOrder(block) == 0
        
        Eyelink('message', 'TRIALID %d', eyeLinkTrial);
        realTimeExpected = realTimeExpected + fixationTime;
        cpuTimeExpected = cpuTimeExpected + fixationTime;
        
        
        fixate(w);
        while GetSecs <= cpuTimeExpected - taskPresentationTime
            % Waiting for fixation
        end;
        
        if block ~= numBlocks
            DrawFormattedText(w, taskNames{conditionOrder(block + 1)}, 'center', 'center');
            Screen('Flip', w);
        end;
        
        while GetSecs <= cpuTimeExpected
            % finish waiting
        end;
        WaitSecs(0.1);
        Eyelink('StopRecording');
        Eyelink('message', 'TRIAL_RESULT 0');
        eyeLinkTrial = eyeLinkTrial + 1;
    else
        
        tLstart = GetSecs;
        
        imageMatrix = randsample(totalImages, imagesPerBlock);
        
        for trial = 1:imagesPerBlock
            Eyelink('Command', 'set_idle_mode');
            WaitSecs(0.05);  
            Eyelink('StartRecording');    
            WaitSecs(0.1);
            Eyelink('message', 'TRIALID %d', eyeLinkTrial);
            Eyelink('message', '!V CLEAR 128 128 128');
            Eyelink('command', 'record_status_message "TRIAL %d / %d"', eyeLinkTrial, numBlocks*(imagesPerBlock+1));
            WaitSecs(0.05);
            Eyelink('message', '!V DRAWBOX 0 0 0 %d %d %d %d', fixBox(1), fixBox(2), fixBox(3), fixBox(4));
            Eyelink('message', '!V IAREA RECTANGLE 1 %d %d %d %d FIXAREA', fixBox(1), fixBox(2), fixBox(3), fixBox(4));
            
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
            
%             Eyelink('Message', '!V IMGLOAD CENTER ../images/%s %d %d %d %d', STIMNAMES{TRIMAT(imageMatrix(trial), 3)}{TRIMAT(imageMatrix(trial), 2)}, round(width/2), round(height/2), stimSize, stimSize);
            Eyelink('message', 'BEGIN IMAGE PRESENTATION');
            Eyelink('Message', '!V IMGLOAD CENTER ./images/%s %d %d %d %d', STIMNAMES{TRIMAT(imageMatrix(trial), 3)}{TRIMAT(imageMatrix(trial), 2)}, round(width/2), round(height/2), round(stimSize.horizontal), round(stimSize.vertical));
            
            
            while GetSecs <= stimEnd   %checks for keypress during stim presentation
                [touch, ~, keyCode] = KbCheck(-1);
                if touch && ~keyCode(triggerKey)
                    UserAns = find(keyCode);
                    break;
                end;
            end;
            
            while GetSecs <= stimEnd
            end;
            
            Eyelink('message', 'BEGIN FIXATION TIME');
            
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
                    Eyelink('message', '!V TRIAL_VAR IMG_NAME %s', STIMNAMES{TRIMAT(imageMatrix(trial), 3)}{TRIMAT(imageMatrix(trial), 2)});
                    ANSMAT{block}(trial, 1) = conditionOrder(block); % Task number
                    ANSMAT{block}(trial, 2) = imageMatrix(trial); %trial number
                    ANSMAT{block}(trial, 3) = TRIMAT(imageMatrix(trial), 4); % direction
                    ANSMAT{block}(trial, 4) = TRIMAT(imageMatrix(trial), 3); % category
                    ANSMAT{block}(trial, 5) = TRIMAT(imageMatrix(trial), 5); % color
                    
                    switch TRIMAT(imageMatrix(trial), 4);
                        case 1
                            Eyelink('message', '!V TRIAL_VAR DIRECTION LEFT');
                        case 2
                            Eyelink('message', '!V TRIAL_VAR DIRECTION CENTER');
                        case 3
                            Eyelink('message', '!V TRIAL_VAR DIRECTION RIGHT');
                        otherwise
                            Eyelink('message', '!V TRIAL_VAR DIRECTION ERROR');
                    end;
                    
                    switch TRIMAT(imageMatrix(trial), 3);
                        case 1
                            Eyelink('message', '!V TRIAL_VAR CATEGORY BEDROOM');
                        case 2
                            Eyelink('message', '!V TRIAL_VAR CATEGORY KITCHEN');
                        case 3
                            Eyelink('message', '!V TRIAL_VAR CATEGORY LIVINGROOM');
                        otherwise
                            Eyelink('message', '!V TRIAL_VAR CATEGORY ERROR');
                    end;
                    
                    switch TRIMAT(imageMatrix(trial), 3);
                        case 1
                            Eyelink('message', '!V TRIAL_VAR COLOR ONE');
                        case 2
                            Eyelink('message', '!V TRIAL_VAR COLOR TWO');
                        case 3
                            Eyelink('message', '!V TRIAL_VAR COLOR THREE');
                        otherwise
                            Eyelink('message', '!V TRIAL_VAR COLOR ERROR');
                    end;
                    
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
                                SCTASK.ANSMAT = ANSMAT;
                                save(PATH, 'SCTASK', 'timeLogger');
                                Screen('CloseAll');
                                return;
                            otherwise
                                ANSMAT{block}(trial, 6) = 0;
                        end;
                    end;
                    
                    fprintf('ANSMAT Log: %d\n', ANSMAT{block}(trial, 6));
                    
                    switch conditionOrder(block) % correct/incorrect
                        case 1 % Direction
                            Eyelink('message', '!V TRIAL_VAR TASK NAVIGATION');
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 3) == ANSMAT{block}(trial, 6);
                        case 2 % Category
                            Eyelink('message', '!V TRIAL_VAR TASK CATEGORIZATION');
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 4) == ANSMAT{block}(trial, 6);
                        case 3 % Color
                            Eyelink('message', '!V TRIAL_VAR TASK COLOR');
                            ANSMAT{block}(trial, 7) = ANSMAT{block}(trial, 5) == ANSMAT{block}(trial, 6);
                    end;
                                
                    UserAns = 0;
                    loggingIsDone = 1;
                end;
            end;
            
            timeLogger.block(block).trial(trial).end = GetSecs-tLstart;
            timeLogger.block(block).trial(trial).imageLength = timeLogger.block(block).trial(trial).imageEnd - timeLogger.block(block).trial(trial).start;
            timeLogger.block(block).trial(trial).blankLength = timeLogger.block(block).trial(trial).end - timeLogger.block(block).trial(trial).imageEnd;
            timeLogger.block(block).trial(trial).trialLength = timeLogger.block(block).trial(trial).end - timeLogger.block(block).trial(trial).start;
            WaitSecs(0.1);
            Eyelink('StopRecording');
            WaitSecs(0.001);
            Eyelink('message', 'TRIAL_RESULT 0');
            WaitSecs(.01);
            eyeLinkTrial = eyeLinkTrial + 1;
        end;
    end;
    
    SCTASK.ANSMAT = ANSMAT;
    timeLogger.block(block).blockLength = GetSecs - blockStart;
    timeLogger.block(block).blockEnd = realTimeExpected;
    fprintf('Block Time: %1.4f\n', timeLogger.block(block).blockLength);    
    save(PATH, 'SCTASK', 'timeLogger');
end;

%% Logging & Cleanup
% Create Covariate Files
% cnt = [1 1 1 1 1];
% 
% for i = 1:numBlocks
%     
%     switch TIME_MAT(i, 1) %condition no
%         case 1 % Faces
%             COVAR{1}(cnt(1), :) = TIME_MAT(i, 2:4);
%             cnt(1) = cnt(1) + 1;
%         case 2 % Objects
%             COVAR{2}(cnt(2), :) = TIME_MAT(i, 2:4);
%             cnt(2) = cnt(2) + 1;
%         case 3 % ObjectsScram
%             COVAR{3}(cnt(3), :) = TIME_MAT(i, 2:4);
%             cnt(3) = cnt(3) + 1;
%         case 4 % Places
%             COVAR{4}(cnt(4), :) = TIME_MAT(i, 2:4);
%             cnt(4) = cnt(4) + 1;
%         case 5 % Fixation
%             COVAR{5}(cnt(5), :) = TIME_MAT(i, 2:4);
%             cnt(5) = cnt(5) + 1;
%     end;
% end
% 
% for j = 1:5
%     dlmwrite(sprintf('FOSS_Sub0%d_Run%d_Cov%d.txt', sub, cbl, j), COVAR{j}, 'delimiter', '\t', 'precision', 4);
% end;

%% Shutdown Procedures
ShowCursor;
Screen('CloseAll');

Eyelink('Command', 'set_idle_mode');
WaitSecs(0.5);
Eyelink('CloseFile');

% download data file
try
    fprintf('Receiving data file ''%s''\n', edfFile );
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edfFile );
end

Eyelink('Shutdown');

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
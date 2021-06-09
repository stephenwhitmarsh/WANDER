% function SSEF_experiment_tactile_2
% Entrainment experiment
% 28-01-2016
% stephen.whitmarsh@ens.fr
% ShowCursor

addpath C:\MATLAB_toolboxes\ParallelPort
addpath C:\MATLAB_toolboxes\Lumina
addpath C:\MATLAB_toolboxes\Psychtoolbox
addpath C:\manips\STEPHEN

close all
clear all

HideCursor;
AssertOpenGL;

subID = input('Please type in the subject number ','s');
if isempty(subID)
    subID = '0';
end
practice_str = input('Practice (1) or for real ','s');

if strcmp(practice_str,'1')
    % For practice
    nrTrials            = 10;
    nrBlocks            = 1;
    TargetFraction      = 0.5;              % 0.1 = 10%
    setup.Eye = false;
else
    % For real
    nrTrials            = 50;
    nrBlocks            = 4;
    TargetFraction      = 0.1;              % 0.1 = 10%
    setup.Eye = true;   
end

%% EXPERIMENT PARAMETERS

StimFreq            = 16;
totalTime           = 1/StimFreq;
onTime              = .014;
offTime             = totalTime - onTime;
BreaksEveryNtrials  = 10;
TrialLengthMean     = 16;               % mean (actual mean in final distribution)
TrialLengthMax      = 30;               % cut off
TrialLengthMin      = 10;               % min
TargetDuration      = 3;                % number of skipped stimulations

scaleWidth          = 400;
scaleHeight         = 135;
scaleNr             = 8;
RatingResponse      = zeros(1,nrTrials);
CorrectResponse     = zeros(1,nrTrials);

FixationDuration    = 1;                % before tactile stimulation onset, in seconds
AlertDuration       = 1.5;              % before tactile stimulation onset, in seconds
AlertDurationJitter = 1;
EndFixationDuration = 1;                % time after stim before rating

fixationouter       = 50;               % radius in pixels
fixationinner       = 8;                % radius in pixels
fixationcolor1      = [200 200 200];    % before steadystate
fixationcolor2      = [200 200 200];    % during steadystate
fixationcolor3      = [200 200 200];    % after steadystate
backgroundColor     = 0;
textColor           = 255;
centerX             = 0;                % offset X of stimulus
centerY             = 100;              % offset X of stimulus
W                   = 1400;             % screen width
H                   = 1050;             % screen height

%% Setup screen
clear screen
whichScreen                 = 1;
oldResolution               = Screen('Resolution', whichScreen, W, H, 60);
[window1, rect]             = Screen('Openwindow', whichScreen, backgroundColor,[],[],2);
slack                       = Screen('GetFlipInterval', window1)/2;
[oldFontName,oldFontNumber] = Screen(window1,'TextFont','Helvetica');
oldTextSize                 = Screen('TextSize', window1, 24);

%% Create fixation screens
[AlertFixScreen,rect]  = Screen('OpenOffscreenWindow', window1, backgroundColor);
Screen('FillOval',  AlertFixScreen,  fixationcolor1,[W/2-fixationinner+centerX H/2-fixationinner+centerY W/2+fixationinner+centerX H/2+fixationinner+centerY]);

[BeforeFixScreen,rect] = Screen('OpenOffscreenWindow', window1, backgroundColor);
Screen('FrameOval', BeforeFixScreen, fixationcolor1,[W/2-fixationouter+centerX H/2-fixationouter+centerY W/2+fixationouter+centerX H/2+fixationouter+centerY], 4);
Screen('FillOval',  BeforeFixScreen, fixationcolor1,[W/2-fixationinner+centerX H/2-fixationinner+centerY W/2+fixationinner+centerX H/2+fixationinner+centerY]);

[DuringFixScreen,rect] = Screen('OpenOffscreenWindow',window1, backgroundColor, []);
Screen('FrameOval', DuringFixScreen, fixationcolor2,[W/2-fixationouter+centerX H/2-fixationouter+centerY W/2+fixationouter+centerX H/2+fixationouter+centerY], 4);
Screen('FillOval',  DuringFixScreen, fixationcolor2,[W/2-fixationinner+centerX H/2-fixationinner+centerY W/2+fixationinner+centerX H/2+fixationinner+centerY]);

[AfterFixScreen,rect]  = Screen('OpenOffscreenWindow', window1, backgroundColor, []);
Screen('FrameOval', AfterFixScreen,  fixationcolor3,[W/2-fixationouter+centerX H/2-fixationouter+centerY W/2+fixationouter+centerX H/2+fixationouter+centerY], 4);
Screen('FillOval',  AfterFixScreen,  fixationcolor3,[W/2-fixationinner+centerX H/2-fixationinner+centerY W/2+fixationinner+centerX H/2+fixationinner+centerY]);

%% Initialize rand
rand('state', sum(100*clock));

%% Initialize parallel port
OpenParPort;
WriteParPort(0);
OpenLumina;

%% Screen priority
% Priority(MaxPriority(window1));
% Priority(2);

for iblock = 1:nrBlocks
    
    FlipRating          = [zeros(1,nrTrials/2), ones(1,nrTrials/2)];
    FlipRating          = FlipRating(randperm(nrTrials));
    FlipDetection       = [zeros(1,nrTrials/2), ones(1,nrTrials/2)];
    FlipDetection       = FlipDetection(randperm(nrTrials));
    TargetTrial         = zeros(1,nrTrials);
    TargetTrial(1:ceil(nrTrials*TargetFraction)) = 1;
    TargetTrial         = TargetTrial(randperm(nrTrials));
    
    % Exponential TrialDuration
    TrialDuration = exprnd(TrialLengthMean-TrialLengthMin,1,nrTrials);          % draw from exp. dist.
    for i = find(TrialDuration>TrialDuration-TrialLengthMin)                    % replace values above maxISI-minISI
        while TrialDuration(i)>TrialLengthMax-TrialLengthMin
            TrialDuration(i) = exprnd(TrialLengthMean-TrialLengthMin);
        end
    end
    TrialDuration = TrialDuration+TrialLengthMin;                       % add minISI
    
    %% Set target onset between 1 and total-1 trialduration
    MinTargetTime = ones(size(TrialDuration))*TrialLengthMin+1;
    MaxTargetTime = TrialDuration-1;
    TargetTime    = round(MinTargetTime*StimFreq + (MaxTargetTime-MinTargetTime) .* StimFreq .* rand);
% 
%     for iTrial = 1:nrTrials
%         %     TargetTime(iTrial)  = randi([1*StimFreq, (TrialDuration(iTrial)-1)*StimFreq],1);
%         TargetTime(iTrial) = MinTargetTime(iTrial)*StimFreq + (MaxTargetTime(iTrial)-MinTargetTime(iTrial))*StimFreq.*rand
%     end
    
    if setup.Eye == true,
        
        %% setup the Eyelink initialization
        dummymode = 0; % set to 1 to run in dummymode (using mouse as pseudo-eyetracker)
        el = EyelinkInitDefaults(window1);
        
        % Initialization of the connection with the Eyelink Gazetracker. Exit program if this fails.
        if ~EyelinkInit(dummymode, 1)
            fprintf('Eyelink Init aborted.\n');
            cleanup(useTrigger);  % cleanup function
            return
        end
        
        [v vs]=Eyelink('GetTrackerVersion');
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        
        % make sure that we get event data from the Eyelink
        Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS'); %"file_sample_data" specifies what type of samples will be wrtten to the EDF file
        Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS'); %"link_sample_data" specifies what type of samples will be sent to link
        Eyelink('command', 'drift_correct_cr_disable = OFF'); % To enable the drift correction procedure to adjust the calibration rather than simply allowing a drift check
        
        % open edf file for recording data from Eyelink
        EDFname = sprintf('S%s_B%s', num2str(subID), num2str(iblock));
        edfFile = [EDFname '.edf'];
        Eyelink('Openfile', edfFile);
        
        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(el);
        
        % start recording eye position
        Eyelink('StartRecording');
        
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        
        % mark zero-plot time in data file
        Eyelink('message', 'start recording Eyelink');
    end
    
    % Start screen experimentor
    Screen('FillRect', window1, backgroundColor)
    DrawFormattedText(window1, ['Experimenter: Press the space bar to begin block ' num2str(iblock) ' of ' num2str(nrBlocks)], 'center', 'center', textColor);
    Screen('Flip',window1);
    HideCursor(window1);
    WaitSecs(0.01);
    KbStrokeWait;
    
    % Start screen subject
    Screen('FillRect', window1, backgroundColor)
    DrawFormattedText(window1, 'Soyez prêt!', 'center', 'center', textColor);
    Screen('Flip',window1);
    WaitSecs(3);
    
    %% TRIAL LOOP
    itrial = 1;
    
    while (itrial <= nrTrials)
        disp(itrial);
        WaitSecs(0.100);
        
        % Start with alerting
        Screen('CopyWindow', AlertFixScreen, window1)
        tAlert = Screen('Flip', window1);
        WriteParPort(1);
        WaitSecs(0.003);
        WriteParPort(0);
        
        %% ADDED 30-5-2016
        if itrial == 1 & ~strcmp(practice_str,'1')
            AlertDuration = 60;
        else
            AlertDuration = 1.5;
        end
        %%
        
        if setup.Eye == true,
            Eyelink ('Message', [num2str(itrial) 'ALERT' ]);
        end
        
        % Start with before fixation
        Screen('CopyWindow', BeforeFixScreen, window1)
        tFixation = Screen('Flip', window1, tAlert + AlertDuration + rand*AlertDurationJitter - slack);
        WriteParPort(2);
        WaitSecs(0.003);
        WriteParPort(0);
        
        if setup.Eye == true,
            Eyelink ('Message', [num2str(itrial) 'FIX' ]);
        end
        
        % Switch to during fixation
        Screen('CopyWindow', DuringFixScreen, window1)
        flipTime = Screen('Flip', window1, tFixation + FixationDuration - slack);
        HideCursor(window1);
        
        % Start stimulus
        for iTrigger = 1 : TrialDuration(itrial)*StimFreq
            t0 = GetSecs;
            if ~TargetTrial(itrial) 
                % NON-TARGET
                if iTrigger == 1
                    WriteParPort(129);
                else
                    WriteParPort(128);
                end
                WaitSecs('UntilTime', t0 + onTime);
                WriteParPort(0);
                pause(0.01);
                WaitSecs('UntilTime', t0 + onTime + offTime);
            else
                % TARGET
                if iTrigger < TargetTime(itrial) | iTrigger >= TargetTime(itrial)+TargetDuration
                    if iTrigger == 1
                        WriteParPort(131);
                    else
                        WriteParPort(130);
                    end
                    WaitSecs('UntilTime', t0 + onTime);
                    WriteParPort(0);
                    pause(0.01);
                    WaitSecs('UntilTime', t0 + onTime + offTime);
                else
                    WriteParPort(127);
                    WaitSecs('UntilTime', t0 + onTime);
                    WriteParPort(0);
                    pause(0.01);
                    WaitSecs('UntilTime', t0 + onTime + offTime);
                end
            end
        end
        
        % Switch to after fixation
        Screen('CopyWindow', AfterFixScreen, window1)
        tTrialOffset = Screen('Flip', window1);
        WriteParPort(3);
        WaitSecs(0.003);
        WriteParPort(0);
        
        if setup.Eye == true
            Eyelink ('Message', [num2str(itrial) 'ENDSTIM' ]);
        end
        
        % Draw attention rating screen after EndFixation duration
        Screen('FillRect', window1, backgroundColor)
        rating = randi(scaleNr);
        if FlipRating(itrial)
            for i = 1 : scaleNr
                colorArray(:,i) = [150, 150, 150];
                rectArray(:,i) = [W/2-scaleWidth/2 + (scaleWidth/scaleNr)*(i-1),  H/2+centerY-scaleHeight/scaleNr*i,           W/2-scaleWidth/2 + (scaleWidth/scaleNr)*i, H/2+centerY];
            end
        else
            for i = 1 : scaleNr
                colorArray(:,i) = [150, 150, 150];
                rectArray(:,i) = [W/2-scaleWidth/2 + (scaleWidth/scaleNr)*(i-1),  H/2+centerY-scaleHeight/scaleNr*(scaleNr-i+1), W/2-scaleWidth/2 + (scaleWidth/scaleNr)*i, H/2+centerY];
            end
        end
        
        Screen('FillRect',window1, colorArray, rectArray);
        Screen('FillRect',window1, [255 255 255], rectArray(:,rating), 6);
        AttentionRatingTime = Screen('Flip',window1, tTrialOffset + EndFixationDuration - slack);
        
        WriteParPort(4);
        WaitSecs(0.003);
        WriteParPort(0);
        
        if setup.Eye == true,
            Eyelink ('Message', [num2str(itrial) 'RATING' ]);
        end
        
        % Get attention rating response
        exitResponse = false;
        moved = false;
        while exitResponse == false
            
            % wait for buttonpress, left: 49, right:50, validate: 52
            Buttons = WaitForLumina;
            
            if Buttons == 52 & moved == true
                exitResponse = true;
            elseif Buttons == 49
                moved = true;
                if rating > 1
                    rating = rating - 1;
                    Screen('FillRect',window1, backgroundColor)
                    Screen('FillRect',window1, colorArray, rectArray);
                    Screen('FillRect',window1, [255 255 255], rectArray(:,rating));
                    Screen('Flip',window1);
                end
            elseif Buttons == 50
                moved = true;
                if rating < scaleNr
                    rating = rating + 1;
                    Screen('FillRect',window1, backgroundColor)
                    Screen('FillRect',window1, colorArray, rectArray);
                    Screen('FillRect',window1, [255 255 255], rectArray(:,rating));
                    Screen('Flip',window1);
                end
            end
        end
        
        if FlipRating(itrial)
            rating = scaleNr-rating+1;
        end
        
        RatingResponse(itrial) = rating;
        SafeWaitSecs(0.1); % to prevent trigger timing conflict
        WriteParPort(10+rating);
        SafeWaitSecs(0.003);
        WriteParPort(0);
        
        % Draw target detection screen
        textY           = -100;
        textX           = -160;
        sizeBox = 50;
        reponse = randi(2);
        clear colorArray rectArray
        rectArray  = [W/2-100-sizeBox, H/2+centerY-sizeBox, W/2-100+sizeBox, H/2+centerY+sizeBox; W/2+100-sizeBox, H/2+centerY-sizeBox, W/2+100+sizeBox, H/2+centerY+sizeBox]';
        Screen('FillRect', window1, backgroundColor)
        Screen('FrameRect', window1, [150, 150, 150], rectArray+[-5 -5 5 5;-5 -5 5 5]',6);
        if FlipDetection(itrial)
            DrawFormattedText(window1, 'TARGET', W/2+35, H/2+centerY+textY, [255 255 255]);
        else
            DrawFormattedText(window1, 'TARGET', W/2-167, H/2+centerY+textY, [255 255 255]);
        end
        Screen('Flip',window1);
        
        % Get discrimination response
        exitResponse = false;
        moved = false;
        while exitResponse == false
            
            % wait for buttonpress, left: 49, right:50, validate: 52
            Buttons = WaitForLumina;
            
            if Buttons == 52 & moved == true
                exitResponse = true;
            elseif Buttons == 49
                moved = true;
                response = 1;
                Screen('FillRect', window1, backgroundColor)
                Screen('FillRect', window1, [255, 255, 255], rectArray(:,response))
                Screen('FrameRect', window1, [150, 150, 150], rectArray+[-5 -5 5 5;-5 -5 5 5]',6);
                if FlipDetection(itrial)
                    DrawFormattedText(window1, 'TARGET', W/2+35, H/2+centerY+textY, [255 255 255]);
                else
                    DrawFormattedText(window1, 'TARGET', W/2-167, H/2+centerY+textY, [255 255 255]);
                end
                Screen('Flip',window1);
            elseif Buttons == 50
                moved = true;
                response = 2;
                Screen('FillRect', window1, backgroundColor)
                Screen('FillRect', window1, [255, 255, 255], rectArray(:,response))
                Screen('FrameRect', window1, [150, 150, 150], rectArray+[-5 -5 5 5;-5 -5 5 5]',6);
                if FlipDetection(itrial)
                    DrawFormattedText(window1, 'TARGET', W/2+35, H/2+centerY+textY, [255 255 255]);
                else
                    DrawFormattedText(window1, 'TARGET', W/2-167, H/2+centerY+textY, [255 255 255]);
                end
                Screen('Flip',window1);
            end
        end
        if FlipDetection(itrial)
            response = 2-response+1;
        end
        
        if TargetTrial(itrial) == 1 & response == 1
            CorrectResponse(itrial) = 1; % hit
        elseif TargetTrial(itrial) == 1 & response == 2
            CorrectResponse(itrial) = 2; % miss
        elseif TargetTrial(itrial) == 0 & response == 1
            CorrectResponse(itrial) = 3; % false alarm
        elseif TargetTrial(itrial) == 0 & response == 2
            CorrectResponse(itrial) = 4; % correct rejection
        end
        
        SafeWaitSecs(0.01); % to prevent trigger timing conflict
        WriteParPort(20+CorrectResponse(itrial));
        SafeWaitSecs(0.003);
        WriteParPort(0);
        
        % get a break every N trials
        if mod(itrial,BreaksEveryNtrials) == 0 && itrial < nrTrials
            
            if ~isempty(find(CorrectResponse(itrial-BreaksEveryNtrials+1:itrial)==2))
                misses = size(find(CorrectResponse(itrial-BreaksEveryNtrials+1:itrial)==2),2);
            else
                misses = 0;
            end
            if ~isempty(find(CorrectResponse(itrial-BreaksEveryNtrials+1:itrial)==3))
                FA     = size(find(CorrectResponse(itrial-BreaksEveryNtrials+1:itrial)==3),2);
            else
                FA = 0;
            end
            
            Screen('FillRect', window1, backgroundColor)
            DrawFormattedText(window1, ['You missed ', num2str(misses), ' target(s) \n \n You gave ', num2str(FA), ' false alarm(s) \n \n Take a short break (but do not move) and press any button to continue.'], 'center', H/2 + centerY, textColor);
            Screen('Flip',window1)
            
            WriteParPort(5);
            WaitSecs(0.003);
            WriteParPort(0);
            
            if setup.Eye == true,
                Eyelink ('Message', [num2str(itrial) 'PAUSE' ]);
            end

            SafeWaitSecs(1);
            
            % wait for buttonpress
            WaitForLumina
        end % trial
        
        itrial = itrial + 1;
    end  % block
    
    if setup.Eye == true,
        Screen('FillRect', window1, backgroundColor)
        EndOfBlockTime = DrawFormattedText(window1, 'End of this block. Please do not move yet', 'center', H/2 + centerY, textColor);
        Screen('Flip',window1)
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        fprintf('Receiving data file ''%s''\n', edfFile );
        setup.eyefilename = sprintf('eye_S%s_block%s_%s.edf', num2str(subID), num2str(iblock), datestr(now,'ddmmyy-HHMM-SSFFF'));
        cd('Eyelink');
        status=Eyelink('ReceiveFile'); %this collects the file from the eyelink
        movefile(edfFile, setup.eyefilename);
        % go one folder upwards
        cd ..   
    end
    
    if ~strcmp(practice_str,'1')
        Screen('FillRect', window1, backgroundColor)
        DrawFormattedText(window1, 'OK', 'center', H/2 + centerY, textColor);
        Screen('Flip',window1,EndOfBlockTime+60);
    end
    %%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RestrictKeysForKbCheck([]);
Screen(window1,'Close');
close all
sca;
return
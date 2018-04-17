function foragingChoiceworldNew(t, events, parameters, visStim, inputs, outputs, audio)
% basicChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% 2017-03 - AP created
% 2018-01 - MW updated: automatic reward reduction, L-R performance
% 2018-02 - MM updated: flipping to 2AFC blending
% 2018-04 - ED updated: change target to world coords, seperate de-bias from task mix
%
% Choice world that adapts with behavior
%
% Task structure: 
% Start trial
% Resetting pre-stim quiescent period
% Stimulus onset
% Fixed cue interactive delay
% Infinite time for response, fix stim azimuth on response
% Short ITI on reward, long ITI on punish, then turn stim off
% End trial


    %% --- Fixed parameters ---
    % move truely fixed paramters here (e.g. not updated when the task is
    % deployed during learning or production
    
    % hardware related
    % audioSampleRate = 192000;
    audioSampleRate = 44100;
    encoderRes = 1024; % Resolution of the rotary encoder
    % TODO: why is this a 'signal' and not conventional math? is this not a
    % static param?#
    millimetersFactor = events.newTrial.map2(31*2*pi/(encoderRes*4), @times); % convert the wheel gain to a value in mm/deg

    %% --- Session level parameters ---
    % TODO: should target be pased in to initializeTrialData? or is it only modified params?
    targetAzimuthLeft = 90;
    targetAzimuthRight = 90;

    
    %% --- Trial level parameters ---
    %% Trial meta parameters
    
    % TODO: these all need to be set in defaultSessionParams

    %sigmoid parameter from (performanceTracer-targetPerformance) -> deltaFlipRate
    sigmoidParamsFlip = [0.7, 0.25, 0.01, -0.01]; %targetPerformance, slope, max, min
    performanceTracker = sigmoidParamsFlip(1); %variable to track performance
    alphaPerformanceTracker = 0.05; %weight on current trial to update performance
    % deltaswitchProb = 0.0005; %how much to increase 

    % sigmoidParamsSideHistory = [0.5, 0.25, 1, 0];
    sigmoidParamsSideHistory = [0.5, 0.1, 1, 0]; %Ines changed
    SideHistoryTracker = sigmoidParamsSideHistory(1);
    alphaSideHistoryTracker = 0.05; 
    %slopeSigmoidSideHistory = 0.25; %min/max is 0/1 for SideHistory sigmoid, midpoint = 0.5

    % sigmoidParamsBias = [0.5, 0.25, 1, 0];
    sigmoidParamsBias = [0.5, 0.1, 1, 0]; %Ines changed
    BiasTracker = sigmoidParamsBias(1);
    alphaBiasTracker = 0.05; 
    %slopeSigmoidBias = 0.25; %min/max is 0/1 for Bias sigmoid, midpoint = 0.5


    %% Trial choice parameters
    rewardProb = 1;
    % flipRate = 0.1; %Ines changed
    flipRate = 0.25;
    % Starting reward size
    rewardSize = 3;


    %% Stimulus/target
    % (which contrasts to use)
%     contrasts = [1,0.5,0.25,0.125,0.06,0];
    contrasts = [1,0.5];
    sigma = [20,20];
    spatialFreq = 1/15;
    startingAzimuthLeft = 90;
    startingAzimuthRight = 90;
    responseCenter = 0;
    responseWidth = 0;

    %% Timing
    % prestimQuiescentTime = 0.5;
    prestimQuiescentTime = 0.2;
    % cueInteractiveDelay = 0.5;
    cueInteractiveDelay = 0.1;
    itiHit = 1;
    itiMiss = 2;

    %% Sounds

    onsetToneAmplitude = 0.2;
    onsetToneFreq = 6000;
    onsetToneDuration = 0.1;
    onsetToneRampDuration = 0.01;
    audioChannels = 2;
    
    % TODO: why is this a 'signal' and not conventional math? is this not a static param?
    toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
        aud.pureTone(onsetToneFreq,onsetToneDuration,audioSampleRate, ...
        onsetToneRampDuration,audioChannels));

    missNoiseDuration = 0.5;
    missNoiseAmplitude = 0.02;
    % TODO: why is this a 'signal' and not conventional math? is this not a static param?
    missNoiseSamples = missNoiseAmplitude*events.expStart.map(@(x) ...
        randn(2, audioSampleRate*missNoiseDuration));

    %% Wheel parameters
    quiescThreshold = 1;
    % TODO: this should only be used in the very first trials to encourage
    % wheel interaction
    wheelGain = 5;

    %% --- Initialize trial data ---
    % 
    trialData = events.expStart.mapn( ...
        contrasts, ...
        sigmoidParamsFlip, performanceTracker, alphaPerformanceTracker, ...
        sigmoidParamsSideHistory, SideHistoryTracker, alphaSideHistoryTracker, ...
        sigmoidParamsBias, BiasTracker, alphaBiasTracker, ...
        rewardProb, flipRate, rewardSize, wheelGain, startingAzimuthLeft, startingAzimuthRight, ...
        @initializeTrialData).subscriptable;

    % trialDataInit = events.expStart.mapn( ...
    %     contrasts, rewardProb, switchProb, rewardSize, wheelGain,...
    %     @initializeTrialData).subscriptable;
    
    

    %% Set up wheel 
    wheel = inputs.wheel.skipRepeats();


    % startingAzimuth = 90;

    %% Trial event times
    % (this is set up to be independent of trial conditon, that way the trial
    % condition can be chosen in a performance-dependent manner)

    % Resetting pre-stim quiescent period
    prestimQuiescentPeriod = at(prestimQuiescentTime,events.newTrial.delay(0)); 
    preStimQuiescence = sig.quiescenceWatch(prestimQuiescentPeriod, t, wheel, quiescThreshold); 

    % Stimulus onset
    %stimOn = sig.quiescenceWatch(preStimQuiescPeriod, t, wheel, quiescThreshold); 
    stimOn = at(true,preStimQuiescence); 
%     stimOn = false;
    
    % Fixed cue interactive delay
    interactiveOn = stimOn.delay(cueInteractiveDelay); 

    % Play tone at interactive onset
    audio.onsetTone = toneSamples.at(interactiveOn);

    % Response
    % (wheel displacement zeroed at interactiveOn)
    % TODO: verify replacement of old azimuth cond
    % TODO: fix so this loads for all trial, then adjust
    wheelDisplacement = cond(...
        interactiveOn, wheelGain*millimetersFactor*(wheel - wheel.at(interactiveOn)), ...
        true, 0);
        % wheel is current angle signal, subtracting from interactive on
    
    % stimulus is located at 
    stimLocation = cond(interactiveOn, wheelDisplacement + trialData.startingAzimuth, ...
                        true, trialData.startingAzimuth); 

    % TODO: define new displacement algorithm
    stimInResponseWindow = interactiveOn.setTrigger( ...
        stimLocation >= responseCenter - responseWidth |  ...
        stimLocation <= responseCenter + responseWidth);
    
    % stimDisplacement
    
    
    
    % threshold = interactiveOn.setTrigger(abs(stimDisplacement) ...
    %     >= responseDisplacement);
%     response = at(-sign(stimDisplacement), threshold);
    response = stimInResponseWindow;
    
    %% Update performance at response

    responseData = stimLocation;
    % Update performance
    trialData = responseData.at(response).scan(@updateTrialData,trialData).subscriptable;
        % maps @updateTrialData on responseData when the response occurs; inits with trialData, then uses t-1
    trialContrast = trialData.trialContrast;

    %% Give feedback and end trial

    % Give reward on hit
    % NOTE: there is a 10ms delay for water output, because otherwise water and
    % stim output compete and stim is delayed
    % TODO: understand why there is this competition, how it's handled and the tolerance
    % TODO: why is this trialDataInit and not trialData
    water = at(trialData.rewardSize,trialData.hit.delay(0.01));  
    outputs.reward = water;

    % Play noise on miss
    audio.missNoise = missNoiseSamples.at(trialData.miss.delay(0.01));

    % ITI defined by outcome
    iti = iff(trialData.hit==1, itiHit, itiMiss);
    % iti = iff(eq(trialData.miss, true), abs(response)*itiHit, abs(response)*itiMiss).at(response);

    % Stim stays on until the end of the ITI
    stimOff = stimInResponseWindow.delay(iti);

    %% Visual stimulus

    % Azimuth control
    % 1) stim fixed in place until interactive on [implicit now]
    % 2) wheel-conditional during interactive
    % 3) fixed at target center after response [implicit now]
    azimuth = cond( ...
        events.newTrial.to(events.newTrial), stimLocation);

    stim = vis.grating(t, 'square', 'gaussian');
    stim.sigma = sigma;
    stim.spatialFreq = spatialFreq;
    stim.phase = 2*pi*events.newTrial.map(@(v)rand);
    stim.azimuth = azimuth;
        % this is the dynamic azimoth signal
    stim.contrast = trialContrast.at(events.newTrial);
    stim.show = stimOn.to(stimOff);
        % timing

    visStim.stim = stim;

    %% --- events are time stamped and logged ---
    % Display and save

    % Wheel and stim
    events.azimuth = azimuth;

    % Trial times
    events.stimOn = stimOn;
    events.stimOff = stimOff;
    events.interactiveOn = interactiveOn;
    events.response = response;
    % events.endTrial = at(~trialData.repeatTrial, stimOff);
    events.endTrial = stimOff;

    % Performance
    events.contrasts = trialData.contrasts;
    events.trialContrast = trialData.trialContrast;
    events.trialSide = trialData.trialSide;
    events.repeatTrial = trialData.repeatTrial;
    events.hit = trialData.hit.at(response);
    events.totalWater = water.scan(@plus, 0).map(fun.partial(@sprintf, '%.1fµl'));

    events.flipRate = trialData.flipRate;
    events.performanceTracker = trialData.performanceTracker;
    events.SideHistoryTracker = trialData.SideHistoryTracker;
    events.BiasTracker = trialData.BiasTracker;

end

function trialDataInit = initializeTrialData(expRef, ...
    contrasts, ...
    sigmoidParamsFlip, performanceTracker,alphaPerformanceTracker, ...
    sigmoidParamsSideHistory, SideHistoryTracker, alphaSideHistoryTracker, ...
    sigmoidParamsBias, BiasTracker, alphaBiasTracker, ...
    rewardProb, flipRate, rewardSize, wheelGain, startingAzimuthLeft, startingAzimuthRight)

    % function trialDataInit = initializeTrialData(expRef, ...
    %     contrasts, rewardProb, switchProb, rewardSize, wheelGain)

    %%%% Get the subject
    % (from events.expStart - derive subject from expRef)
    subject = dat.parseExpRef(expRef);

    %%%% Initialize all of the session-independent performance values
    trialDataInit = struct;

    trialDataInit.contrasts = contrasts;
    trialDataInit.trialContrast = randsample(contrasts,1);
    trialDataInit.sigmoidParamsFlip = sigmoidParamsFlip;
    trialDataInit.performanceTracker = performanceTracker; %don't need it here?
    trialDataInit.alphaPerformanceTracker = alphaPerformanceTracker;
    trialDataInit.sigmoidParamsSideHistory = sigmoidParamsSideHistory;
    trialDataInit.SideHistoryTracker = SideHistoryTracker;
    trialDataInit.alphaSideHistoryTracker = alphaSideHistoryTracker;
    trialDataInit.sigmoidParamsBias = sigmoidParamsBias;
    trialDataInit.BiasTracker = BiasTracker;
    trialDataInit.alphaBiasTracker = alphaBiasTracker;
    trialDataInit.flipRate = flipRate; %don't need it here?
    trialDataInit.rewardProb = rewardProb;
    trialDataInit.rewardSize = rewardSize; %don't need it here?
    trialDataInit.wheelGain = wheelGain; %don't need it here?
    trialDataInit.startingAzimuthLeft = startingAzimuthLeft;
    trialDataInit.startingAzimuthRight = startingAzimuthRight;

    % Set the first trial side randomly
    trialDataInit.trialSide = randsample([-1,1],1);
    % Setup starting azimuth
    % TODO: should this be a 'cond' outside? or here? or inline ifelse?
    if trialDataInit.trialSide < 0
        trialDataInit.startingAzimuth = trialDataInit.startingAzimuthLeft;
    else
        trialDataInit.startingAzimuth = trialDataInit.startingAzimuthRight;
    end
    % Set up the flag for repeating incorrect
    trialDataInit.repeatTrial = false;
    % Initialize hit/miss
    trialDataInit.hit = false;
    trialDataInit.miss = false;

    %%%% Load the last experiment for the subject if it exists
    % (note: MC creates folder on initilization, so start search at 1-back)
    expRef = dat.listExps(subject);
    useOldParams = false;
    if length(expRef) > 1
        % Loop through blocks from latest to oldest, if any have the relevant
        % parameters then carry them over
        for check_expt = length(expRef)-1:-1:1
            previousBlockFilename = dat.expFilePath(expRef{check_expt}, 'block', 'master');
            if exist(previousBlockFilename,'file')
                previousBlock = load(previousBlockFilename);
                if ~isfield(previousBlock.block, 'outputs')||isempty(previousBlock.block.outputs.rewardValues)
                    lastRewardSize = rewardSize;
                else
                    lastRewardSize = previousBlock.block.outputs.rewardValues(end);
                end
        
    %             previousBlock.block
        
                if isfield(previousBlock.block,'events')
                    previousBlock = previousBlock.block.events;
                else
                    previousBlock = [];
                end
        
    %             previousBlock
        
                if isempty(previousBlock) || ~isfield(previousBlock, 'flipRateValues')
                    lastFlipRate = flipRate;
                else
                    lastFlipRate = previousBlock.flipRateValues(end);
    % lastFlipRate = flipRate;
                end
        
                if isempty(previousBlock) || ~isfield(previousBlock, 'performanceTrackerValues')
                    lastPerformanceTracker = performanceTracker;
                else
                    lastPerformanceTracker = previousBlock.performanceTrackerValues(end);
                end
        
                if isempty(previousBlock) || ~isfield(previousBlock, 'SideHistoryTrackerValues')
                    lastSideHistoryTracker = SideHistoryTracker;
                else
                    lastSideHistoryTracker = previousBlock.SideHistoryTrackerValues(end);
                end

                if isempty(previousBlock) || ~isfield(previousBlock, 'BiasTrackerValues')
                    lastBiasTracker = BiasTracker;
                else
                    lastBiasTracker = previousBlock.BiasTrackerValues(end);
                end
            end
    
            % Check if the relevant fields exist
            if exist('previousBlock','var') && ...
                    length(previousBlock.newTrialValues) > 5 
                % Break the loop and use these parameters
                useOldParams = true;
                break
            end       
        end        
    end

    if useOldParams
        % If the last experiment file has the relevant fields, set up performance

        % If the subject did over 200 trials last session, reduce the reward by
        % 0.1, unless it is 2ml
        if length(previousBlock.newTrialValues) > 200 && lastRewardSize > 2
            trialDataInit.rewardSize = lastRewardSize-0.1;
        else
            trialDataInit.rewardSize = lastRewardSize;
        end

        trialDataInit.flipRate = lastFlipRate;
        trialDataInit.performanceTracker = lastPerformanceTracker;
        trialDataInit.SideHistoryTracker = lastSideHistoryTracker;
        trialDataInit.BiasTracker = lastBiasTracker;

    else
        % If this animal has no previous experiments, initialize performance 
        % Initialize water reward size & wheel gain
        trialDataInit.flipRate = flipRate;
        trialDataInit.performanceTracker = performanceTracker;
        trialDataInit.SideHistoryTracker = SideHistoryTracker;
        trialDataInit.BiasTracker = BiasTracker;
        trialDataInit.rewardSize = rewardSize;
        trialDataInit.wheelGain = wheelGain;
    end
end

function trialData = updateTrialData(trialData,responseData)
% Update the performance 

    stimDisplacement = responseData;

    % Next contrast
    trialData.trialContrast = randsample(trialData.contrasts,1); 

    %%%% Define response type based on trial condition
    r = rand;
    trialData.hit = stimDisplacement*trialData.trialSide < 0 && r <= trialData.rewardProb; %correct side and reward
    trialData.miss = stimDisplacement*trialData.trialSide > 0 || ... %wrong side
        (stimDisplacement*trialData.trialSide < 0 && r > trialData.rewardProb); %correct side and no reward

    %update flipRate
    trialData.performanceTracker = trialData.alphaPerformanceTracker*trialData.hit ...
        + (1-trialData.alphaPerformanceTracker)*trialData.performanceTracker;
    trialData.sigmoidParamsFlip
    % deltaFlipRate = trialData.sigmoidParamsFlip(4) ...
    %     + (trialData.sigmoidParamsFlip(3)-trialData.sigmoidParamsFlip(4)) ...
    %     /(1+exp(-(trialData.performanceTracker - trialData.sigmoidParamsFlip(1))/trialData.sigmoidParamsFlip(2)))
    % trialData.flipRate = max(min(trialData.flipRate + deltaFlipRate,1),0.1);
    % %Ines changed
    deltaFlipRate = trialData.sigmoidParamsFlip(4) ...
        + (trialData.sigmoidParamsFlip(3)-trialData.sigmoidParamsFlip(4)) ...
        /(1+exp(-(trialData.performanceTracker - trialData.sigmoidParamsFlip(1))/trialData.sigmoidParamsFlip(2)))
    trialData.flipRate = max(min(trialData.flipRate + deltaFlipRate,1),0.25);
    % trialData.flipRate = 1;

    %update SideHistory
    %SideHistoryTracker: 0,left"ward" SideHistory; 1,right"ward" SideHistory, 
    trialData.SideHistoryTracker = trialData.alphaSideHistoryTracker*(stimDisplacement>0) ...
        + (1-trialData.alphaSideHistoryTracker)*trialData.SideHistoryTracker;
    trialData.sigmoidParamsSideHistory
    sigmoidOutSideHistory = trialData.sigmoidParamsSideHistory(4) ...
        + (trialData.sigmoidParamsSideHistory(3)-trialData.sigmoidParamsSideHistory(4))...
        /(1+exp(-(trialData.SideHistoryTracker-trialData.sigmoidParamsSideHistory(1))/trialData.sigmoidParamsSideHistory(2)))

    %update Bias
    %BiasTracker: 0,left"ward" Bias; 1,right"ward" Bias, 
    trialData.BiasTracker = trialData.alphaBiasTracker*(stimDisplacement>0) ... % TODO: fix for correct/error
        + (1-trialData.alphaBiasTracker)*trialData.BiasTracker;
    trialData.sigmoidParamsBias
    sigmoidOutBias = trialData.sigmoidParamsBias(4) ...
        + (trialData.sigmoidParamsBias(3)-trialData.sigmoidParamsBias(4))...
        /(1+exp(-(trialData.BiasTracker-trialData.sigmoidParamsBias(1))/trialData.sigmoidParamsBias(2)))

    %%%% Set flag to repeat - skip trial choice if so
    % if trialData.miss && ...
    %         ismember(trialData.trialContrast,trialData.contrasts(trialData.repeatOnMiss))

    trialData.repeatTrial = false;

    %%%% Pick next side (this is done at random)
    randomFlip = rand
    randomSideHistory = rand
    if randomFlip < trialData.flipRate
        %now can possibly change side
        nextSide = ((randomSideHistory - sigmoidOutSideHistory)<0)*2-1 %if rightward SideHistory, sigmoidOutSideHistory~=1, nextSide<-1 -> L
        trialData.trialSide = nextSide
    end

end
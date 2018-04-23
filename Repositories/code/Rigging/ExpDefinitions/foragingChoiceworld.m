function foragingChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% basicChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% 2017-03 - AP created
% 2018-01 - MW updated: automatic reward reduction, L-R performance
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


	%% Fixed parameters

	% Trial meta parameters

	%sigmoid parameter from (performanceTracer-targetPerformance) -> deltaFlipRate
	sigmoidParamsFlip = [0.7, 0.25, 0.01, -0.01]; %targetPerformance, slope, max, min
	% targetPerformance = 0.7; %if performance is above this, increase switch prob and vice versa
	% maxDeltaFlipRate = 0.01;
	% minDeltaFlipRate = -maxDeltaFlipRate; %symmetric for now
	% slopeSigmoidFlipRate = 0.25;
	performanceTracker = sigmoidParamsFlip(1); %variable to track performance
	alphaPerformanceTracker = 0.05; %weight on current trial to update performance
	% deltaswitchProb = 0.0005; %how much to increase 


	% sigmoidParamsBias = [0.5, 0.25, 1, 0];
	sigmoidParamsBias = [0.5, 0.1, 1, 0]; %Ines changed
	biasTracker = sigmoidParamsBias(1);
	alphaBiasTracker = 0.05; 
	%slopeSigmoidBias = 0.25; %min/max is 0/1 for bias sigmoid, midpoint = 0.5

	% Trial choice parameters
	rewardProb = 0.7;
	% flipRate = 0.1; %Ines changed
	flipRate = 0.25;


	% Stimulus/target
	% (which contrasts to use)
	contrasts = [1,0.5,0.25,0.125,0.06,0];
	% contrasts = [1,0.5];
	sigma = [20,20];
	spatialFreq = 1/15;
	startingAzimuth = 90;
	responseDisplacement = 90;
	% Starting reward size
	rewardSize = 3;

	% Timing
	% prestimQuiescentTime = 0.5;
	prestimQuiescentTime = 0.2;
	% cueInteractiveDelay = 0.5;
	cueInteractiveDelay = 0.1;
	itiHit = 1;
	itiMiss = 2;

	% Sounds
	% audioSampleRate = 192000;
	audioSampleRate = 44100;

	onsetToneAmplitude = 0.2;
	onsetToneFreq = 6000;
	onsetToneDuration = 0.1;
	onsetToneRampDuration = 0.01;
	audioChannels = 2;
	toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
	    aud.pureTone(onsetToneFreq,onsetToneDuration,audioSampleRate, ...
	    onsetToneRampDuration,audioChannels));

	missNoiseDuration = 0.5;
	missNoiseAmplitude = 0.02;
	missNoiseSamples = missNoiseAmplitude*events.expStart.map(@(x) ...
	    randn(2, audioSampleRate*missNoiseDuration));

	% Wheel parameters
	quiescThreshold = 1;
	encoderRes = 1024; % Resolution of the rotary encoder
	millimetersFactor = events.newTrial.map2(31*2*pi/(encoderRes*4), @times); % convert the wheel gain to a value in mm/deg
	wheelGain = 15;

	%% Initialize trial data

	trialDataInit = events.expStart.mapn( ...
	    contrasts, sigmoidParamsFlip, performanceTracker,alphaPerformanceTracker, ...
	    sigmoidParamsBias, biasTracker, alphaBiasTracker, ...
	    rewardProb, flipRate, rewardSize, wheelGain,...
	    @initializeTrialData).subscriptable;

	% trialData = trialDataInit.mapn( ...
	%     NaN, ...
	%     @updateTrialData).subscriptable;

	% trialDataInit = events.expStart.mapn( ...
	%     contrasts, rewardProb, switchProb, rewardSize, wheelGain,...
	%     @initializeTrialData).subscriptable;
	%% Set up wheel 

	wheel = inputs.wheel.skipRepeats();

	%% Trial event times
	% (this is set up to be independent of trial conditon, that way the trial
	% condition can be chosen in a performance-dependent manner)

	% Resetting pre-stim quiescent period
	prestimQuiescentPeriod = at(prestimQuiescentTime,events.newTrial.delay(0)); 
	preStimQuiescence = sig.quiescenceWatch(prestimQuiescentPeriod, t, wheel, quiescThreshold); 

	% Stimulus onset
	%stimOn = sig.quiescenceWatch(preStimQuiescPeriod, t, wheel, quiescThreshold); 
	stimOn = at(true,preStimQuiescence); 

	% Fixed cue interactive delay
	interactiveOn = stimOn.delay(cueInteractiveDelay); 

	% Play tone at interactive onset
	audio.onsetTone = toneSamples.at(interactiveOn);

	% Response
	% (wheel displacement zeroed at interactiveOn)
	stimDisplacement = wheelGain*millimetersFactor*(wheel - wheel.at(interactiveOn));

	threshold = interactiveOn.setTrigger(abs(stimDisplacement) ...
	    >= responseDisplacement);
	response = at(-sign(stimDisplacement), threshold);

	%% Update performance at response

	responseData = stimDisplacement;
	% Update performance
	trialData = responseData.at(response).scan(@updateTrialData,trialDataInit).subscriptable;
	trialContrast = trialData.trialContrast;

	%% Give feedback and end trial

	% Give reward on hit
	% NOTE: there is a 10ms delay for water output, because otherwise water and
	% stim output compete and stim is delayed
	% water = at(trialDataInit.rewardSize,trialData.hit.delay(0.01));  
	water = at(trialDataInit.rewardSize,trialData.isRewarded.delay(0.01));  
	outputs.reward = water;

	% Play noise on miss
	audio.missNoise = missNoiseSamples.at(trialData.miss.delay(0.01));

	% ITI defined by outcome
	iti = iff(trialData.hit==1, itiHit, itiMiss);
	% iti = iff(eq(trialData.miss, true), abs(response)*itiHit, abs(response)*itiMiss).at(response);

	% Stim stays on until the end of the ITI
	stimOff = threshold.delay(iti);

	%% Visual stimulus

	% Azimuth control
	% 1) stim fixed in place until interactive on
	% 2) wheel-conditional during interactive  
	% 3) fixed at response displacement azimuth after response
	azimuth = cond( ...
	    events.newTrial.to(interactiveOn), startingAzimuth*trialData.trialSide, ...
	    interactiveOn.to(response), startingAzimuth*trialData.trialSide + stimDisplacement, ...
	    response.to(events.newTrial), ...
	    startingAzimuth*trialData.trialSide.at(interactiveOn) + sign(stimDisplacement.at(response))*responseDisplacement);

	stim = vis.grating(t, 'square', 'gaussian');
	stim.sigma = sigma;
	stim.spatialFreq = spatialFreq;
	stim.phase = 2*pi*events.newTrial.map(@(v)rand);
	stim.azimuth = azimuth;
	stim.contrast = trialContrast.at(events.newTrial);
	stim.show = stimOn.to(stimOff);

	visStim.stim = stim;

	%% Display and save

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
	events.isRewarded = trialData.isRewarded.at(response);
	events.totalWater = water.scan(@plus, 0).map(fun.partial(@sprintf, '%.1fµl'));

	events.flipRate = trialData.flipRate;
	events.performanceTracker = trialData.performanceTracker;
	events.biasTracker = trialData.biasTracker;

end

function trialDataInit = initializeTrialData(expRef, ...
    contrasts, sigmoidParamsFlip, performanceTracker,alphaPerformanceTracker, ...
    sigmoidParamsBias, biasTracker, alphaBiasTracker, ...
    rewardProb, flipRate, rewardSize, wheelGain)

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
	trialDataInit.sigmoidParamsBias = sigmoidParamsBias;
	trialDataInit.biasTracker = biasTracker;
	trialDataInit.alphaBiasTracker = alphaBiasTracker;
	trialDataInit.flipRate = flipRate; %don't need it here?
	trialDataInit.rewardProb = rewardProb;
	trialDataInit.rewardSize = rewardSize; %don't need it here?
	trialDataInit.wheelGain = wheelGain; %don't need it here?

	% Set the first trial side randomly
	trialDataInit.trialSide = randsample([-1,1],1)
	% trialDataInit.trialSide = 1
	% Set up the flag for repeating incorrect
	trialDataInit.repeatTrial = false;
	% Initialize hit/miss
	trialDataInit.isRewarded = false;
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
            
	            if isempty(previousBlock) || ~isfield(previousBlock, 'biasTrackerValues')
	                lastBiasTracker = biasTracker;
	            else
	                lastBiasTracker = previousBlock.biasTrackerValues(end);
	            end
            
	            if isempty(previousBlock) || ~isfield(previousBlock, 'trialSideValues')
	                lastTrialSide = trialDataInit.trialSide;
	            else
	                lastTrialSide = previousBlock.trialSideValues(end);
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
	    trialDataInit.biasTracker = lastBiasTracker;
    
	    sigmoidOutBias = trialDataInit.sigmoidParamsBias(4) ...
	    + (trialDataInit.sigmoidParamsBias(3)-trialDataInit.sigmoidParamsBias(4))...
	    /(1+exp(-(trialDataInit.biasTracker-trialDataInit.sigmoidParamsBias(1))/trialDataInit.sigmoidParamsBias(2)));
	    randomFlip = rand
	    randomBias = rand
	    if randomFlip < trialDataInit.flipRate
	        %now can possibly change side
	        nextSide = ((randomBias - sigmoidOutBias)<0)*2-1 %if rightward bias, sigmoidOutBias~=1, nextSide<-1 -> L
	        trialDataInit.trialSide = nextSide;
	    else
	        trialDataInit.trialSide = lastTrialSide;
	    end
    
	else
	    % If this animal has no previous experiments, initialize performance 
	    % Initialize water reward size & wheel gain
	    trialDataInit.flipRate = flipRate;
	    trialDataInit.performanceTracker = performanceTracker;
	    trialDataInit.biasTracker = biasTracker;
	    trialDataInit.rewardSize = rewardSize;
	    trialDataInit.wheelGain = wheelGain;
	end
end

function trialData = updateTrialData(trialData,responseData)
% Update the performance 
% if isnan(responseData)
%     return;
% end

	stimDisplacement = responseData;

	% Next contrast
	trialData.trialContrast = randsample(trialData.contrasts,1); 

	%%%% Define response type based on trial condition
	r = rand;
	isCorrectChoice = stimDisplacement*trialData.trialSide < 0;
	% trialData.hit = isCorrectChoice && r <= trialData.rewardProb; %correct side and reward
	% trialData.miss = stimDisplacement*trialData.trialSide > 0 || ... %wrong side
	%     (stimDisplacement*trialData.trialSide < 0 && r > trialData.rewardProb); %correct side and no reward
	trialData.isRewarded = isCorrectChoice && r <= trialData.rewardProb;
	trialData.hit = isCorrectChoice;
	trialData.miss = ~isCorrectChoice;

	%update flipRate
	trialData.performanceTracker = trialData.alphaPerformanceTracker*isCorrectChoice ...
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

	%update Bias
	%biasTracker: 0,left"ward" bias; 1,right"ward" bias, 
	trialData.biasTracker = trialData.alphaBiasTracker*(stimDisplacement>0) ...
	    + (1-trialData.alphaBiasTracker)*trialData.biasTracker;
	trialData.sigmoidParamsBias
	sigmoidOutBias = trialData.sigmoidParamsBias(4) ...
	    + (trialData.sigmoidParamsBias(3)-trialData.sigmoidParamsBias(4))...
	    /(1+exp(-(trialData.biasTracker-trialData.sigmoidParamsBias(1))/trialData.sigmoidParamsBias(2)))

	%%%% Set flag to repeat - skip trial choice if so
	% if trialData.miss && ...
	%         ismember(trialData.trialContrast,trialData.contrasts(trialData.repeatOnMiss))

	trialData.repeatTrial = false;



	%%%% Pick next side (this is done at random)
	randomFlip = rand
	randomBias = rand
	if randomFlip < trialData.flipRate
	    %now can possibly change side
	    nextSide = ((randomBias - sigmoidOutBias)<0)*2-1 %if rightward bias, sigmoidOutBias~=1, nextSide<-1 -> L
	    trialData.trialSide = nextSide
	end

end

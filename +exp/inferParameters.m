function parsStruct = inferParameters(expdef)
%EXP.INFERPARAMETERS Infers the parameters required for experiment
%   Detailed explanation goes here

% create some signals just to pass to the definition function and track
% which parameter names are used

if ischar(expdef) && file.exists(expdef)
  expdeffun = fileFunction(expdef);
else
  expdeffun = expdef;
  expdef = which(func2str(expdef));
end

net = sig.Net;
e = struct;
e.t = net.origin('t');
e.events = net.subscriptableOrigin('events');
e.pars = net.subscriptableOrigin('pars');
e.pars.CacheSubscripts = true;
e.visual = net.subscriptableOrigin('visual');
e.audio = net.subscriptableOrigin('audio');
e.audio.SampleRate = 44100;
e.audio.NChannels = 2;
e.inputs = net.subscriptableOrigin('inputs');
e.outputs = net.subscriptableOrigin('outputs');

try
  expdeffun(e.t, e.events, e.pars, e.visual, e.inputs , e.outputs);
  % paramNames will be the strings corresponding to the fields of e.pars
  % that the user tried to reference in her expdeffun.
  paramNames = e.pars.Subscripts.keys';
  %The paramValues are signals corresponding to those parameters and they
  %will all be empty, except when they've been given explicit numerical
  %definitions right at the end of the function - and in that case, we'll
  %take those values (extracted into matlab datatypes, from the signals,
  %using .Node.CurrValue) to be the desired default values.
  paramValues = e.pars.Subscripts.values';
  parsStruct = cell2struct(cell(size(paramNames)), paramNames);
  for i = 1:size(paramNames,1)
      parsStruct.(paramNames{i}) = paramValues{i}.Node.CurrValue;
  end
  sz = iff(isempty(fieldnames(parsStruct)), 1,... % if there are no paramters sz = 1
      structfun(@(a)size(a,2), parsStruct)); % otherwise get number of columns
  parsStruct.numRepeats = ones(1,max(sz)); % add 'numRepeats' parameter
  parsStruct.defFunction = expdef;
  parsStruct.type = 'custom';
  % Define the ExpPanel to use (automatically by name convention for now)
  [path, name, ext] = fileparts(expdef);
  ExpPanel_name = [name 'ExpPanel'];
  ExpPanel_fn = [path filesep ExpPanel_name ext];
  if exist(ExpPanel_fn,'file'); parsStruct.expPanelFun = ExpPanel_name; end
catch ex
  net.delete();
  rethrow(ex)
end

net.delete();


end
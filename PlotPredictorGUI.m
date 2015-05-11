function PlotPredictorGUI(varargin)
% PLOTPREDICTORSGUI
if nargin<1, data = []; end
easydefaults('TR', 1, 'nCond', 1, 'nTrials', 1, 'trialDur', 0, 'ITI', 20, 'pad', 20, 'HPF', 128, 'JITTER', 0, 'AUTOSAVEPLOT', 0);
easyparse(varargin, {'TR' 'nCond' 'nTrials' 'trialDur' 'ITI' 'HPF' 'pad' 'JITTER' 'AUTOSAVEPLOT'});
data.screensize = get(0,'Screensize');  % Get screensize.
data.plotW2Hratio = 2.5; 
data.TR = TR; 
data.nCond = nCond;
data.nTrials = nTrials;
data.trialDur = trialDur;
data.ITI = ITI;
data.HPF = HPF;
data.pad = 20; 
data.fontax = 16;
data.fontlabels = 16; 
data.fontbutton = 15;
data.osrate = 1; 
S.fh = figure('numbertitle','off',...
              'menubar','none',...
              'units','normal',...
              'position',[.025 .50 .25 .35],...
              'name','View Predictor GUI',...
              'resize','on');
guidata(S.fh, data); 
% =========================================================================
% MENUS
% =========================================================================
S.menu(1) = uimenu(S.fh, 'Label', 'Menu');
S.menu(2) = uimenu(S.menu(1), 'Label', 'Export Timeseries Data', 'Enable', 'off', 'Callback', {@menu_exporttoworkspace, S});
S.menu(3) = uimenu(S.fh, 'Label', 'Font Size', 'Separator', 'on');
S.menu(4) = uimenu(S.menu(3), 'Label', 'Increase', 'Accelerator', 'i', 'Callback', {@menu_changefontsize, S});
S.menu(5) = uimenu(S.menu(3), 'Label', 'Decrease', 'Accelerator', 'd', 'Separator', 'on', 'Callback', {@menu_changefontsize, S});
S.menu(6) = uimenu(S.menu(1), 'Label', 'Exit', 'Separator', 'on', 'Callback', {@pb_closegui, S});

% =========================================================================
% PANELS
% =========================================================================
S.inputpanel = uipanel('parent', S.fh, 'units', 'norm', 'pos', [.025 .025 .425 .95], 'backg', [1 1 1], 'foreg', [0 0 0], 'borderw', 1);
S.buttonpanel = uipanel('Parent', S.fh, 'units', 'norm', 'pos', [.500 .025 .475 .95], 'backg', [1 1 1], 'foreg', [0 0 0], 'borderw', 1);

% =========================================================================
% INPUT FIELDS
% =========================================================================
W1 = {'parent', S.inputpanel,'style','text','units','normal', 'clip', 'off', 'visible', 'on','backg', [1 1 1], 'fontsize',data.fontlabels,'fontweight','bold','position'};
W2 = {'parent', S.inputpanel,'style','edit','units','normal', 'clip', 'off', 'visible', 'on','fontsize',20, 'position'};
inputname = { 'TR (s)' 'Number of Trials' 'Trial Duration (s)' 'Inter-Trial Interval (s)' 'Begin/End Padding (s)' 'High-Pass Filter (s)'};
defvalues = [data.TR data.nTrials data.trialDur data.ITI data.pad data.HPF];
eh = 1/(length(inputname)*2);
ew = .80;
pos = 1-eh:-eh:0;
pos = reshape(pos, 2, length(inputname));
yadj = [-.010 .020];
if regexp(computer, 'PCWIN'), yadj = [0 0]; end
for i = 1:length(inputname)
    data.os = computer;
    S.tx(i) = uicontrol(W1{:}, [(1-ew)/2 pos(1,i)+yadj(1) ew eh], 'string', inputname{i});
    S.ed(i) = uicontrol(W2{:}, [(1-ew)/2 pos(2,i)+yadj(2) ew eh], 'string', num2str(defvalues(i))); 
end

% =========================================================================
% PUSH BUTTONS
% =========================================================================
T = {'parent', S.buttonpanel,'style','toggle','units','normal', 'clip', 'off', 'visible', 'on','fontsize',data.fontbutton,'fontweight','bold','position'};
togglename = {'Plot Raw' 'Plot Convolved' 'Plot Filtered' 'Jitter Onsets' 'Show Filtered Design Matrix' 'QUIT'};
eh = 1/(length(togglename)*2);
ew = .80;
pos = 1-eh:-eh:0;
pos = pos(1:2:end) - eh/2;
eh = eh*1.25;
for i = 1:length(togglename)
    S.tg(i) = uicontrol(T{:}, [(1-ew)/2 pos(i) ew eh], 'String', togglename{i});
end
set(S.tg(1), 'callback', {@pb_plot, S});
set(S.tg(4), 'enable', 'off', 'callback', {@pb_jitter,S});
set(S.tg(2), 'enable', 'off', 'callback', {@pb_conv,S});
set(S.tg(3), 'enable', 'off', 'callback', {@pb_filter,S});
set(S.tg(5), 'enable', 'off', 'callback', {@pb_xmatrix,S});
set(S.tg(6), 'enable', 'on', 'callback', {@pb_closegui,S});

if AUTOSAVEPLOT
    
    time = strtrim(datestr(now,'HHMMSS'));
    day = strtrim(datestr(now,'mmm_DD'));
    name = ['screen' time '_guiplot_' day '.png'];
    export_fig(name, '-nocrop', S.fh);
    name = ['screen' time '_predictplot_' day '.png'];
    pb_plot([],[],S); export_fig(name, '-nocrop');
    if JITTER
        name = ['screen' time '_predictjitterplot_' day '.png'];
        pb_jitter([],[],S); export_fig(name, '-nocrop');
    end
    name = ['screen' time '_convplot_' day '.png'];
    pb_conv([],[],S); export_fig(name, '-nocrop');
    name = ['screen' time '_filterplot_' day '.png'];
    pb_filter([],[],S); export_fig(name, '-nocrop');
    name = ['screen' time '_xmatrixplot_' day '.png'];
    pb_xmatrix([],[],S); export_fig(name, '-nocrop');
    
end
    

end 
% =========================================================================
% CALLBACKS
% =========================================================================
function pb_plot(varargin)
    S = varargin{3};
    params = get(S.ed(:), 'String');
    if any(cellfun('length', params)==0), errordlg('Missing Inputs!', 'Error'); return; end
    params = cell2mat(cellfun(@str2num, params, 'Unif', false));
    data = guidata(S.fh);
    data.nCond = 1;
    data.TR = params(1); 
    data.nTrials = params(2);
    data.trialDur = params(3);
    data.ITI = params(4);
    data.pad = params(5);
    data.HPF = params(6);
    if data.nTrials==1
        data.scan_length = (2*data.pad)+(data.trialDur);
    else
        data.scan_length = (data.nTrials-1)*(data.trialDur+data.ITI)+(2*data.pad);
    end
    numTR=ceil(data.scan_length/data.TR);
    data.numTR = numTR; 
    unisample = repmat(data.ITI + data.trialDur, data.nTrials*data.nCond, 1);
    unisample(1) = data.pad; 
    onsvector = cumsum(unisample);
    if data.nCond > 1
        tmp = randperm(length(unisample));
        for i = 1:data.nCond
            onsets(:,i) = sort(onsvector(tmp(i:data.nCond:end)));
            [X(:,i), X0(:,i), FX(:,i)] = make_regressor(numTR, data.TR, 30, onsets(:,i), data.trialDur, [], 1, 0, data.HPF);
        end
    else
        onsets = onsvector;
        [X,X0,FX] = make_regressor(numTR, data.TR, 16, onsets, data.trialDur, [], 1, 0, data.HPF);
    end
    data.unisample = unisample; 
    data.X = X - X(1);
    data.X0 = X0;
    data.FX = FX;
    data.axis = [0 data.numTR+1 -.15 1.05];
    tmppos = get(S.fh, 'position');
    figh = tmppos(4)*.75; 
    figw = figh*data.plotW2Hratio; 
    lr = tmppos(1)+tmppos(3); 
    tb = tmppos(2) + tmppos(4) - figh; 
    if lr + figw > 1, lr = tmppos(1)-figw; end
    if lr < 0
        lr = tmppos(1) - (figw/2);
        tb = tmppos(2) - figh; 
        if tb < 0, tb = tmppos(2) + tmppos(4); end
    end
    figpos = [lr tb figw figh];
    data.figh = figure('color', 'white', 'units', 'normal', 'position', figpos,'menubar','none', 'name', 'Predictor Plot', 'visible', 'off'); data.plotax = gca;
    if data.nCond > 1
        for i = 1:data.nCond
            subplot(2, 1, i); 
            plot(X0(:,i), '-', 'LineWidth', 2, 'Color', [0 0 0], 'tag', 'raw'); hold on;
            set(gca, 'FontName', 'Arial', 'FontSize', data.fontax);
            xlabel(sprintf('Condition %d', i), 'fontsize', data.fontlabels);
            ylabel('Activation', 'fontsize', data.fontlabels);
            axis([0 length(X) + 1 -.15 1.05]);
            box('off'); 
        end    
    else
        plot(data.X0, '-', 'LineWidth', 2, 'Color', [0 0 0], 'tag', 'raw'); hold on;
        set(gca, 'FontName', 'Arial', 'FontSize', data.fontax);
        xlabel(data.plotax, 'TR', 'fontname', 'Arial', 'fontsize', data.fontlabels, 'fontweight', 'bold');
        ylabel(data.plotax, 'Predicted Response', 'fontsize', data.fontlabels, 'fontweight', 'bold');
        mytitle = sprintf('Single Predictor - %d Trial(s), Trial Duration = %d s, ITI = %d s, TR = %2.3f s', data.nTrials, data.trialDur, data.ITI, data.TR);
        title(data.plotax, mytitle, 'fontsize', ceil(data.fontlabels*1.10), 'fontweight', 'bold');
        caxis = data.axis; 
        caxis(2) = caxis(2);
        axis(data.plotax, caxis);
        box(data.plotax, 'off');
        set(data.figh, 'visible', 'on');
    end
    figure(data.figh)
    guidata(S.fh, data);
    set(S.tg(:),'enable', 'on', 'backg', [.93 .93 .93]);
    if data.nTrials==1, set(S.tg([4]), 'enable', 'off'); end
    set(S.menu(2), 'enable', 'on');
end
function pb_conv(varargin)
    S = varargin{3};
    data = guidata(S.fh);
    ch = findobj(data.plotax, 'tag', 'conv');
    if isempty(ch)
        axes(data.plotax)
        plot(data.X, '--', 'LineWidth', 4, 'Color', 'r', 'tag', 'conv');
    else
        set(ch, 'ydata', data.X);
    end
    figure(data.figh)
end
function pb_filter(varargin)
    S = varargin{3};
    data = guidata(S.fh);
    HPF = str2num(get(S.ed(6), 'String'));
    FX = bspm_filter(data.X, data.TR, HPF);
    FX = scalematrix(FX);
    data.FX = FX; 
    ch = findobj(data.plotax, 'tag', 'filt');
    if isempty(ch)
        axes(data.plotax)
        plot(FX, '--', 'LineWidth', 4, 'Color', 'b', 'tag', 'filt');
    else
        set(ch, 'ydata', FX);
    end
    guidata(S.fh, data);
    figure(data.figh)
end
function pb_jitter(varargin)
    S = varargin{3};
    data = guidata(S.fh);
    
    jitsample = getjit(ceil(data.ITI*.5), data.ITI, data.nTrials*data.nCond);
    jitsample(1) = data.pad; 
    [JX, JX0] = make_regressor(data.numTR, data.TR, 16, cumsum(jitsample+data.trialDur), data.trialDur, [], 1, 0, data.HPF);
    JX = scalematrix(JX);
    data.X = JX - JX(1);
    data.X0 = JX0;
    tags = {'raw' 'conv' 'filt'};
    for i = 1:length(tags)
        ch = findobj(data.plotax, 'tag', tags{i});
        if ~isempty(ch), set(ch, 'ydata', [], 'xdata', []); end
    end
    ch = findobj(data.plotax, 'tag', 'raw');
    set(ch, 'ydata', data.X0);
    guidata(S.fh, data);
    figure(data.figh)
end
function pb_xmatrix(varargin)

    S = varargin{3};
    data = guidata(S.fh);
    FX = data.FX;
    X = [scalematrix(FX) ones(size(FX,1), 1)];
    tmppos = get(S.fh, 'position');
    figh = tmppos(4)*.75; 
    figw = tmppos(3)*.75; 
    lr = tmppos(1)+tmppos(3); 
    tb = tmppos(2) + tmppos(4) - figh; 
    if lr + figw > 1, lr = tmppos(1)-figw; end
    if lr < 0
        lr = tmppos(1) - (figw/2);
        tb = tmppos(2) - figh; 
        if tb < 0, tb = tmppos(2) + tmppos(4); end
    end
    figpos = [lr tb figw figh];
    data.figxmat = figure('color', 'white', 'units', 'normal', 'position', figpos,'menubar','none', 'name', 'Design Matrix'); data.xmatax = gca;
    imagesc(X); colormap('gray');
    set(gca, 'FontName', 'Arial');
    cn = {'Predictor' 'Constant'};
    set(gca, 'FontName', 'Arial', 'FontSize', data.fontax);
    ylabel(data.plotax, 'TR', 'fontsize', data.fontlabels, 'fontweight', 'bold');
    title('Design Matrix', 'fontsize', ceil(data.fontlabels*1.10), 'fontweight', 'bold');
    set(gca, 'XTick', 1:length(cn));
    set(gca, 'XTickLabel', cn, 'fontsize', data.fontlabels);
    ylabel('TR', 'fontsize', data.fontlabels);
    set(S.tg(:),'backg', [.93 .93 .93]);
    guidata(S.fh, data);
    figure(S.fh)
    figure(data.figxmat)
end
function pb_closegui(varargin)

if isempty(gcbf)
   if length(dbstack) == 1
      warning(message('MATLAB:closereq:ObsoleteUsage'));
   end
   close('force');
else
   delete(gcbf);
   delete(get(0,'Children'))
end
%     if length(varargin)==3
%         S = varargin{3};
%         h = S.fh;
%     else
%         h = varargin{1};
%     end
%     data = guidata(S.fh);
%     if isfield(data, 'plotax'), delete(get(data.plotax, 'parent')); end
%     if isfield(data, 'xmatax'), delete(get(data.xmatax, 'parent')); end
%     delete(h); % Bye-bye figure
end
function menu_changefontsize(varargin)
    S = varargin{3};
    data = guidata(S.fh);
    if strcmp(get(varargin{1}, 'Label'), 'Increase')
        F = 1.1;

    else
        F = 0.9;
    end
    data.fontax = data.fontax*F;
    data.fontlabels = data.fontlabels*F; 
    data.fontbutton = data.fontbutton*F; 
    h = findall(S.fh, 'type', 'uicontrol');
    fs = get(h, 'fontsize');
    for i = 1:length(h)
        set(h(i), 'fontsize', fs{i}*F);
    end
    guidata(S.fh, data);
end
function menu_exporttoworkspace(varargin)
    S = varargin{3};
    tmp = guidata(S.fh);
    output.TR = tmp.TR;
    output.X0 = tmp.X0;
    output.X = tmp.X;
    putvar('output');
    msgbox('Timeseries data exported to data structure ''output''', 'Success');
    fprintf('\n\nExported variable ''output'' has the following fields''\n');
    structdisp(output);
end
% =========================================================================
% SUB-FUNCTIONS
% =========================================================================
function [X, X0, FX] = make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag, cutoff)
% BOB_SPM_MAKE_REGRESSOR
%
%   USAGE: X = bob_spm_make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag)
%   
%   ARGUMENTS
%       nvols = # of volumes
%       TR = TR (in secs)
%       TRbin = # of time bins per scan (for oversampling)
%       ons = onsets (in secs)
%       dur = durations (in secs) - if 0, model all as event
%       pm = parameters of interest (each column is a different parameter)
%       TRons = bin to use as onset
%       TDtag = include temporal derivative
%
%   OUTPUT
%       X       = convolved X matrix
%       X0      = un-convolved X matrix
%
% =========================================================
if nargin<8, TDtag = 0; end
if nargin<7, TRons = 1; end   
if nargin<6, pm = []; end
if nargin<5, error('USAGE: X = bspm_make_regressor(nvols, TR, TRbin, ons, dur, pm, TRons, TDtag)'); end
if size(ons,1)==1 && size(ons,2)~=1, ons = ons'; end
if size(dur,1)==1 && size(dur,2)~=1, dur = dur'; end
if mod(TRbin, TR), TRbin = ceil(TRbin/TR)*TR; end
osrate = TR/(TR/TRbin);
xlength = nvols*osrate;
onso = fix((ons/TR)*osrate)+TRons;
duro = fix((dur/TR)*osrate);
if duro==0, duro = 1; end
if length(duro)==1, duro = repmat(duro,length(onso),1); end
hrf = spm_hrf(TR/TRbin);
xo = zeros(xlength,1);
for i = 1:length(onso)
    xo(onso(i):onso(i)+duro(i)-1) = 1;
end
xoc = conv(xo, hrf);
xoc = xoc(1:xlength);

X = scalematrix(xoc(1:osrate:xlength));
X0 = scalematrix(xo(1:osrate:xlength));
FX = scalematrix(bspm_filter(X, TR, cutoff));
if ~isempty(pm)
    pm = scalematrix(pm);
    po = repmat(xo,1,size(pm,2));
    for p = 1:size(pm,2)
        cpm = pm(:,p);
        cpm = cpm - mean(cpm);
        for i = 1:length(onso)
            po(onso(i):onso(i)+duro(i)-1,p) = 1*cpm(i);
        end
        poc(:,p) = conv(po(:,p),hrf);
    end
    poc = poc(1:xlength,:);
    X(:,2:2+size(pm,2)-1) = poc(1:osrate:xlength,:);
    X0(:,2:2+size(pm,2)-1) = po(1:osrate:xlength,:);
end

if TDtag
    X1 = X;
    X01 = X0;
    clear xo xoc X X0 po cpm po poc 
    hrf = spm_hrf_td(TR/TRbin);
    xo = zeros(xlength,1);
    for i = 1:length(onso)
        xo(onso(i):onso(i)+duro(i)-1) = 1;
    end
    xoc = conv(xo, hrf);
    xoc = xoc(1:xlength);
    X = xoc(1:osrate:xlength);
    X0 = xo(1:osrate:xlength);
    if ~isempty(pm)
        pm = scalematrix(pm);
        po = repmat(xo,1,size(pm,2));
        for p = 1:size(pm,2)
            cpm = pm(:,p);
            cpm = cpm - mean(cpm);
            for i = 1:length(onso)
                po(onso(i):onso(i)+duro(i)-1,p) = 1*cpm(i);
            end
            poc(:,p) = conv(po(:,p),hrf);
        end
        poc = poc(1:xlength,:);
        X(:,2:2+size(pm,2)-1) = poc(1:osrate:xlength,:);
        X0(:,2:2+size(pm,2)-1) = po(1:osrate:xlength,:);
    end
    X2 = X;
    X02 = X0;
    X = zeros(size([X1 X2]));
    X0 = zeros(size([X01 X02]));
    X(:,1:2:end) = X1;
    X(:,2:2:end) = X2;
    X0(:,1:2:end) = X01;
    X0(:,2:2:end) = X02;
    
end
end
function SX = scalematrix(X)
    SX = X;
    [r, c] = size(X);
    mn = min(X);
    mx = max(X);
    for i = 1:c
        SX(:,i) = (X(:,i) - mn(i))/(mx(i)-mn(i));
    end
end
function FX = bspm_filter(X, TR, cutoff)
    K.RT = TR;
    K.HParam = cutoff;
    K.row = 1:length(X);
    K = spm_filter(K);
    FX = spm_filter(K, X);
end
function jitsample = getjit(minSOA, meanSOA, nTrials)
    goodjit = 0; 
    while ~goodjit
        jitsample = minSOA + poissrnd(meanSOA-minSOA, nTrials, 1);
        if round(mean(jitsample)*100)==meanSOA*100, goodjit = 1; end
    end 
end
function easydefaults(varargin)
% easydefaults  Set many default arguments quick and easy.
%
%   - For input arguments x1,x2,x3, set default values x1def,x2def,x3def
%     using easydefaults as parameter-value pairs:
%       easydefaults('x1',x1def,'x2',x2def,'x3',x3def);
%   
%   - Defaults can be set for any input argument, whether explicit or as 
%     part of a parameter-value pair:
%       function dummy_function(x,varargin)
%           easydefaults('x',1,'y',2);
%           ...   
%       end
%
%   - easydefaults and easyparse can in principle be used in either order, 
%     but it is usually better to parse first and fill in defaults after:
%       function dummy_function(x,varargin)
%           easyparse(varargin,'y')
%           easydefaults('x',1,'y',2);
%           ...   
%       end
%
%   CAVEAT UTILITOR: this function relies on evals and assignin statements.
%   Input checking is performed to limit potential damage, but use at your 
%   own risk.
%
%   Author: Jared Schwede 
%   Last update: Jan 14, 2013

    % Check that all inputs come in parameter-value pairs.
    if mod(length(varargin),2)
        error('Default arguments must be specified in pairs!');
    end
    
    for i=1:2:length(varargin)
        if ~ischar(varargin{i})
            error('Variables to easydefaults must be written as strings!');
        end
        
        % We'll check that the varargin is a valid variable name. This
        % should hopefully avoid any nasty code...
        if ~isvarname(varargin{i})
            error('Invalid variable name!');
        end
        
        if exist(varargin{i},'builtin') || (exist(varargin{i},'file') == 2) || exist(varargin{i},'class')
            warning('MATLAB:defined_function',['''' varargin{i} ''' conflicts with the name of a function, m-file, or class along the MATLAB path and will be ignored by easydefaults.' ...
                                        ' Please rename the variable, or use a temporary variable with easydefaults and explicitly define ''' varargin{i} ...
                                        ''' within your function.']);
        else
            if ~evalin('caller',['exist(''' varargin{i} ''',''var'')'])
                % We assign the arguments to a struct, s, which allows us to
                % check that the evalin statement will not either throw an 
                % error or execute some nasty code.
                s.(varargin{i}) = varargin{i+1};
                assignin('caller',varargin{i},varargin{i+1});
            end
        end
    end
end
function s = easyparse(caller_varargin,allowed_names)
% easyparse    Parse parameter-value pairs without using inputParser
%   easyparse is called by a function which takes parameter value pairs and
%   creates individual variables in that function. It can also be used to
%   generate a struct like inputParser.
%
%   - To create variables in the function workspace according to the
%     varargin of parameter-value pairs, use this syntax in your function:
%       easyparse(varargin)
%
%   - To create only variables with allowed_names, create a cell array of
%     allowed names and use this syntax:
%       easyparse(varargin, allowed_names);
%
%   - To create a struct with fields specified by the names in varargin,
%     (similar to the output of inputParser) ask for an output argument:
%       s = easyparse(...);
%  
%   CAVEAT UTILITOR: this function relies on assignin statements. Input
%   checking is performed to limit potential damage, but use at your own 
%   risk.
%
%   Author: Jared Schwede
%   Last update: January 14, 2013

    % We assume all inputs come in parameter-value pairs. We'll also assume
    % that there aren't enough of them to justify using a containers.Map. 
    for i=1:2:length(caller_varargin)
        if nargin == 2 && ~any(strcmp(caller_varargin{i},allowed_names))
            error(['Unknown input argument: ' caller_varargin{i}]);
        end
        
        if ~isvarname(caller_varargin{i})
            error('Invalid variable name!');
        end
        
        
        % We assign the arguments to the struct, s, which allows us to
        % check that the assignin statement will not either throw an error 
        % or execute some nasty code.
        s.(caller_varargin{i}) = caller_varargin{i+1};
        % ... but if we ask for the struct, don't write all of the
        % variables to the function as well.
        if ~nargout
            if exist(caller_varargin{i},'builtin') || (exist(caller_varargin{i},'file') == 2) || exist(caller_varargin{i},'class')
                warning('MATLAB:defined_function',['''' caller_varargin{i} ''' conflicts with the name of a function, m-file, or class along the MATLAB path and will be ignored by easyparse.' ...
                                            ' Please rename the variable, or use a temporary variable with easyparse and explicitly define ''' caller_varargin{i} ...
                                            ''' within your function.']);
            else
                assignin('caller',caller_varargin{i},caller_varargin{i+1});
            end
        end
    end
end
function putvar(varargin)
% Assigns variables from the current workspace down into the base MATLAB workspace
% usage: putvar(var1)
% usage: putvar(var1,var2, ...)
%
% putvar moves variables from the current matlab
% workspace down to the base matlab workspace.
% putvar is something that can be used while
% in a debugging session, to retain the value
% of a variable, saving it into the base matlab
% workspace. putvar is also a way to return a
% specific variable, avoiding the use of a
% return argument (for whatever reason you might
% have.)
%
% putvar cannot assign variables that are not
% already in existence in the caller workspace.
%
%
% arguments: (input)
%  var1, var2 ... - Matlab variables in the
%       current workspace, that will then be
%       assigned into the base matlab workspace.
%
%       Alternatively, you can supply the names
%       of these variables, as strings.
%
%       If a variable with that name already
%       exists in the base workspace, it will
%       be overwritten, with a warning message
%       generated. That warning message can be
%       disabled by the advance command:
%
%       warning('off','PUTVAR:overwrite')
%
%
% Example:
% % First, save the function testputvar.m on your search path.
% % in MATLAB itself, try this.
%
% ==================
% function testputvar
% A = 3;
% B = 23;
% C = pi/2;
% D = 'The quick brown fox';
% putvar(A,'C',D)
% ==================
% 
% % Next, clear your workspace. Clear ensures that no
% % variables exist initially in the base workspace. Then
% % run the function testputvar, and finally execute the
% % who command, all at the command line.
%
% >> clear
% >> testputvar
% >> who
% 
% % Your variables are:
% % A  C  D
%
% % The output from who tells it all. Inside the function
% % testputvar, we had defined four variables, A, B, C and D.
% % But after testputvar terminates, its variables will fall
% % into the bit bucket, disappearing into limbo. The putvar
% % call inside testputvar ensures that the variables A, C
% % and D are returned into the base workspace, yet no return
% % arguments were provided.
%
%
% See also: uigetvar, uigetdir, uigetfile, who, whos
%
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 2.0
% Release date: 4/08/2010

if nargin < 1
  % no variables requested for the assignment,
  % so this is a no-op
  return
end

% how many variables do we need to assign?
nvar = numel(varargin);

% get the list of variable names in the caller workspace.
% callervars will be a cell array that lists the names of all
% variables in the caller workspace.
callervars = evalin('caller','who');

% likewise, basevars is a list of the names of all variables
% in the base workspace.
basevars = evalin('base','who');

% loop over the variables supplied
for i = 1:nvar
  % what was this variable called in the caller workspace?
  varname = inputname(i);
  vari = varargin{i};
  
  if ~isempty(varname)
    % We have a variable name, so assign this variable
    % into the base workspace
    
    % First though, check to see if the variable is
    % already there. If it is, we will need to set
    % a warning.
    if ismember(varname,basevars)
      warning('PUTVAR:overwrite', ...
        ['Input variable #',num2str(i),' (',varname,')', ...
        ' already exists in the base workspace. It will be overwritten.'])
    end
    
    % do the assign into the indicated name
    assignin('base',varname,varargin{i})
    
  elseif ischar(vari) && ismember(vari,callervars)
    % the i'th variable was a character string, that names
    % a variable in the caller workspace. We can assign
    % this variable into the base workspace.
    
    % First though, check to see if the variable is
    % already there. If it is, we will need to set
    % a warning.
    varname = vari;
    if ismember(varname,basevars)
      warning('PUTVAR:overwrite', ...
        ['Input variable #',num2str(i),' (',varname,')', ...
        ' already exists in the base workspace. It will be overwritten.'])
    end
    
    % extract the indicated variable contents from
    % the caller workspace.
    vari = evalin('caller',varname);
    
    % do the assign into the indicated name
    assignin('base',varname,vari)
    
  else
    % we cannot resolve this variable
    warning('PUTVAR:novariable', ...
      ['Did not assign input variable #',num2str(i), ...
      ' as no caller workspace variable was available for that input.'])
    
  end
  
end
end
function strucdisp(Structure, depth, printValues, maxArrayLength, fileName)
%STRUCDISP  display structure outline
%
%   STRUCDISP(STRUC, DEPTH, PRINTVALUES, MAXARRAYLENGTH, FILENAME) displays
%   the hierarchical outline of a structure and its substructures. 
%
%   STRUC is a structure datatype with unknown field content. It can be 
%   either a scalar or a vector, but not a matrix. STRUC is the only
%   mandatory argument in this function. All other arguments are optional.
%
%   DEPTH is the number of hierarchical levels of the structure that are
%   printed. If DEPTH is smaller than zero, all levels are printed. Default
%   value for DEPTH is -1 (print all levels).
%
%   PRINTVALUES is a flag that states if the field values should be printed
%   as well. The default value is 1 (print values)
%
%   MAXARRAYLENGTH is a positive integer, which determines up to which
%   length or size the values of a vector or matrix are printed. For a
%   vector holds that if the length of the vector is smaller or equal to
%   MAXARRAYLENGTH, the values are printed. If the vector is longer than
%   MAXARRAYLENGTH, then only the size of the vector is printed.
%   The values of a 2-dimensional (m,n) array are printed if the number of
%   elements (m x n) is smaller or equal to MAXARRAYLENGTH.
%   For vectors and arrays, this constraint overrides the PRINTVALUES flag.
%
%   FILENAME is the name of the file to which the output should be printed.
%   if this argument is not defined, the output is printed to the command
%   window.
%
%   Contact author: B. Roossien <roossien@ecn.nl>
%   (c) ECN 2007-2008
%
%   Version 1.3.0


%% Creator and Version information
% Created by B. Roossien <roossien@ecn.nl> 14-12-2006
%
% Based on the idea of 
%       M. Jobse - display_structure (Matlab Central FileID 2031)
%
% Acknowledgements:
%       S. Wegerich - printmatrix (Matlab Central FileID 971)
%
% Beta tested by: 
%       K. Visscher
%
% Feedback provided by:
%       J. D'Errico
%       H. Krause
%       J.K. Kok
%       J. Kurzmann
%       K. Visscher
%
%
% (c) ECN 2006-2007
% www.ecn.nl
%
% Last edited on 08-03-2008



%% Version History
%
% 1.3.0 : Bug fixes and added logicals
% 1.2.3 : Buf fix - Solved multi-line string content bug
% 1.2.2 : Bug fix - a field being an empty array gave an error
% 1.2.1 : Bug fix
% 1.2.0 : Increased readability of code
%         Makes use of 'structfun' and 'cellfun' to increase speed and 
%         reduce the amount of code
%         Solved bug with empty fieldname parameter
% 1.1.2 : Command 'eval' removed with a more simple and efficient solution
% 1.1.1 : Solved a bug with cell array fields
% 1.1.0 : Added support for arrayed structures
%         Added small matrix size printing
% 1.0.1 : Bug with empty function parameters fixed
% 1.0.0 : Initial release



%% Main program
%%%%% start program %%%%%

    % first argument must be structure
    if ~isstruct(Structure)
        error('First input argument must be structure');
    end
    
    % first argument can be a scalar or vector, but not a matrix
    if ~isvector(Structure)
        error('First input argument can be a scalar or vector, but not a matrix');
    end
    
    % default value for second argument is -1 (print all levels)
    if nargin < 2 || isempty(depth)
        depth = -1;
    end

    % second argument must be an integer
    if ~isnumeric(depth)
        error('Second argument must be an integer');
    end

    % second argument only works if it is an integer, therefore floor it
    depth = floor(depth);
    
    % default value for third argument is 1
    if nargin < 3 || isempty(printValues)
        printValues = 1;
    end

    % default value for fourth argument is 10
    if nargin < 4 || isempty(maxArrayLength)
        maxArrayLength = 10;
    end

    
    % start recursive function   
    listStr = recFieldPrint(Structure, 0, depth, printValues, ... 
                            maxArrayLength);

    
    % 'listStr' is a cell array containing the output
    % Now it's time to actually output the data
    % Default is to output to the command window
    % However, if the filename argument is defined, output it into a file

    if nargin < 5 || isempty(fileName)
        
        % write data to screen
        for i = 1 : length(listStr)
            disp(cell2mat(listStr(i, 1)));
        end
        
    else
        
        % open file and check for errors
        fid = fopen(fileName, 'wt');
        
        if fid < 0
            error('Unable to open output file');
        end
        
        % write data to file
        for i = 1 : length(listStr)
            fprintf(fid, '%s\n', cell2mat(listStr(i, 1)));
        end
        
        % close file
        fclose(fid);
        
    end
    







%% FUNCTION: recFieldPrint
function listStr = recFieldPrint(Structure, indent, depth, printValues, ...
                                 maxArrayLength)


% Start to initialiase the cell listStr. This cell is used to store all the
% output, as this is much faster then directly printing it to screen.

listStr = {};


% "Structure" can be a scalar or a vector.
% In case of a vector, this recursive function is recalled for each of
% the vector elements. But if the values don't have to be printed, only
% the size of the structure and its fields are printed.

if length(Structure) > 1

    if (printValues == 0)

        varStr = createArraySize(Structure, 'Structure');

        listStr = [{' '}; {['Structure', varStr]}];

        body = recFieldPrint(Structure(1), indent, depth, ...
                             printValues, maxArrayLength);

        listStr = [listStr; body; {'   O'}];

    else

        for iStruc = 1 : length(Structure)

            listStr = [listStr; {' '}; {sprintf('Structure(%d)', iStruc)}];

            body = recFieldPrint(Structure(iStruc), indent, depth, ...
                                 printValues, maxArrayLength);

            listStr = [listStr; body; {'   O'}];

        end

    end

    return

end


%% Select structure fields
% The fields of the structure are distinguished between structure and
% non-structure fields. The structure fields are printed first, by
% recalling this function recursively.

% First, select all fields.

fields = fieldnames(Structure);

% Next, structfun is used to return an boolean array with information of
% which fields are of type structure.

isStruct = structfun(@isstruct, Structure);

% Finally, select all the structure fields

strucFields = fields(isStruct == 1);


%% Recursively print structure fields 
% The next step is to select each structure field and handle it
% accordingly. Each structure can be empty, a scalar, a vector or a matrix.
% Matrices and long vectors are only printed with their fields and not with
% their values. Long vectors are defined as vectors with a length larger
% then the maxArrayLength value. The fields of an empty structure are not
% printed at all.
% It is not necessary to look at the length of the vector if the values
% don't have to be printed, as the fields of a vector or matrix structure
% are the same for each element.

% First, some indentation calculations are required.

strIndent = getIndentation(indent + 1);
listStr = [listStr; {strIndent}];

strIndent = getIndentation(indent);

% Next, select each field seperately and handle it accordingly

for iField = 1 : length(strucFields)

    fieldName = cell2mat(strucFields(iField));
    Field =  Structure.(fieldName);
    
    % Empty structure
    if isempty(Field)

        strSize = createArraySize(Field, 'Structure');

        line = sprintf('%s   |--- %s :%s', ...
                       strIndent, fieldName, strSize);

        listStr = [listStr; {line}];

    % Scalar structure
    elseif isscalar(Field)

        line = sprintf('%s   |--- %s', strIndent, fieldName);

        % Recall this function if the tree depth is not reached yet
        if (depth < 0) || (indent + 1 < depth)
            lines = recFieldPrint(Field, indent + 1, depth, ...
                                  printValues, maxArrayLength);

            listStr = [listStr; {line}; lines; ...
                       {[strIndent '   |       O']}];
        else
            listStr = [listStr; {line}];
        end

    % Short vector structure of which the values should be printed    
    elseif (isvector(Field)) &&  ...
           (printValues > 0) && ...
           (length(Field) < maxArrayLength) && ...
           ((depth < 0) || (indent + 1 < depth))

        % Use a for-loop to print all structures in the array
        for iFieldElement = 1 : length(Field)

            line = sprintf('%s   |--- %s(%g)', ...
                           strIndent, fieldName, iFieldElement);

            lines = recFieldPrint(field(iFieldElement), indent + 1, ...
                                 depth, printValues, maxArrayLength);

            listStr = [listStr; {line}; lines; ...
                       {[strIndent '   |       O']}];

            if iFieldElement ~= length(Field)
                listStr = [listStr; {[strIndent '   |    ']}];
            end

        end

    % Structure is a matrix or long vector
    % No values have to be printed or depth limit is reached
    else

        varStr = createArraySize(Field, 'Structure');

        line = sprintf('%s   |--- %s :%s', ...
                       strIndent, fieldName, varStr);

        lines = recFieldPrint(Field(1), indent + 1, depth, ...
                              0, maxArrayLength);

        listStr = [listStr; {line}; lines; ...
                   {[strIndent '   |       O']}];

    end

    % Some extra blank lines to increase readability
    listStr = [listStr; {[strIndent '   |    ']}];

end % End iField for-loop


%% Field Filler
% To properly align the field names, a filler is required. To know how long
% the filler must be, the length of the longest fieldname must be found.
% Because 'fields' is a cell array, the function 'cellfun' can be used to
% extract the lengths of all fields.
maxFieldLength = max(cellfun(@length, fields));

%% Print non-structure fields without values
% Print non-structure fields without the values. This can be done very
% quick.
if printValues == 0
    
    noStrucFields = fields(isStruct == 0);

    for iField  = 1 : length(noStrucFields)

        Field = cell2mat(noStrucFields(iField));

        filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);

        listStr = [listStr; {[strIndent '   |' filler ' ' Field]}];

    end

    return

end


%% Select non-structure fields (to print with values)
% Select fields that are not a structure and group them by data type. The
% following groups are distinguished:
%   - characters and strings
%   - numeric arrays
%   - logical
%   - empty arrays
%   - matrices
%   - numeric scalars
%   - cell arrays
%   - other data types

% Character or string (array of characters)
isChar = structfun(@ischar, Structure);
charFields = fields(isChar == 1);

% Numeric fields
isNumeric = structfun(@isnumeric, Structure);

% Numeric scalars
isScalar = structfun(@isscalar, Structure);
isScalar = isScalar .* isNumeric;
scalarFields = fields(isScalar == 1);

% Numeric vectors (arrays)
isVector = structfun(@isvector, Structure);
isVector = isVector .* isNumeric .* not(isScalar);
vectorFields = fields(isVector == 1);

% Logical fields
isLogical = structfun(@islogical, Structure);
logicalFields = fields(isLogical == 1);

% Empty arrays
isEmpty = structfun(@isempty, Structure);
emptyFields = fields(isEmpty == 1);

% Numeric matrix with dimension size 2 or higher
isMatrix = structfun(@(x) ndims(x) >= 2, Structure);
isMatrix = isMatrix .* isNumeric .* not(isVector) ...
                    .* not(isScalar) .* not(isEmpty);
matrixFields = fields(isMatrix == 1);

% Cell array
isCell = structfun(@iscell, Structure);
cellFields = fields(isCell == 1);

% Datatypes that are not checked for
isOther = not(isChar + isNumeric + isCell + isStruct + isLogical + isEmpty);
otherFields = fields(isOther == 1);



%% Print non-structure fields
% Print all the selected non structure fields
% - Strings are printed to a certain amount of characters
% - Vectors are printed as long as they are shorter than maxArrayLength
% - Matrices are printed if they have less elements than maxArrayLength
% - The values of cells are not printed


% Start with printing strings and characters. To avoid the display screen 
% becoming a mess, the part of the string that is printed is limited to 31 
% characters. In the future this might become an optional parameter in this
% function, but for now, it is placed in the code itself.
% if the string is longer than 31 characters, only the first 31  characters
% are printed, plus three dots to denote that the string is longer than
% printed.

maxStrLength = 31;

for iField = 1 : length(charFields)

    Field = cell2mat(charFields(iField));

    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    if (size(Structure.(Field), 1) > 1) && (size(Structure.(Field), 2) > 1)
        
        varStr = createArraySize(Structure.(Field), 'char');
        
    elseif length(Field) > maxStrLength
        
        varStr = sprintf(' ''%s...''', Structure.(Field(1:maxStrLength)));
        
    else
        
        varStr = sprintf(' ''%s''', Structure.(Field));
        
    end

    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' :' varStr]}];
end


% Print empty fields

for iField = 1 : length(emptyFields)
    
    
    Field = cell2mat(emptyFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' : [ ]' ]}];

end


% Print logicals. If it is a scalar, print true/false, else print vector
% information

for iField = 1 : length(logicalFields)
    
    Field = cell2mat(logicalFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    if isscalar(Structure.(Field))
        
        logicalValue = {'False', 'True'};
        
        varStr = sprintf(' %s', logicalValue{Structure.(Field) + 1});

    else

        varStr = createArraySize(Structure.(Field), 'Logic array');
                
    end

    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' :' varStr]}];
    
end


% Print numeric scalar field. The %g format is used, so that integers,
% floats and exponential numbers are printed in their own format.

for iField = 1 : length(scalarFields)
    
    Field = cell2mat(scalarFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    varStr = sprintf(' %g', Structure.(Field));

    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' :' varStr]}];

end


% Print numeric array. If the length of the array is smaller then
% maxArrayLength, then the values are printed. Else, print the length of
% the array.

for iField = 1 : length(vectorFields)
    
    Field = cell2mat(vectorFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    if length(Structure.(Field)) > maxArrayLength
        
        varStr = createArraySize(Structure.(Field), 'Array');
        
    else

        varStr = sprintf('%g ', Structure.(Field));

        varStr = ['[' varStr(1:length(varStr) - 1) ']'];
                    
    end
    
    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' : ' varStr]}];

end


% Print numeric matrices. If the matrix is two-dimensional and has more
% than maxArrayLength elements, only its size is printed.
% If the matrix is 'small', the elements are printed in a matrix structure.
% The top and the bottom of the matrix is indicated by a horizontal line of
% dashes. The elements are also lined out by using a fixed format
% (%#10.2e). Because the name of the matrix is only printed on the first
% line, the space is occupied by this name must be filled up on the other
% lines. This is done by defining a 'filler2'.
% This method was developed by S. Wegerich.

for iField = 1 : length(matrixFields)
    
    Field = cell2mat(matrixFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    if numel(Structure.(Field)) > maxArrayLength
        
        varStr = createArraySize(Structure.(Field), 'Array');

        varCell = {[strIndent '   |' filler ' ' Field ' :' varStr]};
        
    else

        matrixSize = size(Structure.(Field));
        
        filler2 = char(ones(1, maxFieldLength + 6) * 32);

        dashes = char(ones(1, 12 * matrixSize(2) + 1) * 45);

        varCell = {[strIndent '   |' filler2 dashes]};
        
        % first line with field name
        varStr = sprintf('%#10.2e |', Structure.(Field)(1, :));

        varCell = [varCell; {[strIndent '   |' filler ' ' ...
                              Field ' : |' varStr]}];

        % second and higher number rows
        for j = 2 : matrixSize(1)

            varStr = sprintf('%#10.2e |', Structure.(Field)(j, :));
            
            varCell = [varCell; {[strIndent '   |' filler2 '|' varStr]}];
        end

        varCell = [varCell; {[strIndent '   |' filler2 dashes]}];
                    
    end
    
    listStr = [listStr; varCell];

end


% Print cell array information, i.e. the size of the cell array. The
% content of the cell array is not printed.

for iField = 1 : length(cellFields)

    Field = cell2mat(cellFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    varStr = createArraySize(Structure.(Field), 'Cell');

    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' :' varStr]}];

end


% Print unknown datatypes. These include objects and user-defined classes

for iField = 1 : length(otherFields)

    Field = cell2mat(otherFields(iField));
    
    filler = char(ones(1, maxFieldLength - length(Field) + 2) * 45);
    
    varStr = createArraySize(Structure.(Field), 'Unknown');

    listStr = [listStr; {[strIndent '   |' filler ' ' Field ' :' varStr]}];

end

end



%% FUNCTION: getIndentation
% This function creates the hierarchical indentations

function str = getIndentation(indent)
    x = '   |    ';
    str = '';
    
    for i = 1 : indent
        str = cat(2, str, x);
    end
end



%% FUNCTION: createArraySize
% This function returns a string with the array size of the input variable
% like: "[1x5 Array]" or "[2x3x5 Structure]" where 'Structure' and 'Array'
% are defined by the type parameter

function varStr = createArraySize(varName, type)
    varSize = size(varName);

    arraySizeStr = sprintf('%gx', varSize);
    arraySizeStr(length(arraySizeStr)) = [];
    
    varStr = [' [' arraySizeStr ' ' type ']'];
end

end
function x = easy_resample(x, s, r, tol)
%EASY_RESAMPLE  upsample or downsample a 1-d signal
%
% Provides an easy to use wrapper for interp and resample.
%
% USAGE:
% xr = easy_resample(x, s, r, [tol])
%
%INPUTS:
%
%        x: the signal (a 1-D vector)
%        s: sampling rate of the original signal (arbitrary units)
%        r: desired sampling rate of the resampled signal (arbitrary units)
%      tol: optional argument that controls precision of rounding (to
%           nearest tol units)
%
%OUTPUTS:
%
%       xr: the resampled signal
%
%EXAMPLE USAGE:
%       x = 1:10;
%       y = sin(x);
%       hold on;
%       plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
%       plot(easy_resample(x, 1, 10), easy_resample(y, 1, 10), 'k.-', 'LineWidth', 2);
%       plot(easy_resample(x, 1, 10), sin(easy_resample(x, 1, 10)), 'b--', 'LineWidth', 2);
%       hold off;
%       legend({'Original' 'Easy Resample' 'Ground Truth'});
%
%
%SEE ALSO: INTERP, RESAMPLE

%CHANGELOG
% 9-1-13   JRM   Wrote it.

assert(isscalar(s) && s > 0, 'sampling rate must be a positive scalar')
assert(isscalar(r) && r > 0, 'resampling rate must be a positive scalar')
assert(~exist('tol','var') || (isscalar(tol) && tol > 0 && tol < 1), 'tolerance must be a positive scalar between 0 and 1')

if ~exist('tol','var'), tol = 0.1; end

p = s/r;

if p < 1
    x = interp(x,1/p);
elseif p > 1       
    x = resample(x,1/tol,round(p*1/tol));
end
end
function ux = upsamplex(x, osrate)

    tmp = findseq(x);
    if sum(tmp(:,1))>0
        on = tmp(tmp(:,1)==1, 2);
        off = tmp(tmp(:,1)==1,3);
    else
        on = find(x); off = on; 
    end
    ux = zeros(length(x)*osrate, 1);
    on2 = fix(on*osrate); 
    off2 = fix(off*osrate);
    for i = 1:length(on2)
        ux(on2(i):off2(i)) = 1; 
    end

end
function [d, id] = getchunks(a, opt)

%GETCHUNKS Get the number of repetitions that occur in consecutive chunks.
%   C = GETCHUNKS(A) returns an array of n elements, where n is the number
%   of consecutive chunks (2 or more repetitions) in A, and each element is
%   the number of repetitions in each chunk. A can be LOGICAL, any
%   numeric vector, or CELL array of strings. It can also be a character
%   array (see below, for its special treatment).
%
%   [C, I] = GETCHUNKS(A) also returns the indices of the beginnings of the
%   chunks.
%
%   If A is a character array, then it finds words (consecutive
%   non-spaces), returning the number of chararcters in each word and the
%   indices to the beginnings of the words.
%
%   GETCHUNKS(A, OPT) accepts an optional argument OPT, which can be any of
%   the following three:
%
%       '-reps'  : return repeating chunks only. (default)
%       '-full'  : return chunks including single-element chunks.
%       '-alpha' : (for CHAR arrays) only consider alphabets and numbers as
%                  part of words. Punctuations and symbols are regarded as
%                  spaces.
%
%   Examples:
%     A = [1 2 2 3 4 4 4 5 6 7 8 8 8 8 9];
%     getchunks(A)
%       ans =
%           2   3   4
%
%
%     B = 'This is a generic (simple) sentence';
%     [C, I] = getchunks(B)
%       C =
%            4     2     1     7     8     8
%       I =
%            1     6     9    11    19    28
%
%
%     [C, I] = getchunks(B, '-alpha')
%       C =
%            4     2     1     7     6     8
%       I =
%            1     6     9    11    20    28
%
%   See also HIST, HISTC.
%
%   VERSIONS:
%     v1.0 - first version
%     v1.0.1 - added option '-alpha'
%

% Copyright 2009 The MathWorks, Inc.

%--------------------------------------------------------------------------
% Error checking
%--------------------------------------------------------------------------
error(nargchk(1, 2, nargin));
if ndims(a) > 2 || min(size(a)) > 1
  error('Input must be a 2-D vector');
end

alphanumeric = false;
fullList     = false;

%--------------------------------------------------------------------------
% Process options
%--------------------------------------------------------------------------
if nargin == 2
  if ~ischar(opt)
    error('Additional argument must be a string array');
  end
  
  % Allow for partial arguments
  possibleOptions = ['-full '; '-reps '; '-alpha'];
  iOpt = strmatch(lower(opt), possibleOptions);
  
  if isempty(iOpt) || length(iOpt) > 1
    error('Invalid option. Allowed option: ''-full'', ''-reps'', ''-alpha''');
  else
    switch iOpt
      
      case 1  % '-full'
        % Include single-element chunks
        fullList = true;
        if ischar(a)
          fprintf('''-full'' option not applicable to CHAR arrays.\n');
        end
        
      case 2  % '-reps'
        % Only find 2 or more repeating blocks
        fullList = false;
        
      case 3  % '-alpha'
        % For char arrays, only consider alphabets and numbers as part of
        % words. Punctuations and symbols are regarded as space.
        alphanumeric = true;
        if ~ischar(a)
          fprintf('''-alpha'' option only applicable to CHAR arrays.\n');
        end
        
    end
  end
end

%--------------------------------------------------------------------------
% Convert to a row vector for STRFIND
%--------------------------------------------------------------------------
a = a(:)';

%--------------------------------------------------------------------------
% Deal with differet classes
%--------------------------------------------------------------------------
switch class(a)
  
  case 'double'
    % Leave as is
    
  case {'logical', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single'}
    % Convert to DOUBLE
    a = double(a);
    
  case 'char'
    if alphanumeric % Get alphabet and number locations
      try % call C-helper function directly (slightly faster)
        a = isletter(a) | ismembc(a, 48:57);
      catch %#ok<CTCH>
        a = isletter(a) | ismember(a, 48:57);
      end
      
    else  % Get non-space locations
      a = ~isspace(a);  
    end
  
  case 'cell'
    % Convert cell array of strings into unique numbers
    if all(cellfun('isclass', a, 'char'))
      [tmp, tmp, a] = unique(a); %#ok<ASGLU>
    else
      error('Cell arrays must be array of strings.');
    end
    
  otherwise
    error('Invalid type. Allowed type: CHAR, LOGICAL, NUMERIC, and CELL arrays of strings.');
end

%--------------------------------------------------------------------------
% Character arrays (now LOGICAL) are dealt differently
%--------------------------------------------------------------------------
if islogical(a)
  % Pad the array
  a  = [false, a, false];

  % Here's a very convoluted engine
  b  = diff(a);
  id = strfind(b, 1);
  d  = strfind(b, -1) - id;

%--------------------------------------------------------------------------
% Everything else (numeric arrays) are processed here
else
  % Pad the array
  a                 = [NaN, a, NaN];

  % Here's more convoluted code
  b                 = diff(a);
  b1                = b;  % to be used in fullList (below)
  ii                = true(size(b));
  ii(strfind(b, 0)) = false;
  b(ii)             = 1;
  c                 = diff(b);
  id                = strfind(c, -1);
  
  % Get single-element chunks also
  if fullList
  
    % And more convoluted code
    b1(id)          = 0;
    ii2             = find(b1(1:end-1));
    d               = [strfind(c, 1) - id + 1, ones(1, length(ii2))];
    id              = [id,ii2];
    [id,tmp]        = sort(id);
    d               = d(tmp);
    
  else
    
    d               = strfind(c, 1) - id + 1;
    
  end
end
end
function varargout = findseq(A,dim)

% FINDSEQ Find sequences of repeated (adjacent/consecutive) numeric values
%
%   FINDSEQ(A) Find sequences of repeated numeric values in A along the
%              first non-singleton dimension. A should be numeric.
%
%   FINDSEQ(...,DIM) Look for sequences along the dimension specified by the 
%                    positive integer scalar DIM.
%
%   OUT = findseq(...)
%       OUT is a "m by 4" numeric matrix where m is the number of sequences found.
%       
%       Each sequence has 4 columns where:
%           - 1st col.:  the value being repeated
%           - 2nd col.:  the position of the first value of the sequence
%           - 3rd col.:  the position of the last value of the sequence
%           - 4th col.:  the length of the sequence
%       
%   [VALUES, INPOS, FIPOS, LEN] = findseq(...)
%       Get OUT as separate outputs. 
%
%       If no sequences are found no value is returned.
%       To convert positions into subs/coordinates use IND2SUB
%
% 
% Examples:
%
%     % There are sequences of 20s, 1s and NaNs (column-wise)
%     A   =  [  20,  19,   3,   2, NaN, NaN
%               20,  23,   1,   1,   1, NaN
%               20,   7,   7, NaN,   1, NaN]
%
%     OUT = findseq(A)
%     OUT =  
%            20        1          3        3
%             1       14         15        2
%           NaN       16         18        3
%     
%     % 3D sequences: NaN, 6 and 0
%     A        = [  1, 4
%                 NaN, 5
%                   3, 6];
%     A(:,:,2) = [  0, 0
%                 NaN, 0
%                   0, 6];
%     A(:,:,3) = [  1, 0
%                   2, 5
%                   3, 6];
%     
%     OUT = findseq(A,3)
%     OUT = 
%             6     6    18     3
%             0    10    16     2
%           NaN     2     8     2
%
% Additional features:
% - <a href="matlab: web('http://www.mathworks.com/matlabcentral/fileexchange/28113','-browser')">FEX findseq page</a>
% - <a href="matlab: web('http://www.mathworks.com/matlabcentral/fileexchange/6436','-browser')">FEX rude by us page</a>
%
% See also: DIFF, FIND, SUB2IND, IND2SUB

% Author: Oleg Komarov (oleg.komarov@hotmail.it) 
% Tested on R14SP3 (7.1) and on R2012a. In-between compatibility is assumed.
% 02 jul 2010 - Created
% 05 jul 2010 - Reorganized code and fixed bug when concatenating results
% 12 jul 2010 - Per Xiaohu's suggestion fixed bug in output dimensions when A is row vector
% 26 aug 2010 - Cast double on logical instead of single
% 28 aug 2010 - Per Zachary Danziger's suggestion reorganized check structure to avoid bug when concatenating results
% 22 mar 2012 - Per Herbert Gsenger's suggestion fixed bug in matching initial and final positions; minor change to distribution of OUT if multiple outputs; added 3D example 
% 08 nov 2013 - Fixed major bug in the sorting of Final position that relied on regularity conditions not always verified

% NINPUTS
error(nargchk(1,2,nargin));

% NOUTPUTS
error(nargoutchk(0,4,nargout));

% IN
if ~isnumeric(A)
    error('findseq:fmtA', 'A should be numeric')
elseif isempty(A) || isscalar(A)
    varargout{1} = [];
    return
elseif islogical(A)
    A = double(A);
end

% DIM
szA = size(A);
if nargin == 1 || isempty(dim)
    % First non singleton dimension
    dim = find(szA ~= 1,1,'first');
elseif ~(isnumeric(dim) && dim > 0 && rem(dim,1) == 0) || dim > numel(szA)
    error('findseq:fmtDim', 'DIM should be a scalar positive integer <= ndims(A)');
end

% Less than two elements along DIM
if szA(dim) == 1
    varargout{1} = [];
    return
end

% ISVECTOR
if nnz(szA ~= 1) == 1
    A = A(:);
    dim = 1;
    szA = size(A);
end

% Detect 0, NaN, Inf and -Inf
OtherValues    = cell(1,4);
OtherValues{1} = A ==    0;
OtherValues{2} = isnan(A) ;
OtherValues{3} = A ==  Inf;
OtherValues{4} = A == -Inf;
Values         = [0,NaN, Inf,-Inf];

% Remove zeros
A(OtherValues{1}) = NaN;                             

% Make the bread
bread = NaN([szA(1:dim-1),1,szA(dim+1:end)]);

% [1] Get chunks of "normal" values
Out = mainengine(A,bread,dim,szA);

% [2] Get chunks of 0, NaN, Inf and -Inf
for c = 1:4
    if nnz(OtherValues{c}) > 1
        % Logical to double and NaN padding
        OtherValues{c} = double(OtherValues{c});                        
        OtherValues{c}(~OtherValues{c}) = NaN;                          
        % Call mainengine and concatenate results
        tmp = mainengine(OtherValues{c}, bread,dim,szA);
        if ~isempty(tmp)
            Out = [Out; [repmat(Values(c),size(tmp,1),1) tmp(:,2:end)]];  %#ok
        end
    end
end

% Distribute output
if nargout < 2 
    varargout = {Out};
else
    varargout = num2cell(Out(:,1:nargout),1);
end

end
% MAINENGINE This functions uses run length encoding and retrieve positions 
function Out = mainengine(meat,bread,dim,szMeat)

% Make a sandwich  
sandwich    = cat(dim, bread, meat, bread);

% Find chunks (run length encoding engine)
IDX         = diff(diff(sandwich,[],dim) == 0,[],dim);

% Initial and final row/col subscripts
[rIn, cIn]  = find(IDX  ==  1);
[rFi, cFi]  = find(IDX  == -1);

% Make sure row/col subs correspond (relevant if dim > 1)
[In, idx]   = sortrows([rIn, cIn],1);
Fi          = [rFi, cFi];
Fi          = Fi(idx,:);

% Calculate length of blocks
if dim < 3
    Le = Fi(:,dim) - In(:,dim) + 1;
else
    md = prod(szMeat(2:dim-1));
    Le = (Fi(:,2) - In(:,2))/md + 1;
end

% Convert to linear index
InPos       = sub2ind(szMeat,In(:,1),In(:,2));
FiPos       = sub2ind(szMeat,Fi(:,1),Fi(:,2));

% Assign output
Out         = [meat(InPos),...    % Values
               InPos      ,...    % Initial positions 
               FiPos      ,...    % Final   positions
               Le         ];      % Length of the blocks
end
function tilefigs(handles,resize,nRows,nCols,leftRightSpacing,topBottomSpacing,...
    border,monitor,monitorLocation,monitorSize)
% TILEFIGS Tile figures (spread them out). Tile figs has a number of
% arguments, which may either be specified, or entered as [] to use the
% default
%   TILEFIGS(handles,resize,nRows,nCols,leftRightSpacing,...
%       topBottomSpacing,border,monitor,monitorLocation,monitorSize)
%   TILEFIGS , by itself, resizes and tiles all open figures on to the primary monitor
%   TILEFIGS(handles) Resizes and tiles the figures specified in handles,
%   with the first handle in the upper left corner
%   TILEFIGS(...,resize) If resize is set to false, the figures will be
%   moved but not resized. Effort has been made to have minimal overlap,
%   but some may still occur.
%   TILEFIGS(...,nRows) Specifies a number of rows for the tile grid
%   TILEFIGS(...,nCols) Specifies a number of column for the tile grid.
%   Note that the figures are still included row by row, so if one has 8
%   open figures and tilefigs([],[],[],5) is called, the figures will be
%   put into a 4 x 2 grid with space left for a 5th column.
%   TILEFIGS(...,leftRightSpacing) leaves a horizontal space between the
%   figures. Spacing is specified in pixels. This may have no effect if
%   resize is set to false.
%   TILEFIGS(...,topBottomSpacing) Same as above, except for vertical
%   spacing
%   TILEFIGS(...,border) Leaves a border around the tiling of figures.
%   border is a 1 x 4 matrix, which is [leftSpace bottomSpace rightSpace
%   topSpace]. If the primary monitor is used, the default is [0 35 0 0] 
%   which is approximately the size of the taskbar in Windows 7. If a
%   secondary monitor is used, the default is [0 0 0 0]
%   TILEFIGS(...,monitor , monitorLocation,monitorSize) Specifies a monitor to use, 
%   which is the row of the call get(0,'MonitorPositions'). MATLAB does not
%   officially support dual monitors. The author has found that the
%   MonitorPositions call correctly gets the size of the monitor (in
%   pixels), but not the location of the monitor. monitorLocation and
%   monitorSize overrides the location and size of the monitor.

% Written by Brendan Tracey October 11th, 2012.

% Inspirational credit to:
%Charles Plum                    Nichols Research Corp.
%<cplum@nichols.com>             70 Westview Street
%Tel: (781) 862-9400             Kilnbrook IV
%Fax: (781) 862-9485             Lexington, MA 02173

% Nomenclature -- first number is width, second is length
% x is left to right, y is up and down

if ~exist('topBottomSpacing','var') || isempty(topBottomSpacing)
    topBottomSpacing = 0;
end
if ~exist('leftRightSpacing','var') || isempty(leftRightSpacing)
    leftRightSpacing = 0;
end

%% Select the monitor and get the location and size of the monitor

monitorPositions = get(0,'MonitorPositions');
matlabMonitorLocation = monitorPositions(:,1:2);
matlabMonitorSize = monitorPositions(:,3:4) - monitorPositions(:,1:2) + 1;

if ~exist('monitor','var') || isempty(monitor)
    % No monitor specified, use the primary
    monitor = 1;
end

if ~exist('monitorSize','var') || isempty(monitorSize)
    % No monitor size specified, use what Matlab says
    monitorSize = matlabMonitorSize(monitor,:);
end

if ~exist('monitorLocation','var') || isempty(monitorLocation)
    % No monitor location specified, use what Matlab says
    monitorLocation = matlabMonitorLocation(monitor,1:2);
end

% Set the default border if none is specified [left, bottom, right, top]
if monitor > 1
    if ~exist('border','var') || isempty(border)
        border = [0 0 0 0];
    end
else
    if ~exist('border','var') || isempty(border)
        border = [0 30 0 0];
    end
end
% Modify the monitor size and location to account for the border
monitorLocation = monitorLocation + border(1:2);
monitorSize(1) = monitorSize(1) - border(1) - border(3);
monitorSize(2) = monitorSize(2) - border(2) - border(4);

%% Select figures to use
if ~exist('handles','var') || isempty(handles)
    % No figure handles specified, select all figures
    handles = get (0,'Children'); % locate all open figure handles
    handles = handles(end:-1:1); % Re-order so that first created is in upper left
end
nFigures = length(handles);

% Put all of the figures in pixel units
units = get(handles, 'units');
set(handles,'units','pixels');

%% Determine the grid for figures
if (~exist('nRows','var') || isempty(nRows)) && (~exist('nCols','var') || isempty(nCols))
    % No grid specified, choose the grid which roughly matches the monitor
    % aspect ratio
    monitorAspectRatio = monitorSize(1) / monitorSize(2);
    if monitorAspectRatio < 1
       nCols = round(sqrt(nFigures));
       nRows = ceil(nFigures/nCols);
    else
       nRows = round(sqrt(nFigures));
       nCols = ceil(nFigures/nRows);
    end
elseif (exist('nRows','var') && ~isempty(nRows)) && (~exist('nCols','var') || isempty(nCols))
    nCols = ceil(nFigures/nRows);
elseif (~exist('nRows','var') || isempty(nRows)) && (exist('nCols','var') && ~isempty(nCols))
    nRows = ceil(nFigures/nCols);
elseif (exist('nRows','var') && ~isempty(nRows)) && (exist('nCols','var') && ~isempty(nCols))
    if nRows*nCols < nFigures
        error('Grid size not big enough')
    end
else
    error('Should not be here')
end

%% Calculate the grid sizing
if ~exist('resize','var') || isempty(resize)
    % Resize not set, default is to resize
    resize = 1;
end
if resize
    width = (monitorSize(1) - leftRightSpacing*(nCols - 1))/nCols;
    height = (monitorSize(2) - topBottomSpacing*(nRows - 1))/nRows;
else
    % Make the spacing equal with the constraint that the figures do not go
    % off the edge of the screen.
    figureWidths = zeros(nFigures,1);
    figureHeights = zeros(nFigures,1);
    for ii = 1:nFigures
        figSize = get(handles(ii),'OuterPosition');
        figureWidths(ii) = figSize(3);
        figureHeights(ii) = figSize(4);
    end
    widthEndInds = nCols:nCols:nFigures;
    heightEndInds = 1:nCols;
    maxWidth = max(figureWidths(widthEndInds));
    maxHeight = max(figureHeights(heightEndInds));
    width = (monitorSize(1) - maxWidth)/max(nCols-1,1);
    if nRows ==1
        height = 0;
    else
        height = (monitorSize(2) - maxHeight)/(nRows-1);
    end
end

%% Move and resize the figures
if resize
    pnum = 0;
    for row = 1:nRows
        for col = 1:nCols
            pnum = pnum+1;
            if (pnum>nFigures)
                break
            end
            xLocation = monitorLocation(1) + (col - 1)*width + (col - 1)*leftRightSpacing;
            yLocation = monitorLocation(2) + monitorSize(2) - row*height ...
                - (row - 1)*topBottomSpacing;
            figure(handles(pnum))
            set(handles(pnum),'OuterPosition',[xLocation yLocation width height]);
        end
    end
else
    % Not resizing, set initial locations
    xLocations = inf*ones(nRows,nCols);
    yLocations = inf*ones(nRows,nCols);
    widthMat = zeros(nRows,nCols);
    heightMat = zeros(nRows,nCols);
    pnum = 0;
    for row = 1:nRows
        for col = 1:nCols
            pnum = pnum+1;
            if (pnum>nFigures)
                break
            end
            widthMat(row,col) = figureWidths(pnum);
            heightMat(row,col) = figureHeights(pnum);
        end
    end
    pnum = 0;
    for row = 1:nRows
        for col = 1:nCols
            pnum = pnum+1;
            if (pnum>nFigures)
                break
            end
            xLocation = (col - 1)*width;
            yLocation = monitorSize(2) - (row - 1)*height - heightMat(row,col);
            xLocations(row,col) = xLocation;
            yLocations(row,col) = yLocation;
        end
    end
    
    % Modify the positions to make it look nicer.
    % First, if the figures are too big to not overlap, spread the rows and
    % columns out to use up all of the space at the border to minimize
    % overlap. The subtraction of the width if there are no figures is to 
    % avoid the issue of both the column and the row being moved into the
    % empty spot
    if (sum(sum(widthMat,2) > monitorSize(1)) > 0)
        for row = 1:nRows
            lastNonemptyCol = find(xLocations(row,:) < inf,1,'Last');
            whitespace = monitorSize(1) - xLocations(row,lastNonemptyCol) - widthMat(row,lastNonemptyCol) - width*(nCols-lastNonemptyCol);
            if lastNonemptyCol>1
                for col = 1:lastNonemptyCol
                    xLocations(row,col) = xLocations(row,col) + (col-1)/(lastNonemptyCol-1) * whitespace;
                end
            end
        end
    end
    if sum(sum(heightMat,1) > monitorSize(2)) > 0
        for col = 1:nCols
            lastNonemptyRow = find(xLocations(:,col) < inf,1,'Last');
            whitespace = yLocations(lastNonemptyRow,col) - 1 - height*(nRows - lastNonemptyRow);
            if lastNonemptyRow>1
                for row = 1:lastNonemptyRow
                    yLocations(row,col) = yLocations(row,col) - ((row-1)/(nRows-1)) * whitespace;
                end
            end
        end
    end
    % If the figures can be condensed, slide the rows and columns in.
    for col = 1:nCols - 1
        lastNonemptyRow = find(xLocations(:,col) < inf,1,'Last');
        endOfPlot = xLocations(1:lastNonemptyRow,col) + widthMat(1:lastNonemptyRow,col) + leftRightSpacing;
        beginningOfNextPlot = xLocations(1:lastNonemptyRow,col + 1);
        gap = beginningOfNextPlot - endOfPlot;
        minGap = max(min(gap),0);
        xLocations(1:lastNonemptyRow,col+1) = xLocations(1:lastNonemptyRow,col+1) - minGap;
    end
    for row = 2:nRows
        lastNonemptyCol = find(xLocations(row,:) < inf,1,'Last');
        bottomOfUpperRow = yLocations(row-1,1:lastNonemptyCol) - topBottomSpacing;
        topOfThisRow = yLocations(row,1:lastNonemptyCol) + heightMat(row,1:lastNonemptyCol);
        gap =  bottomOfUpperRow - topOfThisRow;
        minGap = max(min(gap),0);
        yLocations(row,1:lastNonemptyCol) = yLocations(row,1:lastNonemptyCol) + minGap;
    end
    % Finally, actually move and replot the figures.
    pnum = 0;
    for row = 1:nRows
        for col = 1:nCols
            pnum = pnum+1;
            if (pnum>nFigures)
                break
            end
            figure(handles(pnum))
            xloc = xLocations(row,col);
            yloc = yLocations(row,col);
            set(handles(pnum),'OuterPosition',[monitorLocation(1)+xloc monitorLocation(2) + yloc figureWidths(pnum) figureHeights(pnum)]);
        end
    end
end
for ii = 1:nFigures
    set(handles(ii),'units',units{ii});
end
end
function BeamProfileFits_v4_Weighted(varargin)
% Beam Profile Fitting routine
%
% INPUTS: 
%   None
% 
% OUTPUTS: 
%   Plots showing fit parameters and results as well as fit object for
%   vertical and horizontal fits.
%
% NOTE:
%   - Be careful fitting M^2 if you have not profiled the waist. If unsure default to fit_M2 = 0;
%   - Distances are expected in centimeters 
%   - Beam diameter is expected in microns
%   - Data format should order columns as 
%       Col. 1        | Col. 2         | Col. 3         | Col. 4         | Col. 5         |
%       Distance (cm) | Vert spot (um) | Vert Unc. (um) | Horz spot (um) | Horz Unc. (um) |
%
% Fitting Model
%   fittype('sqrt(a^2+((lambda/pi)^2/a^2)*(x-c)^2)' from Kogelnik and Li
% Flags below suppress certain warnings that Matlab was complainging about
%#ok<*UNRCH>
%#ok<*NASGU>
%#ok<*ASGLU>
% CHANGELOG
% v5 - added variable input arguments so you can pass the lambda and
% fitfile in if you specify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variables
% For saving you can enter the direct path to your dataset here. If a path
% is not present then one will be asked for.
lambda  = '';
dir     = './'; % Full path to folder containing fitfile
fitfile = ''; % Name of fitfile containing data

% Decide whether to fix M2 value or allow to vary
fitM2 = 0; %1 for true, 0 for false

% Decide whether to use weights or not
weightFit = 0; %1 for true, 0 for false

% Default to not passing data through varargin, but instead rad from file
data_passed_in = 0;

% Figure properties
factorXBounds = 3; %Factor to zoom out X limits for the propagation figure
%% Allow runtime determination of some behavior
if ~isempty(varargin)
    for varArgInd = 1:2:length(varargin)
        varName = varargin{varArgInd};
        switch varName
            case 'lambda'
                lambda = varargin{varArgInd + 1};
            case 'abs_filepath'
                abs_filepath = varargin{varArgInd + 1};
            case 'fit_M2'
                abs_filepath = varargin{varArgInd + 1};
            case 'weightFit'
                abs_filepath = varargin{varArgInd + 1};
            case 'factorXBounds'
                abs_filepath = varargin{varArgInd + 1};
            case 'rawData'
                data_passed_in = 1;
                input_table = varargin{varArgInd + 1};
                rawData = cell(1,width(input_table));
                for i = 1:width(input_table); 
                    rawData{i} = input_table{:,i}; 
                end
        end
    end
end

%% Get variable input from user
while isempty(lambda)
    % Ask for raw data
    lambda = input('Enter wavelength (nm): ')*10^-7; %wavelength in cm
    if isempty(lambda) || isnan(lambda)
        error('Invalid selection. Please run again.');
    else
        break
    end
end

if (isempty(fitfile) || isempty(dir)) && ~data_passed_in
    % Ask for raw data
    [fitfile,dir] = uigetfile('*.txt','Choose the fit file');
    abs_filepath = [dir fitfile];
end

%% Define function handles
% Fitting function (gaussian beam propagation)
rayLength = @(w0) pi*w0^2/lambda;
if fitM2;
    waistFunc = @(w0,x0,M,x) w0*1e-4*M.*sqrt(1 + ((x - x0)./rayLength(w0*1e-4)).^2); 
else
    waistFunc = @(w0,x0,x)   w0*1e-4.*sqrt(1 + ((x - x0)./rayLength(w0*1e-4)).^2);
end

% Normalize and scale weights
wghtScaleFunc = @(wght) 1/max(wght).*wght;


%% Load and process raw data
% Reserve figures
figFit  = figure;

% Get raw data in matlab
if ~data_passed_in
    dataFileID = fopen(abs_filepath);
    rawData    = textscan(dataFileID,'%f%f%f%f%f','commentstyle','%'); fclose(dataFileID);
end

% Define Limits used for fitting and create function handle for convenience
if fitM2
    lowBnd = [0,-Inf,1];
    uppBnd = [Inf,Inf,Inf];
else
    lowBnd = [0,-Inf];
    uppBnd = [Inf,Inf];
end
fitFuncOpt = @(beta,wghtVec) ...
    fitoptions('Method'     ,   'NonlinearLeastSquares' ,...
               'TolFun'     ,   1e-10                   ,...
               'TolX'       ,   1e-10                   ,...
               'Weights'    ,   wghtVec                 ,...
               'Lower'      ,   lowBnd                  ,...
               'Upper'      ,   uppBnd                  ,...
               'StartPoint' ,   beta                    );

% Convert to Units for fitting
D     = rawData{1};     % Expect distance in centimeters
Vspot = rawData{2}./2;  % divide by 2 for diameter to radius (gaussian fit)
Vstd  = rawData{3}./2; 
Wspot = rawData{4}./2;
Wstd  = rawData{5}./2;

%% Initial guesses for fitting
[minVal,minPos]  = min(Vspot);
guessPosVert     = D(minPos); %vertical waist position best guess for fit (cm)
guessWaistVert   = minVal;    %vertical waist size best guess for fit (um)
guessM2Vert      = 1;

[minVal,minPos]  = min(Wspot);
guessPosHoriz    = D(minPos);   %horizontal waist position best guess for fit (cm)
guessWaistHoriz  = minVal;      %horizontal waist size best guess for fit (um)
guessM2Horiz     = 1;

%% --- Create fit - Vertical Fit
% Build correct weight vector and initial guess vector
beta = [guessWaistVert guessPosVert];
if fitM2; beta(end + 1) = guessM2Vert; end
if weightFit; wghtVec = 1./Vstd.^2; else wghtVec = ones(1,length(D)); end
funcFitType = fittype(waistFunc,'options',fitFuncOpt(beta,wghtScaleFunc(wghtVec)));

% Fit horizontal model using new data
[cfun_Vert,gof_Vert,output_Vert] = fit(D,Vspot*1e-4,funcFitType); 
cf_Vert                          = coeffvalues(cfun_Vert);


%% --- Create fit - Horizontal Fit
% Build correct weight vector and initial guess vector
beta = [guessWaistHoriz guessPosHoriz];
if fitM2; beta(end + 1) = guessM2Horiz; end
if weightFit; wghtVec = 1./Wstd.^2; else wghtVec = ones(1,length(D)); end
funcFitType = fittype(waistFunc,'options',fitFuncOpt(beta,wghtScaleFunc(wghtVec)));

% Fit horizontal model using new data
[cfun_Horz,gof_Horz,output_Horz] =  fit(D,Wspot*1e-4,funcFitType);
cf_Horz                          = coeffvalues(cfun_Horz);

%% Calculate uncertainties
unc_vert     = confint(cfun_Vert);
vert_was_unc = unc_vert(:,1); vert_was_unc = diff(vert_was_unc)/2;
vert_pos_unc = unc_vert(:,2); vert_pos_unc = diff(vert_pos_unc)/2;
if fitM2; vert_M2_unc  = unc_vert(:,3); vert_M2_unc  = diff(vert_M2_unc)/2; end

unc_horz     = confint(cfun_Horz);
horz_was_unc = unc_horz(:,1); horz_was_unc = diff(horz_was_unc)/2;
horz_pos_unc = unc_horz(:,2); horz_pos_unc = diff(horz_pos_unc)/2;
if fitM2; horz_M2_unc  = unc_horz(:,3); horz_M2_unc  = diff(horz_M2_unc)/2; end

%% Calculate the predicted fits and bounds
% create x vector to plots fit points with
xFitLim  = [min(D) max(D)] + [-1 1] * 1.5 * diff([min(D) max(D)]);
x_pnts = linspace(xFitLim(1),xFitLim(end),1e3);

[ci_Vert,pred_Vert] = predint(cfun_Vert,x_pnts,0.95,'observation','on');
[ci_Horz,pred_Horz] = predint(cfun_Horz,x_pnts,0.95,'observation','on');

%% Plot close up of fit and raw data
% This function was automatically generated on 12-Nov-2004 11:15:57
% and subsequently modified... numerours times.

% Set up figure to receive datasets and fits
figure(figFit);

% --- Plot data originally in dataset "vwaist vs. D"
ax_ = subaxis(2,2,1,1,'SH',0.01,'ML',0.05,'MT',0.05); 
set(ax_,'Box','on','Linewidth',2,'FontSize',15);
axes(ax_); hold on; grid on;
if weightFit
    h_ = errorbar(D,Vspot,Vstd);
else
    h_ = plot(D,Vspot);
end
set(h_,'Parent',ax_,'Color','k',...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','o', 'MarkerSize',12,...
    'MarkerFaceColor',[0.847 0.161 0]);
ylabel('Vertical Waist Radius','FontSize',20)

% Plot the Vertical fit
figure(figFit)
ax_ = subaxis(2,2,1,1,'SH',0.01,'ML',0.05,'MT',0.05); 
hold on
h_ = plot(x_pnts,pred_Vert*1e4,x_pnts,ci_Vert*1e4);
set(h_(1),...
     'LineStyle','--', 'LineWidth',3,...
     'Marker','none', 'MarkerSize',6,...
     'Color',[0.502 0.502 0.502]);
set(h_(2:3),'Color',[0.502 0.502 0.502],...
    'LineStyle',':', 'LineWidth',2); 
 
% Nudge axis limits beyond data limits
xlim_ = [min(D) max(D)];
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end

% --- Plot data originally in dataset "wwaist vs. D"
ax_ = subaxis(2,2,1,2,'SH',0.01,'ML',0.05,'MT',0.05); 
set(ax_,'Box','on','Linewidth',2,'FontSize',15);
axes(ax_); hold on; grid on
if weightFit
    h_  = errorbar(D,Wspot,Wstd);
else
    h_  = plot(D,Wspot);
end
set(h_,'Parent',ax_,'Color','k',...
         'LineStyle','none', 'LineWidth',1,...
         'Marker','o', 'MarkerSize',12,...
         'MarkerFaceColor',[0.847 0.161 0]);
ylabel('Horiztonal Waist Radius','FontSize',20)

% Plot this Horizontal fit
figure(figFit)
ax_ = subaxis(2,2,1,2,'SH',0.01,'ML',0.05,'MT',0.05); 
hold on
h_ = plot(x_pnts,pred_Horz*1e4,x_pnts,ci_Horz*1e4);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0 0 1],...
     'LineStyle','--', 'LineWidth',3,...
     'Marker','none', 'MarkerSize',6,...
     'Color',[0.502 0.502 0.502]);
set(h_(2:3),'Color',[0.502 0.502 0.502],...
    'LineStyle',':', 'LineWidth',2);
 
 % Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end
 
% Add textbox showing the fit information 
% Table information
horz_pos_str = sprintf('Horizontal Waist Position (cm): %g +/- %g',cf_Horz(2),horz_pos_unc);
horz_was_str = sprintf('Horizontal Waist Size     (um): %g +/- %g',cf_Horz(1),horz_was_unc);
vert_pos_str = sprintf('Vertical Waist Position   (cm): %g +/- %g',cf_Vert(2),vert_pos_unc);
vert_was_str = sprintf('Vertical Waist Size       (um): %g +/- %g',cf_Vert(1),vert_was_unc);

if fitM2; 
    horz_M2_str  = sprintf('Horizontal M^2: %g +/- %g',cf_Horz(3)^2,2*horz_M2_unc);
    vert_M2_str  = sprintf('Vertical M^2: %g +/- %g',cf_Vert(3)^2,2*vert_M2_unc);
else
    horz_M2_str  = sprintf('Horizontal M^2 fixed to 1');
    vert_M2_str  = sprintf('Vertical M^2 fixed to 1');
end

data = {vert_pos_str; vert_was_str; vert_M2_str; ''; horz_pos_str; horz_was_str; horz_M2_str};

% Add legend and information to fit plot
hold on;
%legend(ax_,legh_, legt_,4,'Location','Best');
annotation('textbox'         ,   [.15,.8,.1,.1]  ,...
           'String'          ,   data            ,...
           'LineStyle'       ,   '-'             ,...
           'BackgroundColor' ,   [1 1 .84]       ,...
           'EdgeColor'       ,   [0 0 0]         ,...
           'FontSize'        ,   16              );
    
xlabel('Distance From Source (cm)','FontSize',20)

%% Show Theoretical beam propagation and data on a large scale
figure(figFit)
ax_ = subaxis(2,2,2,1,1,2,'MR',0.01,'MT',0.05); 
set(ax_,'Box','on','Linewidth',2,'FontSize',15);
hold on;
if weightFit
    errorbar(D,Wspot,Wstd,'b*'); hold on
    errorbar(D,Vspot,Vstd,'r+');
else
    plot(D,Wspot,'b*'); hold on
    plot(D,Vspot,'r+');
end
% Extend plot beyond set propagation bounds
xPredPnt = factorXBounds*xFitLim + [-1 1] * 0.5 * diff(factorXBounds*xFitLim);
xPred    = linspace(xPredPnt(1),xPredPnt(2),5e3); %cm
horzTheo = cfun_Horz(xPred)*1e4;
vertTheo = cfun_Vert(xPred)*1e4;

plot(xPred,horzTheo,'-b',xPred,vertTheo,'-r');
grid on;
xlim(xlim_ + [-1 1]*factorXBounds/2*diff(xlim_));
xlabel('Distance [cm]','Fontsize',20);
ylabel('Spot Size [um]','Fontsize',20);
legend('Horizontal (W)','Vertical (V)','Location','Best');

%% Clean up workspace but keep fit objects
clearvars -except cfun_Vert cfun_Horz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Separate sub-functions

function h=subaxis(varargin)
%SUBAXIS Create axes in tiled positions. (just like subplot)
%   Usage:
%      h=subaxis(rows,cols,cellno[,settings])
%      h=subaxis(rows,cols,cellx,celly[,settings])
%      h=subaxis(rows,cols,cellx,celly,spanx,spany[,settings])
%
% SETTINGS: Spacing,SpacingHoriz,SpacingVert
%           Padding,PaddingRight,PaddingLeft,PaddingTop,PaddingBottom
%           Margin,MarginRight,MarginLeft,MarginTop,MarginBottom
%           Holdaxis
%
%           all units are relative (i.e. from 0 to 1)
%
%           Abbreviations of parameters can be used.. (Eg MR instead of MarginRight)
%           (holdaxis means that it wont delete any axes below.)
%
% Example:
%
%   >> subaxis(2,1,1,'SpacingVert',0,'MR',0); 
%   >> imagesc(magic(3))
%   >> subaxis(2,'p',.02);
%   >> imagesc(magic(4))
%
% 2001-2014 / Aslak Grinsted  (Feel free to modify this code.)

f=gcf;

UserDataArgsOK=0;
Args=get(f,'UserData');
if isstruct(Args) 
    UserDataArgsOK=isfield(Args,'SpacingHorizontal')&isfield(Args,'Holdaxis')&isfield(Args,'rows')&isfield(Args,'cols');
end
OKToStoreArgs=isempty(Args)|UserDataArgsOK;

if isempty(Args)&&(~UserDataArgsOK)
    Args=struct('Holdaxis',0, ...
        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
        'rows',[],'cols',[]); 
end
Args=parseArgs(varargin,Args,{'Holdaxis'},{'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});

if (length(Args.NumericArguments)>2)
    Args.rows=Args.NumericArguments{1};
    Args.cols=Args.NumericArguments{2};
%remove these 2 numerical arguments
    Args.NumericArguments={Args.NumericArguments{3:end}};
end

if OKToStoreArgs
    set(f,'UserData',Args);
end


switch length(Args.NumericArguments)
   case 0
       return % no arguments but rows/cols.... 
   case 1
       if numel(Args.NumericArguments{1}) > 1 % restore subplot(m,n,[x y]) behaviour
           [x1 y1] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(1)); % subplot and ind2sub count differently (column instead of row first) --> switch cols/rows
           [x2 y2] = ind2sub([Args.cols Args.rows],Args.NumericArguments{1}(end));
       else
           x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
           y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
       end
%       x1=mod((Args.NumericArguments{1}-1),Args.cols)+1; x2=x1;
%       y1=floor((Args.NumericArguments{1}-1)/Args.cols)+1; y2=y1;
   case 2
      x1=Args.NumericArguments{1};x2=x1;
      y1=Args.NumericArguments{2};y2=y1;
   case 4
      x1=Args.NumericArguments{1};x2=x1+Args.NumericArguments{3}-1;
      y1=Args.NumericArguments{2};y2=y1+Args.NumericArguments{4}-1;
   otherwise
      error('subaxis argument error')
end
    

cellwidth=((1-Args.MarginLeft-Args.MarginRight)-(Args.cols-1)*Args.SpacingHorizontal)/Args.cols;
cellheight=((1-Args.MarginTop-Args.MarginBottom)-(Args.rows-1)*Args.SpacingVertical)/Args.rows;
xpos1=Args.MarginLeft+Args.PaddingLeft+cellwidth*(x1-1)+Args.SpacingHorizontal*(x1-1);
xpos2=Args.MarginLeft-Args.PaddingRight+cellwidth*x2+Args.SpacingHorizontal*(x2-1);
ypos1=Args.MarginTop+Args.PaddingTop+cellheight*(y1-1)+Args.SpacingVertical*(y1-1);
ypos2=Args.MarginTop-Args.PaddingBottom+cellheight*y2+Args.SpacingVertical*(y2-1);

if Args.Holdaxis
    h=axes('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
else
    h=subplot('position',[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1]);
end


set(h,'box','on');
%h=axes('position',[x1 1-y2 x2-x1 y2-y1]);
set(h,'units',get(gcf,'defaultaxesunits'));
set(h,'tag','subaxis');

if (nargout==0), clear h; 
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ArgStruct=parseArgs(args,ArgStruct,varargin)
% Helper function for parsing varargin. 
%
%
% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
%
% * ArgStruct is the structure full of named arguments with default values.
% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
% * Aliases can be used to map one argument-name to several argstruct fields
%
%
% example usage: 
% --------------
% function parseargtest(varargin)
%
% %define the acceptable named arguments and assign default values
% Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]); 
%
% %The capital letters define abrreviations.  
% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0) 
%
% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
% disp(Args)
%
%
%
%
% Aslak Grinsted 2004

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2double(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (length(varargin)>0) 
    FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
    if length(varargin)>1
        Aliases=varargin{2};
    end
end
 

%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2}); %#ok
FnamesAbbr=strvcat(Fnames{:,3}); %#ok

if length(Aliases)>0  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
end

%--------------get parameters--------------------
l=NumArgCount+1; 
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(a,FnamesFull); 
        end
        if (length(FieldIdx)>1) %shortest fieldname should win 
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
        end
        
        if isempty(FieldIdx) 
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0; 
                    end
                else
                    val=true;
                    paramHasValue=0; 
                end
            end
            if matlabver>=6
                ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
            else
                ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
            end
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end
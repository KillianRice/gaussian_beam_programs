function BeamProfileFits
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
%   - Distances are expected in centimeters 
%   - Beam diameter is expected in microns
%   - Data format should order columns as 
%       Col. 1        | Col. 2          | Col. 3
%       Distance (cm) | Vert waist (um) | Horz waist (um)
%
% Fitting Model
%   fittype('sqrt(a^2+((lambda/pi)^2/a^2)*(x-c)^2)' from Kogelnik and Li

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variables
% For saving you can enter the direct path to your dataset here. If a path
% is not present then one will be asked for.
lambda  = '';
dir     = ''; % Full path to folder containing fitfile
fitfile = ''; % Name of fitfile containing data

% Initial guesses for fitting
guessPosVert     = 1;      %vertical waist position best guess for fit (cm)
guessWaistVert   = 105;    %vertical waist size best guess for fit (um)

guessPosHoriz    = 1;      %horizontal waist position best guess for fit (cm)
guessWaistHoriz  = 105;    %horizontal waist size best guess for fit (um)

lowBoundFitWaist = 0;       %in um
upBoundFitWaist  = 5000;

lowBoundFitPos   = -1e4;    %in cm
upBoundFitPos    =  1e4;


% Set fit options
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[lowBoundFitWaist*1e-4 lowBoundFitPos],...
    'Upper',[upBoundFitWaist*1e-4  upBoundFitPos],...
    'Robust','bisquare');

%% Process raw data and plot
while isempty(lambda)
    % Ask for raw data
    lambda = input('Enter wavelength (nm): ')*10^-7; %wavelength in cm
    if isempty(lambda) || isnan(lambda)
        error('Invalid selection. Please run again.');
    else
        break
    end
end

if isempty(fitfile) || isempty(dir)
    % Ask for raw data
    [fitfile,dir] = uigetfile('*.txt','Choose the fit file');
end

% Get raw data in matlab
dataFileID = fopen([dir fitfile]);
rawData = textscan(dataFileID,'%f%f%f','commentstyle','%'); fclose(dataFileID);

% Convert to Units for fitting
D      = rawData{1};     % Expect distance in centimeters
Vwaist = rawData{2}./2;  % divide by 2 for diameter to radius (gaussian fit)
Wwaist = rawData{3}./2;

%% Plot raw data
figure
plot(D,Wwaist,'k*',D,Vwaist,'r+')
legend('Horizontal','Vertical')
xlabel('Distance From Source(cm)')
ylabel('Beam Waist Size (cm)')
title('Beam Waist Profile')

%% Fitting to gaussian beam equation
% This function was automatically generated on 12-Nov-2004 11:15:57

% Set up figure to receive datasets and fits
figure; grid on
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_   = subplot(1,1,1);
set(ax_,'Box','on');
axes(ax_); hold on;

% --- Plot data originally in dataset "vwaist vs. D"
h_ = line(D,Vwaist,'Parent',ax_,'Color',[0.333333 0 0.666667],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','s', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(D));
xlim_(2) = max(xlim_(2),max(D));
legh_(end+1) = h_;
legt_{end+1} = 'vertical';
 
% --- Plot data originally in dataset "wwaist vs. D"
h_ = line(D,Wwaist,'Parent',ax_,'Color',[0.333333 0.666667 0],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','+', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(D));
xlim_(2) = max(xlim_(2),max(D));
legh_(end+1) = h_;
legt_{end+1} = 'horizontal';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end

% create x vector to plots fit points with
x_pnts = linspace(xlim_(1),xlim_(end),1e3);

%% --- Create fit - Vertical Axis Fit
st_ = [guessWaistVert*1e-4 guessPosVert];
set(fo_,'Startpoint',st_);
ft_ = fittype(strcat('sqrt(waist^2+((',num2str(lambda/pi),')^2/waist^2)*(x-pos)^2)') ,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'waist', 'pos'});

% Fit vertical model using new data
cf_vert      = fit(D,Vwaist*1e-4,ft_ ,fo_);
cf_vert_spot = cf_vert(x_pnts)*1e4;

% Plot Vertical fit
h_ = plot(x_pnts,cf_vert_spot,'r');
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','--', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'vertical fit';

%% --- Create fit - Horizontal Fit
st_ = [guessWaistHoriz*1e-4 guessPosHoriz];
set(fo_,'Startpoint',st_);
ft_ = fittype(strcat('sqrt(waist^2+((',num2str(lambda/pi),')^2/waist^2)*(x-pos)^2)') ,...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'waist', 'pos'});

% Fit horizontal model using new data
cf_horz = fit(D,Wwaist*1e-4,ft_ ,fo_);
cf_horz_spot = cf_horz(x_pnts)*1e4;

% Plot this Horizontal fit
h_ = plot(x_pnts,cf_horz_spot,'b');
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0 0 1],...
     'LineStyle','--', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'horizontal fit';

%% Calculate uncertainties
horz_unc     = confint(cf_horz,.99);
horz_pos_unc = horz_unc(:,2); horz_was_unc = horz_unc(:,1)*1e4;
horz_pos_unc = diff(horz_pos_unc)/2; horz_was_unc = diff(horz_was_unc)/2;

vert_unc     = confint(cf_vert,.99);
vert_pos_unc = vert_unc(:,2); vert_was_unc = vert_unc(:,1)*1e4;
vert_pos_unc = diff(vert_pos_unc)/2; vert_was_unc = diff(vert_was_unc)/2;

% Table information
horz_pos_str = sprintf('Horizontal Waist Position (cm): %g +/- %g',cf_horz.pos,horz_pos_unc);
horz_was_str = sprintf('Horizontal Waist Size     (um): %g +/- %g',cf_horz.waist*1e4,horz_was_unc);
vert_pos_str = sprintf('Vertical Waist Position   (cm): %g +/- %g',cf_vert.pos,vert_pos_unc);
vert_was_str = sprintf('Vertical Waist Size       (um): %g +/- %g',cf_vert.waist*1e4,vert_was_unc);

data = {vert_pos_str; vert_was_str; ''; horz_pos_str; horz_was_str;};

% Add legend and information to fit plot
hold on;
legend(ax_,legh_, legt_,4,'Location','Best');
annotation('textbox',[.15,.8,.1,.1],'String',data,...
    'LineStyle','-','BackgroundColor',[1 1 .84],'EdgeColor',[0 0 0]);
    
xlabel('Distance From Source (cm)','FontSize',14)
ylabel('Beam Waist Size (\mum)','FontSize',14)
title('Beam Waist Profiles','FontSize',14)
plotedit on

%% Show Theoretical beam propagation 
xTheo = -500:500; %cm
horzTheo = sqrt(cf_horz.waist^2 + ((lambda/pi)^2/cf_horz.waist^2).*(xTheo - cf_horz.pos).^2);
vertTheo = sqrt(cf_vert.waist^2 + ((lambda/pi)^2/cf_vert.waist^2).*(xTheo - cf_vert.pos).^2);

figure
plot(xTheo,horzTheo.*1e4,'-b',xTheo,vertTheo.*1e4,'-r');
grid on
xlabel('Distance [cm]');
ylabel('Waist [um]');
legend('Horizontal','Vertical','Location','Best');

%% Clean up workspace but keep fit objects
clearvars -except cf_vert cf_horz
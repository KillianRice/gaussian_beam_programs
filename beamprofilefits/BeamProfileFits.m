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
guessPosVert     = -25;      %vertical waist position best guess for fit (cm)
guessWaistVert   = 50;   %vertical waist size best guess for fit (um)

guessPosHoriz    = -10;      %horizontal waist position best guess for fit (cm)
guessWaistHoriz  = 50;   %horizontal waist size best guess for fit (um)

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
rawData    = textscan(dataFileID,'%f%f%f','commentstyle','%'); fclose(dataFileID);

% Convert to Units for fitting
D      = rawData{1};     % Expect distance in centimeters
Vwaist = rawData{2}./2;  % divide by 2 for diameter to radius (gaussian fit)
Wwaist = rawData{3}./2;

%% Define fit function
rayLength = @(w0) pi*w0^2/lambda;
waistFunc = @(coeffs,z) coeffs(1)*1e-4.*sqrt(1+((z - coeffs(2))./rayLength(coeffs(1)*1e-4)).^2);

%% Plot raw data
dataFig = figure;
plot(D,Wwaist,'b*',D,Vwaist,'r+')
legend('Horizontal','Vertical')
xlabel('Distance From Source(cm)')
ylabel('Beam Waist Size (cm)')
title('Beam Waist Profile')

%% Plotting original data
% This function was automatically generated on 12-Nov-2004 11:15:57

% Set up figure to receive datasets and fits
figure; grid on
ax_     = subplot(2,1,1);
legh_   = []; legt_ = {};   % handles and text for legend
xlim_   = [Inf -Inf];       % limits of x axis
set(ax_,'Box','on');
axes(ax_); hold on;

% --- Plot data originally in dataset "vwaist vs. D"
h_ = line(D,Vwaist,'Parent',ax_,'Color','k',...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','o', 'MarkerSize',12,...
     'MarkerFaceColor',[0.847 0.161 0]);
ylabel('Vertical Waist Radius','FontSize',14)
xlim_(1) = min(xlim_(1),min(D));
xlim_(2) = max(xlim_(2),max(D));
legh_(end+1) = h_;
legt_{end+1} = 'vertical';
 
% --- Plot data originally in dataset "wwaist vs. D"
ax_     = subplot(2,1,2);
h_ = line(D,Wwaist,'Parent',ax_,'Color','k',...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','o', 'MarkerSize',12,...
     'MarkerFaceColor',[0.847 0.161 0]);
ylabel('Horiztonal Waist Radius','FontSize',14)
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

%% --- Create fit - Vertical Fit
beta = [guessWaistVert guessPosVert];

% Fit horizontal model using new data
[cf_Vert,R_Vert,J_Vert] = nlinfit(D,Vwaist*1e-4,waistFunc,beta);

% Plot this Horizontal fit
subplot(2,1,1); hold on
h_ = plot(x_pnts,waistFunc(cf_Vert,x_pnts)*1e4);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0 0 1],...
     'LineStyle','--', 'LineWidth',3,...
     'Marker','none', 'MarkerSize',6,...
     'Color',[0.502 0.502 0.502]);
legh_(end+1) = h_(1);
legt_{end+1} = 'vertical fit';

%% --- Create fit - Horizontal Fit
beta = [guessWaistHoriz guessPosHoriz];

% Fit horizontal model using new data
[cf_Horz,R_Horz,J_Horz] = nlinfit(D,Wwaist*1e-4,waistFunc,beta);

% Plot this Horizontal fit
subplot(2,1,2); hold on
h_ = plot(x_pnts,waistFunc(cf_Horz,x_pnts)*1e4);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0 0 1],...
     'LineStyle','--', 'LineWidth',3,...
     'Marker','none', 'MarkerSize',6,...
     'Color',[0.502 0.502 0.502]);
legh_(end+1) = h_(1);
legt_{end+1} = 'horizontal fit';

%% Calculate uncertainties
horz_unc     = nlparci(cf_Horz,R_Horz,'jacobian',J_Horz);
horz_pos_unc = horz_unc(2,:); horz_was_unc = horz_unc(1,:);
horz_pos_unc = diff(horz_pos_unc)/2; horz_was_unc = diff(horz_was_unc)/2;

vert_unc     = nlparci(cf_Vert,R_Vert,'jacobian',J_Vert);
vert_pos_unc = vert_unc(2,:); vert_was_unc = vert_unc(1,:);
vert_pos_unc = diff(vert_pos_unc)/2; vert_was_unc = diff(vert_was_unc)/2;

% Table information
horz_pos_str = sprintf('Horizontal Waist Position (cm): %g +/- %g',cf_Horz(2),horz_pos_unc);
horz_was_str = sprintf('Horizontal Waist Size     (um): %g +/- %g',cf_Horz(1),horz_was_unc);
vert_pos_str = sprintf('Vertical Waist Position   (cm): %g +/- %g',cf_Vert(2),vert_pos_unc);
vert_was_str = sprintf('Vertical Waist Size       (um): %g +/- %g',cf_Vert(1),vert_was_unc);

data = {vert_pos_str; vert_was_str; ''; horz_pos_str; horz_was_str;};

% Add legend and information to fit plot
hold on;
%legend(ax_,legh_, legt_,4,'Location','Best');
annotation('textbox',[.15,.8,.1,.1],'String',data,...
    'LineStyle','-','BackgroundColor',[1 1 .84],'EdgeColor',[0 0 0]);
    
xlabel('Distance From Source (cm)','FontSize',14)

%% Show Theoretical beam propagation 
xTheo = linspace(-500,500,5e3); %cm
horzTheo = waistFunc(cf_Horz,xTheo)*1e4;
vertTheo = waistFunc(cf_Vert,xTheo)*1e4;

figure(dataFig); hold on
plot(xTheo,horzTheo,'-b',xTheo,vertTheo,'-r');
grid on
xlabel('Distance [cm]');
ylabel('Waist [um]');
legend('Horizontal','Vertical','Location','Best');

%% Clean up workspace but keep fit objects
clearvars -except cf_vert cf_horz
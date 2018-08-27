function sysStruct = GaussBeamProp(sysStruct)
%GaussBeamProp propagates Gaussian beam through optical system of lenses.
%All the equations used can be found in Kogelnik and Li (1968)
%Remember to use the expected units when inputing data. All units expected
%are shown below.
%   [] = GaussBeamProp(sysStruct) - input optical system defined in
%   structure. Final waist information printed to command window and graph
%   shown of the beams propagation through the optical system
%
%   [sysStruct] = GaussBeamProp - If a structure is not specified at call then the
%   user is asked to build the structure. This structure can then be
%   passed to the user to modify and run the program again.
%
% INPUTS
%   sysStruct - structure containing the setup information of the optical
%               system.
%               Structure must contain the following fieldnames
%                    lambda  - wavelength of light [nm]
%                    z_range - [min max] range to be plotted, if one number then the
%                              farthest waist position is taken as min and
%                              the z_range as max
%                    Wst_H_Size   - radius of horizontsl input waist [um]
%                    Wst_H_Pos    - position of horizontal input waist [cm]
%                    h_lens       - cell containing 2 element vectors describing
%                                   each horizontal lens as
%                                   [focal_length(mm),lens_position(cm)]
%                    Wst_V_Size   - (optional) same as Wst_H_Size
%                    Wst_V_Pos    - (optional) same as Wst_H_Pos
%                    v_lens       - (optional) same as horz_lens
%                    full_curves  - (optional) Boolean to plot propagation of each
%                                   beam along full range. Default is 0
%                    show_lenses  - (optional) Boolean to plot a line showing
%                                   positions of lenses in the system.
%                                   Default is 1
%                    fig_handle   - (optional) Plot data on specified figure
%                                   
% OUTPUTS
%   sysStruct - structure containg the setup information of the optical
%               system. When run without inputs this output can be
%               specified to capture the structure that will be created.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Jim Aman - 9.9.13
% Updates
%   9.11.13 - Add input structure to specify system and reorgaanize code
%   9.14.13 - Sort lens, calculate beams outside of specified range
%   9.19.13 - Fix bug to plot waist with no lenses, add ability to specify
%             figure handle to plot data on

% Future Work
% - check that waist are included for both beam types

%% Check input arguments
% if no input arguments create system structure
if nargin ~= 1
    sysStruct = sysStructSetup;
end

%% Setup the workspace and relevant variables
% Initialize default value variables
options.lineColor   = {'b','r'}; % used to identify horizontal and vertical beams
options.lineWidth   = 2.5;       % Thicker lines
options.full_curves = 1;         % full_curve default
options.show_lenses = 1;         % show_lenses default
options.precision   = 1e-2;      % centimeters, precision of grid (space between pts)
options.compPow     = 1e1;       % power to convert to integer for vector comparison (compPow*1/precision)
options.lensHeight  = 50;        % size of lens to show on plot (um)
options.FontSize    = 16;        % Fontsize for axes

% Assign unchanged loop variables from structure
lambda = sysStruct.lambda*1e-9;
h_opt  = isfield(sysStruct,'h_lens');
v_opt  = isfield(sysStruct,'v_lens');

if isfield(sysStruct,'full_curves'); options.full_curves  = sysStruct.full_curves; end
if isfield(sysStruct,'show_lenses'); options.show_lenses  = sysStruct.show_lenses; end
options.h_opt = h_opt; options.v_opt = v_opt; % tells beamTransfer how to write legend

%% Error checking to be sure appropriate fields are provided
if ~isfield(sysStruct,'lambda');
    error('Wavelength not specified in the input structure.')
end

if ~isfield(sysStruct,'z_range')
    error(['No range provided in the input structure. Please provide the ',...
        'upper bound of the calculation and try again.'])
end

% check at least one of h_lens or v_lens is present
if ~(h_opt || v_opt)
    error(['Neither h_lens or v_lens was found in the structure provided. Please ',...
        'check the input structure and try again.'])
end

%% z Limit assignment
if numel(sysStruct.z_range) == 2
    z_lim = sysStruct.z_range*1e-2;
elseif h_opt && v_opt
    z_lim = [min(sysStruct.Wst_H_Pos,sysStruct.Wst_V_Pos) sysStruct.z_range]*1e-2;
elseif h_opt
    z_lim = [sysStruct.Wst_H_Pos sysStruct.z_range]*1e-2;
elseif v_opt
    z_lim = [sysStruct.Wst_V_Pos sysStruct.z_range]*1e-2;
else
    error('An unexpected event occured in the assignment of z_lim')
end

%% Setup Figure
% if fig_handle exists and is not zero then plot on same figure
if isfield(sysStruct,'fig_handle')
    if sysStruct.fig_handle
        figure(sysStruct.fig_handle)
    end
else
    figure 
end
hold on; grid on
    
%% Setup and propogate each waist given
    if h_opt
        lens = sysStruct.h_lens;
        Wst_Size = sysStruct.Wst_H_Size*1e-6;
        Wst_Pos  = sysStruct.Wst_H_Pos*1e-2;
        options.beam_type = 'Horizontal';
        
        beamTransform(lambda,lens,z_lim,Wst_Size,Wst_Pos,options);
    end
    
    if v_opt
        lens = sysStruct.v_lens;
        Wst_Size = sysStruct.Wst_V_Size*1e-6;
        Wst_Pos  = sysStruct.Wst_V_Pos*1e-2;
        options.beam_type = 'Vertical';
        
        beamTransform(lambda,lens,z_lim,Wst_Size,Wst_Pos,options);
    end
    
    %% Set Figure Properties
    xlim(z_lim*1e2)
    xlabel('[cm]','FontSize',options.FontSize);
    ylabel('Beam Radius [\mum]','FontSize',options.FontSize)
    legend('show','Location','Best')
end

function beamTransform(lambda,lens,z_lim,Wst_Size,Wst_Pos,options)
%beamTransform calculates the beam transformation through the system of
%lens with given focal length and position

% Preallocate variables
q = zeros(1,length(lens) + 1);
z = zeros(1,length(lens) + 2);

% Sort lens positions and save into z vector
if ~isempty(lens) % only sort lenses if there are lenses
    lens_tmp   = [lens{:}];               % convert cell to matrix
    lens_mat(:,1) = lens_tmp(1:2:end)*1e-3; lens_mat(:,2) = lens_tmp(2:2:end)*1e-2;
    [~,IX] = sort(lens_mat(:,2)); lens_mat = lens_mat(IX,:); % sort lens in ascending order
    z(2:end-1) = lens_mat(:,2);           % pick out lens positions (convert to meter)
end
z(1:length(z)-1:end) = z_lim;         % assign start and stop position 
% allow for lens to be outside the range specified and still account for them
if find(z == min(z)) ~= 1
    z(1) = min(z) - .1; % Go 10cm in front of first lens
end
if find(z == max(z)) ~= length(z)
    z(end) = max(z) + .1; % Go 10cm after last lens
end

% Initial values for loop through lenses
q(1) = 1i*pi*Wst_Size^2/lambda + (z(1) - Wst_Pos);

%% Loop through optical components and determine spotsize
% This is the meat of the algorithm. Most of the progation mathematics are
% in this loop
for i = 1:length(q)
    if i <= length(lens)
        % Find z values between last object and current object
        zSep       = z(i):options.precision*1e-2:z(i+1);
        % transformation of current beam at lens through freespace
        q_tmp      = (zSep(end) - zSep(1)) + q(i);
        % determine optical system to transform current beam
        opticalSys = [1 0; -1/(lens_mat(i,1)) 1];
        % transform current beam through optical system
        q(i+1)     = (opticalSys(1,:)*[q_tmp; 1])...
                      /(opticalSys(2,:)*[q_tmp; 1]);
    else
        % If last loop then propogate to the end of the specified range
        zSep = linspace(z(i),z(end),abs(z(i+1) - z(i))*1e4);
    end
    
%% Setup plotting variables
    if strcmpi(options.beam_type,'Horizontal')
        plot_color = options.lineColor(1);
    elseif strcmpi(options.beam_type,'Vertical')
        plot_color = options.lineColor(2);
    end
    
    % Before plotting if there are two axes plot point to identify color
    % with each axis
    if i == 1 && options.h_opt && options.v_opt
        h = plot(0,0,char(plot_color));
        set(h,'DisplayName',options.beam_type,'Visible','off','Marker','*')
    end
    
    if options.full_curves
        % if full curves then plot along the entire range
        z_range_calc = z(1):options.precision*1e-2:z(end);
        z_range_plot = z_range_calc; 
    else
        %plot from last object to current object
        z_range_calc = zSep; % account for distance between objects
        z_range_plot = zSep;
    end
    
%% Calculate the beam parameter, spotSize, waist size, and position
    q_freespace  = (z_range_calc - z(i)) + q(i);       % calc q 
    spotSize     = sqrt(-lambda./(pi.*imag(1./q_freespace)));
    R            = 1/real(1/q_freespace(end));         % Radius of curvature is real part of q
    Wst_Out_Size = sqrt(spotSize(end)^2/(1 + (pi*spotSize(end)^2/(lambda*R))^2)); % K&L eq. 24
    Wst_Out_Pos  = z_range_plot(end) - R/(1 + (lambda*R/(pi*spotSize(end)^2))^2); % K&L eq. 25

%% Plotting spotsize vs. position    
    % comparison vectors to convert to integers for accurate comparison
    Z_Range_Plot = round(options.compPow*(1/options.precision)*z_range_plot); 
    ZSep         = round(options.compPow*(1/options.precision)*zSep);
    
    h1 = plot(z_range_plot(ismember(Z_Range_Plot,ZSep))*1e2,...
        spotSize(ismember(Z_Range_Plot,ZSep))*1e6,...
        char(plot_color),'LineWidth',options.lineWidth);
    if i == 1
        set(h1,'DisplayName',['Input    - Radius: ',num2str(round(Wst_Out_Size*1e6)),'\mum',...
            ', Z_0: ',num2str(round(Wst_Out_Pos*1e3)*1e-1),'cm'])
    elseif i > 1
        set(h1,'DisplayName',['Lens ',num2str(i-1),' - f: ',num2str(lens_mat(i-1,1)*1e3),'mm, Z_0: ',...
            num2str(lens_mat(i-1,2)*1e2),'cm - Radius: ',num2str(round(Wst_Out_Size*1e6)),...
            '\mum, Z_0: ',num2str(round(Wst_Out_Pos*1e3)*1e-1),'cm'])
    end
    
    % Use a different linestyle to plot the part of the lines outside the
    % path of the actual beam
    if options.full_curves
        h2 = plot(z_range_plot(~ismember(Z_Range_Plot,ZSep))*1e2,...
            spotSize(~ismember(Z_Range_Plot,ZSep))*1e6,...
            char(plot_color),'LineWidth',options.lineWidth*1e-1);
        set(h2,'LineStyle',':'); %Change linestyle to differentiate from main path
        set(get(get(h2,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off'); % Exclude line from legend
    end
    
    % Plot lines showing lenses
    if options.show_lenses && i > 1
        yLensPos = round(spotSize(single(z_range_plot) == single(z(i)))*1e6);
        h3 = line([z(i) z(i)]*1e2,...
            [(yLensPos + options.lensHeight) (yLensPos - options.lensHeight)]);
        set(h3,'LineStyle','-','Color','k','LineWidth',options.lineWidth);
        set(get(get(h3,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off'); % Exclude line from legend
    end
end
end

function sysStruct = sysStructSetup
%sysStructSetup is a helper function that is called to create the needed
%sysStruct structure if GaussBeamProp is called without input arguements

sysStruct.lambda = input('Wavelength (nm): ');
sysStruct.z_range = input('Range to plot over (cm) [min max]: ');
sysStruct.full_curves = input('Show full curves on plot? Yes (1), No (0): ');
sysStruct.show_lenses = input('Show lens positions on plot? Yes(1), No (0): ');

while 1
    waist_types = input('What types of waists are you dealing with? Horz (1), Vert (2), Both (3): ');
    switch waist_types
        case 1
            sysStruct.Wst_H_Size = input('Horizontal waist size (um): ');
            sysStruct.Wst_H_Pos  = input('Horizontal waist position (cm): ');
            num_h_lens = input('How many horizontal lenses?: ');
            tmp_cell = cell(1,num_h_lens);
            for i = 1:num_h_lens
                focal_length = input(sprintf('Focal length of lens %g (mm): ',i));
                lens_position = input(sprintf('Position of lens %g (cm): ',i));
                tmp_cell{i} = [focal_length lens_position];
            end
            sysStruct.h_lens = tmp_cell;
            
            break
        case 2
            sysStruct.Wst_V_Size = input('Vertical waist size (um): ');
            sysStruct.Wst_V_Pos  = input('Vertical waist position (cm): ');
            num_h_lens = input('How many vertical lenses?: ');
            tmp_cell = cell(1,num_h_lens);
            for i = 1:num_h_lens
                focal_length = input(sprintf('Focal length of lens %g (mm): ',i));
                lens_position = input(sprintf('Position of lens %g (cm): ',i));
                tmp_cell{i} = [focal_length lens_position];
            end
            sysStruct.h_lens = tmp_cell;
            
            break
        case 3
            sysStruct.Wst_H_Size = input('Horizontal waist size (um): ');
            sysStruct.Wst_H_Pos  = input('Horizontal waist position (cm): ');
            sysStruct.Wst_V_Size = input('Vertical waist size (um): ');
            sysStruct.Wst_V_Pos  = input('Vertical waist position (cm): ');
            
            lens_type = input('Are you considering cylindrical(1), spherical(2), or a mixture(3) of lenses?: ');
            if lens_type == 1 || lens_type == 3
                num_h_lens = input('How many horizontal lenses?: ');
                tmp_cell = cell(1,num_h_lens);
                for i = 1:num_h_lens
                    focal_length = input(sprintf('Focal length of lens %g (mm): ',i));
                    lens_position = input(sprintf('Position of lens %g (cm): ',i));
                    tmp_cell{i} = [focal_length lens_position];
                end
                sysStruct.h_lens = tmp_cell;
                num_v_lens = input('How many vertical lenses?: ');
                tmp_cell = cell(1,num_v_lens);
                for i = 1:num_v_lens
                    focal_length = input(sprintf('Focal length of lens %g (mm): ',i));
                    lens_position = input(sprintf('Position of lens %g (cm): ',i));
                    tmp_cell{i} = [focal_length lens_position];
                end
                sysStruct.v_lens = tmp_cell;
            else
                num_h_lens = input('How many lenses?: ');
                tmp_cell = cell(1,num_h_lens);
                for i = 1:num_h_lens
                    focal_length = input(sprintf('Focal length of lens %g (mm): ',i));
                    lens_position = input(sprintf('Position of lens %g (cm): ',i));
                    tmp_cell{i} = [focal_length lens_position];
                end
                sysStruct.h_lens = tmp_cell;
                sysStruct.v_lens = tmp_cell;
            end
            break
        otherwise
            disp('Incorrect option chosen. Please try again'); drawnow
    end

% questdlg('What types of waists are you dealing with?','Beam Type','Horizontal','Vertical','Both','Horizontal')
% inputdlg('Please input the wavelength (nm)','Wavelegth',1)
end
end
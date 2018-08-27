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
%                    h_comp       - cell containing 2 element vectors describing
%                                   each horizontal component as
%                                   [comp specific  , comp_position(cm)]
%                                   Ex. (1) lens - h_comp(1) is the focal length in mm
%                                       (2) afocal magnifcation - h_comp(1) is magnification  
%                    h_comp_type  - 1D vector of integers specifying the particular
%                                   components of the system
%                                       1 - lens
%                                       2 - afocal magnification (prism pair)
%                    Wst_V_Size   - (optional) same as Wst_H_Size
%                    Wst_V_Pos    - (optional) same as Wst_H_Pos
%                    v_comp       - (optional) same as horz_comp
%                    v_comp_type  - (optional) same as h_comp_type
%                    full_curves  - (optional) Boolean to plot propagation of each
%                                   beam along full range. Default is 0
%                    show_comps   - (optional) Boolean to plot a line showing
%                                   positions of components in the system.
%                                   Default is 1
%                    fig_handle   - (optional) Plot data on specified figure
%                    axes_handle  - (optional) Plot data on specified axes of specified figure
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
%   10.9.13 - Added ability to specify transfer matricies other than lenses
%             Built-in options are for lenses or afocal magnification (prism pairs)
%   8.15.15 - Added ability to plot on specified axes_handle

% Future Work
% - check that waist are included for both beam types
% - Figure out how to propagate with a non-unity M^2 value

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
options.show_comps  = 1;         % show_comps default
options.precision   = 1e-2;      % centimeters, precision of grid (space between pts)
options.compPow     = 1e1;       % power to convert to integer for vector comparison (compPow*1/precision)
options.compHeight  = 50;        % size of components to show on plot (um)
options.FontSize    = 16;        % Fontsize for axes
options.compSpecs   = {'f' 'mm'; 'Mag' 'X'}; % Component specific information and units

% Assign unchanged loop variables from structure
lambda = sysStruct.lambda*1e-9;
h_opt  = isfield(sysStruct,'h_comp');
v_opt  = isfield(sysStruct,'v_comp');

if isfield(sysStruct,'full_curves'); options.full_curves  = sysStruct.full_curves; end
if isfield(sysStruct,'show_comps'); options.show_comps  = sysStruct.show_comps; end
options.h_opt = h_opt; options.v_opt = v_opt; % tells beamTransfer how to write legend

%% Error checking to be sure appropriate fields are provided
if ~isfield(sysStruct,'lambda');
    error('Wavelength not specified in the input structure.')
end

if ~isfield(sysStruct,'z_range')
    error(['No range provided in the input structure. Please provide the ',...
        'upper bound of the calculation and try again.'])
end

% check at least one of h_comp or v_comp is present
if ~(h_opt || v_opt)
    error(['Neither h_comp or v_comp was found in the structure provided. Please ',...
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

% if axes_handle exists and is not zero then plot on specified axes
if isfield(sysStruct,'axes_handle')
    if sysStruct.axes_handle
        axes(sysStruct.axes_handle); 
        set(sysStruct.axes_handle,'Box','on','Linewidth',2);
    end
else
    set(gca,'Box','on','Linewidth',2);
end
hold on; grid on
    
%% Setup and propogate each waist given
    if h_opt
        if ~isfield(sysStruct,'h_comp_type')
            error(['Horizontal optical component types not specified in input structure. Please ',...
                'provide the component types and try again.'])
        end
        
        comp      = sysStruct.h_comp;
        comp_type = sysStruct.h_comp_type;
        Wst_Size  = sysStruct.Wst_H_Size*1e-6;
        Wst_Pos   = sysStruct.Wst_H_Pos*1e-2;
        options.beam_type = 'Horizontal';
        
        beamTransform(lambda,comp,comp_type,z_lim,Wst_Size,Wst_Pos,options);
    end
    
    if v_opt
        if ~isfield(sysStruct,'v_comp_type')
            error(['Vertical optical component types not specified in input structure. Please ',...
                'provide the component types and try again.'])
        end
        
        comp      = sysStruct.v_comp;
        comp_type = sysStruct.v_comp_type;
        Wst_Size  = sysStruct.Wst_V_Size*1e-6;
        Wst_Pos   = sysStruct.Wst_V_Pos*1e-2;
        options.beam_type = 'Vertical';
        
        beamTransform(lambda,comp,comp_type,z_lim,Wst_Size,Wst_Pos,options);
    end
    
    %% Set Figure Properties
    xlim(z_lim*1e2)
    xlabel('[cm]','FontSize',options.FontSize);
    ylabel('Beam Radius [\mum]','FontSize',options.FontSize)
    legend('show','Location','Best')
end

function beamTransform(lambda,comp,comp_type,z_lim,Wst_Size,Wst_Pos,options)
%beamTransform calculates the beam transformation through the system of
%componenets with given focal length and position

% Preallocate variables
q = zeros(1,length(comp) + 1);
z = zeros(1,length(comp) + 2);

% Sort comp positions and save into z vector
if ~isempty(comp) % only sort components if there are components
    comp_tmp   = [comp{:}];               % convert cell to matrix
    comp_mat(:,1) = comp_tmp(1:2:end)*1e-3; comp_mat(:,2) = comp_tmp(2:2:end)*1e-2;
    [tmp,IX] = sort(comp_mat(:,2)); comp_mat = comp_mat(IX,:); % sort comp in ascending order
    z(2:end-1) = comp_mat(:,2);           % pick out comp positions (convert to meter)
end
z(1:length(z)-1:end) = z_lim;         % assign start and stop position 
% allow for comp to be outside the range specified and still account for them
if find(z == min(z)) ~= 1
    z(1) = min(z) - .1; % Go 10cm in front of first comp
end
if find(z == max(z)) ~= length(z)
    z(end) = max(z) + .1; % Go 10cm after last comp
end

% Initial values for loop through components
q(1) = 1i*pi*Wst_Size^2/lambda + (z(1) - Wst_Pos);

%% Loop through optical components and determine spotsize
% This is the meat of the algorithm. Most of the progation mathematics are
% in this loop
for i = 1:length(q)
    if i <= length(comp)
        % Find z values between last object and current object
        zSep       = z(i):options.precision*1e-2:z(i+1);
        % transformation of current beam at comp through freespace
        q_tmp      = (zSep(end) - zSep(1)) + q(i);
        % determine optical system to transform current beam
        transformMat = transformMatSelect(comp{i},comp_type(i));
        % transform current beam through optical system
        q(i+1)     = (transformMat(1,:)*[q_tmp; 1])...
                      /(transformMat(2,:)*[q_tmp; 1]);
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
    q_freespace  = (z_range_calc - z(i)) + q(i);     % calc q 
    spotSize     = sqrt(-lambda./(pi.*imag(1./q_freespace)));
    R            = 1/real(1/q_freespace(end));       % Radius of curvature is real part of 1/q
    Wst_Out_Size = sqrt(spotSize(end)^2/(1 + (pi*spotSize(end)^2/(lambda*R))^2)); % K&L eq. 24
    Wst_Out_Pos  = z_range_plot(end) - R/(1 + (lambda*R/(pi*spotSize(end)^2))^2); % K&L eq. 25
        % Eq. 25 gives distance from waist to current point, then reference
        % your known position to get position of waist for current axes.

%% Plotting spotsize vs. position    
    % comparison vectors to convert to integers for accurate comparison
    Z_Range_Plot = round(options.compPow*(1/options.precision)*z_range_plot); 
    ZSep         = round(options.compPow*(1/options.precision)*zSep);
    
    h1 = plot(z_range_plot(ismember(Z_Range_Plot,ZSep))*1e2,...
        spotSize(ismember(Z_Range_Plot,ZSep))*1e6,...
        char(plot_color),'LineWidth',options.lineWidth);
    if i == 1
        set(h1,'DisplayName',['Input     - Radius: ',num2str(round(Wst_Out_Size*1e6)),'\mum',...
            ', Z_0: ',num2str(round(Wst_Out_Pos*1e3)*1e-1),'cm'])
    elseif i > 1
        set(h1,'DisplayName',['Comp ',num2str(i-1),' - ',options.compSpecs{comp_type(i-1),1},...
            ': ',num2str(comp_mat(i-1,1)*1e3),options.compSpecs{comp_type(i-1),2},', Z: ',...
            num2str(comp_mat(i-1,2)*1e2),'cm - Radius: ',num2str(round(Wst_Out_Size*1e6)),...
            '\mum, Z_0: ',num2str(round(Wst_Out_Pos*1e3)*1e-1),'cm'])
    end
    
    % Use a different linestyle to plot the part of the lines outside the
    % path of the actual beam
    if options.full_curves
        h2 = plot(z_range_plot*1e2,spotSize*1e6,...
            char(plot_color),'LineWidth',options.lineWidth*1e-1);
        set(h2,'LineStyle',':'); %Change linestyle to differentiate from main path
        set(get(get(h2,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off'); % Exclude line from legend
    end
    
    % Plot lines showing components
    if options.show_comps && i > 1
        yCompPos = round(spotSize(single(z_range_plot) == single(z(i)))*1e6);
        h3 = line([z(i) z(i)]*1e2,...
            [(yCompPos + options.compHeight) (yCompPos - options.compHeight)]);
        set(h3,'LineStyle','-','Color','k','LineWidth',options.lineWidth);
        set(get(get(h3,'Annotation'),'LegendInformation'),...
         'IconDisplayStyle','off'); % Exclude line from legend
    end
end
end

function transformMat = transformMatSelect(comp,comp_type)
% Select transformation matrix for optical system
%
% Componenet selction guide (comp_type)
%   1 = lens:                 comp = [focal length(mm), position(cm)]
%   2 = afocal magnification: comp = [magnification, position(cm)]

switch comp_type
    case 1 %(Kogelnik and Li)
        transformMat = [    1          0
                        -1e3/comp(1)   1];
    case 2 %(Introduction to Matrix Methods in Optics pg. 78 )
        transformMat = [comp(1)    0
                        0       1/comp(1)];
    otherwise
        error(['Invalid selection of component type. Check the help for ',...
            'supported transfer matricies and try again.'])
end
end

function sysStruct = sysStructSetup
%sysStructSetup is a helper function that is called to create the needed
%sysStruct structure if GaussBeamProp is called without input arguements

sysStruct.lambda  = input('Wavelength (nm): ');
sysStruct.z_range = input('Range to plot over (cm) [min max]: ');
sysStruct.full_curves = input('Show full curves on plot? Yes (1), No (0): ');
sysStruct.show_comps  = input('Show component positions on plot? Yes(1), No (0): ');

while 1
    waist_types = input('What types of waists are you dealing with? Horz (1), Vert (2), Both (3): ');
    switch waist_types
        case 1
            sysStruct.Wst_H_Size = input('Horizontal waist size (um): ');
            sysStruct.Wst_H_Pos  = input('Horizontal waist position (cm): ');
            break
        case 2
            sysStruct.Wst_V_Size = input('Vertical waist size (um): ');
            sysStruct.Wst_V_Pos  = input('Vertical waist position (cm): ');            
            break
        case 3
            sysStruct.Wst_H_Size = input('Horizontal waist size (um): ');
            sysStruct.Wst_H_Pos  = input('Horizontal waist position (cm): ');
            sysStruct.Wst_V_Size = input('Vertical waist size (um): ');
            sysStruct.Wst_V_Pos  = input('Vertical waist position (cm): ');
            break
        otherwise
            disp('Incorrect option chosen. Please try again'); drawnow
    end
end

comp_types = input('What types of components are you dealing with? Lenses (1), Prism pairs (2), Both (3): ');
if comp_types == 1 || comp_types == 3
    switch waist_types
        case 1
            num_h_lens = input('How many horizontal lenses?: ');
            tmp_cell   = cell(1,num_h_lens);
            for i = 1:num_h_lens
                focal_length  = input(sprintf('Focal length of lens %g (mm): ',i));
                lens_position = input(sprintf('Position of lens %g (cm): ',i));
                tmp_cell{i}   = [focal_length lens_position];
            end
            sysStruct.h_comp                    = tmp_cell;
            sysStruct.h_comp_type(1:num_h_lens) = 1;
        case 2
            num_v_lens = input('How many vertical lenses?: ');
            tmp_cell   = cell(1,num_v_lens);
            for i = 1:num_v_lens
                focal_length  = input(sprintf('Focal length of lens %g (mm): ',i));
                lens_position = input(sprintf('Position of lens %g (cm): ',i));
                tmp_cell{i}   = [focal_length lens_position];
            end
            sysStruct.v_comp = tmp_cell;
            sysStruct.v_comp_type(1:num_v_lens) = 1;
        case 3
            while 1
                lens_type = input('Are you considering cylindrical(1) or spherical(2) lenses?: ');
                if lens_type == 1
                    num_h_lens = input('How many horizontal lenses?: ');
                    tmp_cell   = cell(1,num_h_lens);
                    for i = 1:num_h_lens
                        focal_length  = input(sprintf('Focal length of lens %g (mm): ',i));
                        lens_position = input(sprintf('Position of lens %g (cm): ',i));
                        tmp_cell{i}   = [focal_length lens_position];
                    end
                    sysStruct.h_comp                    = tmp_cell;
                    sysStruct.h_comp_type(1:num_h_lens) = 1;
                    
                    num_v_lens = input('How many vertical lenses?: ');
                    tmp_cell = cell(1,num_v_lens);
                    for i = 1:num_v_lens
                        focal_length  = input(sprintf('Focal length of lens %g (mm): ',i));
                        lens_position = input(sprintf('Position of lens %g (cm): ',i));
                        tmp_cell{i}   = [focal_length lens_position];
                    end
                    sysStruct.v_comp                    = tmp_cell;
                    sysStruct.v_comp_type(1:num_v_lens) = 1;
                    break
                elseif lens_type == 2
                    num_h_lens = input('How many lenses?: ');
                    tmp_cell   = cell(1,num_h_lens);
                    for i = 1:num_h_lens
                        focal_length  = input(sprintf('Focal length of lens %g (mm): ',i));
                        lens_position = input(sprintf('Position of lens %g (cm): ',i));
                        tmp_cell{i}   = [focal_length lens_position];
                    end
                    sysStruct.h_comp                    = tmp_cell;
                    sysStruct.v_comp                    = tmp_cell;
                    sysStruct.h_comp_type(1:num_h_lens) = 1;
                    sysStruct.v_comp_type(1:num_h_lens) = 1;
                    break
                else
                    disp('Incorrect option chosen. Please try again'); drawnow
                end
            end
    end
end
if comp_types == 2 || comp_types == 3
    switch waist_types
        case 1
            num_h_comp = input('How many horizontal prism pairs?: ');
            tmp_cell   = cell(1,num_h_comp);
            for i = 1:num_h_comp
                magnification = input(sprintf('Magnification of prism pair %g: ',i));
                pair_position = input(sprintf('Position of prism pair %g (cm): ',i));
                tmp_cell{i}   = [magnification pair_position];
            end
            if isfield(sysStruct,'h_comp')
                sysStruct.h_comp = [sysStruct.h_comp tmp_cell];
                lengthTmp        = length(sysStruct.h_comp_type) + 1;
                sysStruct.h_comp_type(lengthTmp:lengthTmp + num_h_comp - 1) = 2;
            else
                sysStruct.h_comp = tmp_cell;
                sysStruct.h_comp_type(1:num_h_comp) = 2;
            end
        case 2
            num_v_comp = input('How many vertical prism pairs?: ');
            tmp_cell   = cell(1,num_v_comp);
            for i = 1:num_v_comp
                magnification = input(sprintf('Magnification of prism pair %g: ',i));
                pair_position = input(sprintf('Position of prism pair %g (cm): ',i));
                tmp_cell{i}   = [magnification pair_position];
            end
            if isfield(sysStruct,'v_comp')
                sysStruct.v_comp = [sysStruct.v_comp tmp_cell];
                lengthTmp        = length(sysStruct.v_comp_type) + 1;
                sysStruct.v_comp_type(lengthTmp:lengthTmp + num_v_comp - 1) = 2;
            else
                sysStruct.v_comp = tmp_cell;
                sysStruct.v_comp_type(1:num_v_comp) = 2;
            end

        case 3
            num_h_comp = input('How many horizontal prism pairs?: ');
            tmp_cell   = cell(1,num_h_comp);
            for i = 1:num_h_comp
                magnification = input(sprintf('Magnification of prism pair %g: ',i));
                pair_position = input(sprintf('Position of prism pair %g (cm): ',i));
                tmp_cell{i}   = [magnification pair_position];
            end
            if isfield(sysStruct,'h_comp')
                sysStruct.h_comp = [sysStruct.h_comp tmp_cell];
                lengthTmp        = length(sysStruct.h_comp_type) + 1;
                sysStruct.h_comp_type(lengthTmp:lengthTmp + num_h_comp - 1) = 2;
            else
                sysStruct.h_comp = tmp_cell;
                sysStruct.h_comp_type(1:num_h_comp) = 2;
            end
            
            num_v_comp = input('How many vertical prism pairs?: ');
            tmp_cell   = cell(1,num_v_comp);
            for i = 1:num_v_comp
                magnification = input(sprintf('Magnification of prism pair %g: ',i));
                pair_position = input(sprintf('Position of prism pair %g (cm): ',i));
                tmp_cell{i}   = [magnification pair_position];
            end
            if isfield(sysStruct,'v_comp')
                sysStruct.v_comp = [sysStruct.v_comp tmp_cell];
                lengthTmp        = length(sysStruct.v_comp_type) + 1;
                sysStruct.v_comp_type(lengthTmp:lengthTmp + num_v_comp - 1) = 2;
            else
                sysStruct.v_comp = tmp_cell;
                sysStruct.v_comp_type(1:num_v_comp) = 2;
            end
    end
end
end
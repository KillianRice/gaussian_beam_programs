function GaussBeamProp
% Propagates Gaussian beam through optical system of lenses. 
%
% INPUTS
%   lambda - wavelength of light [nm]
%   q      - complex beam parameter describing the beam waist
%   z      - distance from the beam waist to find spot size [m] (can be a
%            vector)
%
% OUTPUTS
%
% Created by Jim Aman - 9.9.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input waist size given as radius in um
% Input waist position given in cm
% Input range to plot over
% Input lens focal length (mm) and position (cm)

% Future Work
% - include ability to specify setup file or run call function with arguments
% if no setup or args then ask for needed info
% - show full curves? yes or no
% - error on z_end less than position of last lens

% Assign variables from variable input if any
lens = {[-100 .1] [200 .2]};
z_end = 50;         %[cm]
lambda = 689;       %[nm]
Waist_1_Size = 250; %[um]
Waist_1_Pos  = 0;   %[cm]

% Convert lengths to meters
lambda       = lambda*1e-9;
Waist_1_Size = Waist_1_Size*1e-6;
Waist_1_Pos  = Waist_1_Pos*1e-2;
z_end        = z_end*1e-2;

% Reserve variable space
q = zeros(1,length(lens) + 1);
z = zeros(1,length(lens) + 2);
figure; hold on

% Initial values for loop through lenses
q(1) = 1i*pi*Waist_1_Size^2/lambda;
z(1) = Waist_1_Pos;

% Save lens positions into z vector
lens_tmp   = [lens{:}];         % convert cell to matrix
z(2:end-1) = lens_tmp(2:2:end); % pick out lens positions
z(end)     = z_end;             % end of z vector is edge of specified range

% Loop through optical components and determine spotsize
for i = 1:length(q)
    if i <= length(lens)
        % Find z values between last object and current object
        zSep       = linspace(z(i),lens{i}(2));
        % transformation of current beam at lens through freespace
        q_tmp      = (zSep(end) - zSep(1)) + q(i);  
        % determine optical system to transform current beam
        opticalSys = [1 0; -1/(lens{i}(1)*1e-3) 1]; 
        % transform current beam through optical system
        q(i+1)     = (opticalSys(1,:)*[q_tmp; 1])...
                      /(opticalSys(2,:)*[q_tmp; 1]); 
    else
        % If last loop then propogate to the end of the specified range
        zSep = linspace(z(i),z(end));
    end
    
    q_freespace = (zSep - z(i)) + q(i); %calc q from last object to current object
    spotSize    = sqrt(-lambda./(pi.*imag(1./q_freespace))); 
    plot(zSep*1e2,spotSize*1e6);
end

% Calculate the waist size and position after the optical system (from
% Kogelnik and Li eq. 24 & 25)
R = 1/real(1/q_freespace(end));
Waist_2_Size = sqrt(spotSize(end)^2/(1 + (pi*spotSize(end)^2/(lambda*R))^2));
Waist_2_Pos  = z(end) - R/(1 + (lambda*R/(pi*spotSize(end)^2))^2);

% Figure Properties
grid on
xlim([Waist_1_Pos*1e6,z_end*1e2])

% Print final waist information
fprintf('Output beam parameters\n\t Waist Position: %g cm \n\t Waist Size:\t %g um\n\n',...
    Waist_2_Pos*1e2,Waist_2_Size*1e6)
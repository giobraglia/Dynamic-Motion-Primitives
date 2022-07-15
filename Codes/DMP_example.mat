%--------------------------------------------------------------------------
%% DMP CARTESIAN MODEL
%--------------------------------------------------------------------------
% This script processes recorded waveforms.mat files demonstrated to the 
% robot. Waveforms are expressed in [m] and are firstly rounded at [10e-4m]
% in order to ensure zero initial velocity.
% After that, waveforms are interpolated to achieve smooth first and second
% derivatives.
% Those values are then used to create [X Y Z] object which is useful to run
% simulation with 'Recorder.slc' in order to retrieve the desire forcing 
% function.
% Once we obtain the desired forcing functions for xyz coordinates, DMP 
% object are ready to be fitted; always run first the 'fit' method, then the
% 'get_forcing_function' one

clear all
clc

%% ------------------------------------------------------------------------
% 1 -Trajectory pre-processing
%--------------------------------------------------------------------------

% Load waveforms that has to be interpolated
%--------------------------------------------------------------------------

wave_x = load('X.mat');
wave_y = load('Y.mat');
wave_z = load('Z.mat');

wave_x = (wave_x.X)' ;
wave_y = (wave_y.Y)' ;
wave_z = (wave_z.Z)' ;

time = wave_x(:,1);
X.x  = wave_x(:,2);
Y.y  = wave_y(:,2);
Z.z  = wave_z(:,2);

clear wave_x wave_y wave_z

Ts   = 0.001                ;  % sampling time
T    = length(wave_x) * Ts  ;  % duration of the desired trajectory
time = ( 0 : Ts : (T-Ts) )' ;  % time array


% Ensure smoothness of derivatives
%--------------------------------------------------------------------------
fx = fit(time, X.x,'smoothingspline', 'SmoothingParam', 0.9 );
[X.dx,X.ddx] = differentiate(fx,time);

fy = fit(time, Y.y,'smoothingspline', 'SmoothingParam', 0.9 );
[Y.dy,Y.ddy] = differentiate(fy,time);

fz = fit(time, Z.z,'smoothingspline', 'SmoothingParam', 0.9 );
[Z.dz,Z.ddz] = differentiate(fz,time);


x = coeffvalues(fx);
y = coeffvalues(fy);
z = coeffvalues(fz);

X.x = x.coefs(1:end,4);     % here xyz trajectories has been substituted by their respective
Y.y = y.coefs(1:end,4);     % spline interpolation array in order to get rid of possible 
Z.z = z.coefs(1:end,4);     % presence of noise in the original waveforms

time = time(1:length(X.x));
T    = time(end) ; 

% Prepare trajectory structures to be fitted in DMP models
%--------------------------------------------------------------------------
X.x     = [ time, X.x   ];                             
X.dx    = [ time, X.dx  ];
X.ddx   = [ time, X.ddx ];
X.nx    = [ time, flip(X.x(1:end,2))   ];         % the 'n' prefix indicates that this is the desired trajectory 
X.ndx   = [ time, flip(-X.dx(1:end,2)) ];         % for the backward execution of the task. Note that, to achieve
X.nddx  = [ time, flip(X.ddx(1:end,2)) ];         % that, DMP's dynamics should be properly inverted

Y.y     = [ time, Y.y   ];
Y.dy    = [ time, Y.dy  ];
Y.ddy   = [ time, Y.ddy ];
Y.ny    = [ time, flip(Y.y(1:end,2))   ];
Y.ndy   = [ time, flip(-Y.dy(1:end,2)) ];
Y.nddy  = [ time, flip(Y.ddy(1:end,2)) ];

Z.z     = [ time, Z.z   ];
Z.dz    = [ time, Z.dz  ];
Z.ddz   = [ time, Z.ddz ];
Z.nz    = [ time, flip(Z.z(1:end,2))   ];
Z.ndz   = [ time, flip(-Z.dz(1:end,2)) ];
Z.nddz  = [ time, flip(Z.ddz(1:end,2)) ];

End = length(Y.ddy);

%%
% Calculate tangential array to the desired path from Frenet-Serret frame
%--------------------------------------------------------------------------

xyz = [ X.x(:,2) Y.y(:,2) Z.z(:,2) ]';             % Tangential direction of the curve is employed
[~,~,ttlist,~,~]=frenet_robust(xyz,50,0.1);        % just during simulations for kinesthetic guidance.
                                                   % This code snippet is not useful to fit DMP model,
X.tp = ttlist(1,:)';                               % it has been reported just for coherence with our project
Y.tp = ttlist(2,:)';
Z.tp = ttlist(3,:)';

X.tp = flip(X.tp);    
Y.tp = flip(Y.tp);
Z.tp = flip(Z.tp);


%% ------------------------------------------------------------------------
% 2 - Model Variables
%--------------------------------------------------------------------------

% Trasformation systems
%--------------------------------------------------------------------------

sys_x = DMP_model        ... 
        (                ...
          100          , ...  N 
          time(end)    , ...  T, phase0
          1            , ...  m 
          8            , ...  a_z
          2            , ...  b_z
          1            , ...  t0 
          Ts             ...  Ts    
        ) ;

sys_y = DMP_model        ... 
        (                ...
          100          , ...  N 
          time(end)    , ...  T, phase0
          1            , ...  m 
          8            , ...  a_z
          2            , ...  b_z
          1            , ...  t0 
          Ts             ...  Ts    
        ) ;

sys_z = DMP_model        ... 
        (                ...
          100          , ...  N 
          time(end)    , ...  T, phase0
          1            , ...  m 
          8            , ...  a_z
          2            , ...  b_z
          1            , ...  t0 
          Ts             ...  Ts    
        ) ;

sys_x.begin = X.x(1,2);
sys_x.goal  = X.x(end,2);
sys_y.begin = Y.y(1,2);
sys_y.goal  = Y.y(end,2);
sys_z.begin = Z.z(1,2);
sys_z.goal  = Z.z(end,2);


%% ------------------------------------------------------------------------
% 3 - Fit models with desired trajectories
%--------------------------------------------------------------------------    

% Estimate weights for forcing term interpolation
%--------------------------------------------------------------------------
clc

tic
sys_x.fit2 ( fd_x, X.x(1:end,2), time, 1 );                 % 'fd_' stands for the desired forcing term in forward direction (1);            
sys_y.fit2 ( fd_y, Y.y(1:end,2), time, 1 );                 % it is of size (lenght(X.x),1) and needs to be previously          
sys_z.fit2 ( fd_z, Z.z(1:end,2), time, 1 );                 % calculated.        

sys_x.fit2 ( fd_nx, X.nx(1:end,2), time, -1 );              % 'nfd_' stands for the desired forcing term in backward direction (-1)          
sys_y.fit2 ( fd_ny, Y.ny(1:end,2), time, -1 );              % it is of size (lenght(X.nx),1)          
sys_z.fit2 ( fd_nz, Z.nz(1:end,2), time, -1 ); 
toc

clear fd_x fd_y fd_z 
clear fd_nx fd_ny fd_nz


% Get forcing function as a 'sym' function 
%--------------------------------------------------------------------------
clc

tic 
sys_x.get_forcing_function( 1 );                      
sys_y.get_forcing_function( 1 );
sys_z.get_forcing_function( 1 );

sys_x.get_forcing_function( -1 );
sys_y.get_forcing_function( -1 );
sys_z.get_forcing_function( -1 );
toc


% Save desired forcing term in order to export it in Simulink simulations
%--------------------------------------------------------------------------
clc

tic
matlabFunction ( sys_x.F  , 'File', '~/path/to/desired/folder/x_forcing_fn');
matlabFunction ( sys_y.F  , 'File', '~/path/to/desired/folder/y_forcing_fn');
matlabFunction ( sys_z.F  , 'File', '~/path/to/desired/folder/z_forcing_fn');

matlabFunction ( sys_x.nF  , 'File', '~/path/to/desired/folder/nx_forcing_fn');
matlabFunction ( sys_y.nF  , 'File', '~/path/to/desired/folder/ny_forcing_fn');
matlabFunction ( sys_z.nF  , 'File', '~/path/to/desired/folder/nz_forcing_fn');
toc

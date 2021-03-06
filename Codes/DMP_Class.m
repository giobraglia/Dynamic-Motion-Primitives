%  This Class enables to create DMP object for cartesian discrete motion. Here in this project we used a linear Canonical
%  system to account for DMP phase and it was implemented indipendently, therefore this class enables just to fit a 
%  trasformation system, providing in particular its forcing function. 
%  Here methods are implemented using parfor for optimizing codes through parallel computation.
%  Two fit methods here are presented based on different Transformation System representation:
%  - fit1 -->  A. J. Ijspeert, J. Nakanishi, H. Hoffmann, P. Pastor and S. Schaal,
%              "Dynamical Movement Primitives: Learning Attractor Models for Motor Behaviors" 
%  - fit2 -->  A. J. Ijspeert, J. Nakanishi and S. Schaal, 
%              "Movement imitation with nonlinear dynamical systems in humanoid robots" 
%  
%  Please refer to DMP_example.mat for a correct implementation of the DMP_model



%--------------------------------------------------------------------------
classdef DMP_model < handle
%--------------------------------------------------------------------------
  

properties  ( GetAccess = public, SetAccess = public )
%--------------------------------------------------------------------------
        N                                                                       % number of radial basis functions (RBF) used for interpolation
        T                                                                       % trajectory length
        m                                                                       % slope of linear Canonical System
        a_z                                                                     %   parameters of
        b_z                                                                     %   Transformation System
        phase0                                                                  % initial value of the phase
        Ts                                                                      % sampling period
        var                                                                     % variance of RBFs
        c                                                                       % centers of RBFs
        begin                                                                   % initial position
        goal                                                                    % final position
end


properties  ( GetAccess = public, SetAccess = private )
%--------------------------------------------------------------------------   
        
        sd                                                                      % desired position trajectory
        fd                                                                      % desired  forcing function (FORWARD)
        nfd                                                                     % desired  forcing function (BACKWARD)
        F                                                                       % computed forcing function (FORWARD)
        nF                                                                      % computed forcing function (BACKWARD)
        w                                                                       % weights (FORWARD)
        nw                                                                      % weights (BACKWARD)
        rbf
        time
        phase                                                                   % phase (FORWARD)    
        nphase                                                                  % phase (BACKWARD)
        x
end


methods  ( Access = private )
%--------------------------------------------------------------------------

    function  get_rbfs ( obj )                                                 % creates a pool of 'N' RBFs 
    %----------------------------------------------------------------------    % centered in 'c' and with variance 'var'
        RBF = cell ( [obj.N, 1] );
        C   = obj.c;
        k   = obj.var;
        n   = obj.N;
        X   = obj.x;

        parfor i = 1:n
            RBF {i} = exp ( - (X - C(i)).^2 / ( 2 * k^2 ) )';
        end   

        obj.rbf = RBF;
    end 
    

    function  get_kernels ( obj )                                             % creates linearly spaced centers 'c'  
    %----------------------------------------------------------------------   % for 'N' RBFs
        obj.c = linspace ( 0, obj.T, obj.N );
        obj.get_rbfs();
    end

   
    function W = get_weights ( obj, phase, fd )                               % calculates weights 'w'
    %----------------------------------------------------------------------   % how weight are calculated is explained 
        W   = zeros ( [obj.N,1] );                                            % on the article relative to
        s   = ( phase );                                                      % fit1 method
        RBF = obj.rbf;
        FD  = fd;
        n   = obj.N;

        parfor i = 1 : n
            rbfFH = matlabFunction(RBF{i});
            psi   = rbfFH(s);
            gamma = diag ( psi );
            
            W (i) =     ( s' * gamma * FD ) / ...
                   ... -------------------------     
                           ( s' * gamma * s );
            if ( isnan ( W(i) ) ) 
                W(i) = 0;
            end
        end
        
%         obj.w = W;
    end



end


methods  ( Access = public )
%--------------------------------------------------------------------------

    function obj = DMP_model ( N, T, m, a_z, b_z, a_py, t0, Ts )
    %----------------------------------------------------------------------
        if ( nargin == 8 )
            syms x
            obj.N       = N         ;   
            obj.T       = T         ; 
            obj.m       = m         ;      
            obj.a_z     = a_z       ;
            obj.b_z     = b_z       ;
            obj.a_py    = a_py      ;
            obj.phase0  = T         ;
            obj.t0      = t0        ;
            obj.Ts      = Ts        ;
            obj.var     = 0.2       ;
            obj.x       = x         ;
        else
            disp ('Wrong number of parameters for DMP model, 9 required');
        end
    end

    
    function fit1 ( obj, s, ds, dds, time )                                   % calculates weights 'w' using method fit1
    %----------------------------------------------------------------------
        obj.phase = reshape ( ( obj.phase0 : -obj.Ts : 0 ), [], 1 );
        obj.time  = time;
        
        obj.sd    = s;
        obj.begin = s(1);
        obj.goal  = s(end);

        obj.fd = 1 / ( obj.a_z * obj.b_z) *         ...
                 (                                  ...
                   dds - obj.a_z *                  ...
                   (                                ...
                     obj.b_z *                      ...
                     (                              ...
                       obj.goal - s -               ...
                       ( obj.goal - obj.begin ) *   ...
                       obj.phase / obj.phase0       ...
                     )                              ...
                     - ds                           ...
                   )                                ...
                 );
        
        obj.fd = [ obj.time, obj.fd ]; 
        obj.get_kernels();
        obj.get_weights();
    end

    
    function fit2 ( obj, fd, s, time, dir )                                 % calculates weights 'w' using method fit2
    %----------------------------------------------------------------------   
        if ( dir == 1 )
            obj.phase = reshape ( ( obj.phase0 : -obj.Ts : 0 ), [], 1 );
            obj.time  = time;

            obj.sd    = s;
            obj.begin = s(1);
            obj.goal  = s(end);

%             obj.fd = [ obj.time, fd ]; 
            obj.fd =  fd ;
            obj.get_kernels();
            obj.w  = obj.get_weights( obj.phase, obj.fd );
        end
        
        if ( dir == -1 )
            obj.nphase = flip( obj.phase );
            obj.sd     = s;
            obj.nfd    =  fd; 
            obj.get_kernels();
            obj.nw = obj.get_weights( obj.nphase, obj.nfd - fd(1) );
        end
    end
    
    
    function get_fit_error ( obj, F, phase, fd )                                 % calculates the mean square error between
    %----------------------------------------------------------------------      % the calculated forcing function and the
        simf = matlabFunction(F);                                                % desired one
        simf = simf( phase );
        simf = reshape ( simf, [], 1 );

        e   = sum ( abs(simf - fd).^2 );
        disp (['Fit error:', num2str(e) ]);
    end

   
    function get_forcing_function ( obj, dir )                                  % creates the non-linear forcing term as a
    %----------------------------------------------------------------------     % 'sym' variable
        rbf_sum    = 0;
        rbf_sum_w  = 0;
        
        RBF = obj.rbf;
        n   = obj.N;
        
        if ( dir == 1  ) 
            W  = obj.w;
        end
        if ( dir == -1 )
            W = obj.nw;
        end

        parfor i = 1:n
            rbf_sum   = rbf_sum + RBF{i};
            rbf_sum_w = rbf_sum_w + W(i) * RBF{i};
        end
         
        if ( dir == 1 )
            obj.F = (rbf_sum_w / rbf_sum * obj.x ); 
            obj.get_fit_error( obj.F, obj.phase, obj.fd );
        end
        if ( dir == -1 )
            obj.nF = (rbf_sum_w / rbf_sum * obj.x ) + obj.nfd(1); 
            obj.get_fit_error( obj.nF, obj.nphase, obj.nfd );
        end
    end



end                                                                          % end methods



end                                                                          % end DMP_model



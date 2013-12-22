% grnd_trjs = ground_truth( mdl, user )
%
% Analytically compute trajectories
%
% Inputs:
% mdl  - struct - the hybrid model
% user - struct - holds all the simulation parameters
%
% Outputs:
% grnd_trjs - list of trajectory structs
% grnd_trj - trajectory struct
%   .ti   - begining time
%   .tf   - final time
%   .sol  - the analytically computed solution
%   .dsol - the analytically computed derivative of the solution
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function grnd_trjs = ground_truth( mdl, user )

% parameters used repeatedly
time_arr = user.time_arr;
dr = user.params( user.mdl.idp.damping_ratio );
nf = user.params( user.mdl.idp.natural_freq );
cr = user.params( user.mdl.idp.coefficient_restitution );
x_max = user.params( user.mdl.idp.x_max );
input_amp = user.input_amp;
input_freq = user.input_freq;
w_tilde = sqrt( nf^2 - dr^2 );
eta = dr/w_tilde;

% deciding search space
nVel = -1e-5;
nStep = 1e-10;

% create a dummy time variable
syms tau t s u A B;

% construct the particular solution
% xp = ( 0.5/w_tilde ) * input_amp * exp( dr * ( u - t ) ) * ...
%     ( ( dr * sin( w_tilde * ( t - u ) + input_freq * u ) + ...
%     ( dr - input_freq ) * cos( w_tilde *  ( t - u ) + input_freq * u ) )/...
%     ( dr^2 + ( w_tilde - input_freq )^2 ) + ...
%     ( dr * sin( w_tilde * t - u * ( w_tilde + input_freq ) ) + ( w_tilde + input_freq ) + ...
%     cos( w_tilde * t - u * ( w_tilde + input_freq ) ) )/...
%     ( dr^2 + ( w_tilde + input_freq )^2 ));
% xp = subs( xp, u, 0 ) - subs( xp, u, t );

xp = int( input_amp * cos( input_freq * tau ) * ...
    exp( -dr * ( t - tau ) ) * sin( w_tilde * ( t - tau ) ), ...
    tau, 0, t ) / w_tilde;
% xpdot = int( input_amp * cos( input_freq * tau ) * ...
%     exp( -nf * ( t - tau ) ) * ( cos( w_tilde * ( t - tau ) ) - ...
%     eta * sin( w_tilde * ( t - tau ) ) ), ...
%     tau, 0, t ) / w_tilde;

% construct the actual solution
x = exp( -dr * t ) * ( A * cos( w_tilde * t ) + B * sin( w_tilde * t ) ) + xp;
xdot = diff( x, 't' );

% construct the trajectory going forward
trj_idx = 1;
grnd_trjs( trj_idx ).ti = time_arr( 1 );
% used for recursive computation (equation 6)
coeff_1 = user.x0( user.mdl.ids0.x );
coeff_2 = user.x0( user.mdl.ids0.xdot )/w_tilde + eta * user.x0( user.mdl.ids0.x );
grnd_trjs( trj_idx ).sol = matlabFunction( subs( subs( x, A, coeff_1 ), B, coeff_2 ) );
grnd_trjs( trj_idx ).dsol = matlabFunction( subs( subs( xdot, A, coeff_1 ), B, coeff_2 ) );

while( grnd_trjs( trj_idx ).ti < time_arr( 2 ) )
    
    % are we zenoing ...
    if( grnd_trjs( trj_idx ).sol( grnd_trjs( trj_idx ).ti ) >= x_max  && ... 
            grnd_trjs( trj_idx ).dsol( grnd_trjs( trj_idx ). ti ) >= nVel && ...
            input_amp * cos( input_freq * ( grnd_trjs( trj_idx ).ti + nStep ) ) >= nf^2 * x_max )
        
        % the solution is trivial during the zeno part of the execution
        grnd_trjs( trj_idx ).sol = @(y) x_max;
        grnd_trjs( trj_idx ).dsol = @(y) 0;
        
        helper = @(y) input_amp * cos( input_freq * y ) - nf^2 * x_max;
        
        search_interval = linspace( grnd_trjs( trj_idx ).ti + nStep, time_arr( 2 ), 10000 );
        search_value = [];
        for i = search_interval
            search_value = [ search_value helper( i ) ];
        end
        
        if( sign( helper( grnd_trjs( trj_idx ).ti + nStep ) )  == 1 )
            t_max = search_interval( find( search_value < 0, 1 ) );
        else
            t_max = search_interval( find( search_value > 0, 1 ) );
        end
        
        if( isempty( t_max ) )
            grnd_trjs( trj_idx ).tf= time_arr( 2 );
            break;
        end
        
        tf = fzero( helper, [ grnd_trjs( trj_idx ).ti + nStep t_max ] );
        grnd_trjs( trj_idx ).tf = tf;
        
        % determine coeff_1 and coeff_2  using the initial conditions to
        % the ode
        helper1 = subs( x, t, tf ) - x_max;
        helper2 = subs( xdot, t, tf );
        sol_coeff = solve( helper1, helper2 );
        
        coeff_1 = sol_coeff.A;
        coeff_2 = sol_coeff.B;
        
        trj_idx = trj_idx + 1;
        grnd_trjs( trj_idx ).ti = tf;
        
        % perform substitution and turn them into matlabFunctions in order to
        % speed up computation
        grnd_trjs( trj_idx ).sol = matlabFunction( subs( subs( x, A, coeff_1 ), B, coeff_2 ) );
        grnd_trjs( trj_idx ).dsol = matlabFunction( subs( subs( xdot, A, coeff_1), B, coeff_2 ) );
    else
        helper = @(y)  x_max - grnd_trjs( trj_idx ).sol( y );
        
        search_interval = linspace( grnd_trjs( trj_idx ).ti + nStep, time_arr( 2 ), 10000 );
        search_value = [];
        for i = search_interval
            search_value = [ search_value helper( i ) ];
        end
        
        if( sign( helper( grnd_trjs( trj_idx ).ti + nStep ) )  == 1 )
            t_max = search_interval( find( search_value < 0, 1 ) );
        else
            t_max = search_interval( find( search_value > 0, 1 ) );
        end
        
        if( isempty( t_max ) )
            grnd_trjs( trj_idx ).tf= time_arr( 2 );
            break;
        end
        
        tf = fzero( helper, [ grnd_trjs( trj_idx ).ti + nStep t_max ] );
        grnd_trjs( trj_idx ).tf = tf;
        
        xdot_tf = grnd_trjs( trj_idx ).dsol( tf );
        coeff_1 = coeff_1 + ( 1 + cr )/w_tilde * exp( dr * tf ) * sin( w_tilde * tf ) * xdot_tf;
        coeff_2 = coeff_2 - ( 1 + cr )/w_tilde * exp( dr * tf ) * cos( w_tilde * tf ) * xdot_tf;
        
        trj_idx = trj_idx + 1;
        grnd_trjs( trj_idx ).ti = tf;
        
        % perform substitution and turn them into matlabFunctions in order to
        % speed up computation
        grnd_trjs( trj_idx ).sol = matlabFunction( subs( subs( x, A, coeff_1 ), B, coeff_2 ) );
        grnd_trjs( trj_idx ).dsol = matlabFunction( subs( subs( xdot, A, coeff_1), B, coeff_2 ) );
    end
    
end












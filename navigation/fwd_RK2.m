% trjs = fwd_RK2( mdl, user, h_in, rx )
%
% Compute Runge-Kutta approximation to relaxed hybrid execution
%
% Inputs:
% mdl  - struct - the hybrid model
% user - struct - holds all the simulation parameters
% h_in - scalar - step size (optional argument)
% rx   - scalar - relaxation size (optional argument)
%
% Outputs:
% trjs - list of trajectory structs
% trj - trajectory struct
%   .t - times
%   .x - states
%   .u - inputs
%   .p - parameters
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function trjs = fwd_RK2( mdl, user, h_in, rx )

% if optional arguments were not set
if ( nargin < 3 )
    h_in = user.step_size;
    rx = user.rx;
end

% simulation parameters
t_in = user.time_arr;
x_in = user.x0;
p_in = user.params;
n = user.max_jumps;
beta = user.euler_beta; % beta = 0.5
kmax = user.euler_kmax; % kmax = 50
debug = user.debug;

% parameter indices
idp = mdl.idp;

% store final time
tf = t_in(2);

% store given stepsize
h = h_in;

% arrays that are continually updated
t = t_in(1);
x = x_in;
params = p_in;

% hybrid system functions
vfield = mdl.F;
guard = mdl.G;
reset = mdl.R;

% Euler loop
idx = 1;
trj_idx = 1;
while( ( t( idx ) < tf ) && ( trj_idx <= n ) && ... 
        params( user.mdl.idp.mode ) < 8  && params( user.mdl.idp.mode ) > -1 ...
        && ( x( idx, 1 ) < user.setx0( 1, 1 ) || x( idx, 1 ) > user.setx0( 1, 2 ) || ...
        x( idx, 2 ) < user.setx0( 2, 1 ) || x( idx, 2 ) > user.setx0( 2, 2 ) || ...
        x( idx, 3 ) < user.setv0( 1, 1 ) || x( idx, 3 ) > user.setv0( 1, 2 ) || ...
        x( idx, 4 ) < user.setv0( 2, 1 ) || x( idx, 4 ) > user.setv0( 2, 2 ) || trj_idx ==  1 ) )
    % check if we are in an actual mode
    if( isnan( params( user.mdl.idp.mode ) ) )
        warning( 'ex:fwd_RK2', '(fwd_RK2) Current mode is NaN!' );
        break;
    end
    
    % do RK2 step (RV: this is a complete hack, but just done to even the
    % playing field)
    hdx1 = h * vfield( x( idx, : ), params, user );
    hdx2 = h * vfield( x( idx, : ) + hdx1/2, params, user );
    t( idx + 1 ) = t( idx ) + h;
    x( idx + 1, : ) = x( idx, : ) + hdx2;
    g = guard( x( idx + 1, : ), params, user );
    
    % halve step size until trajectory doesn't jump over strip
    k = 0;
    while( any( g < -rx ) )
        dx = vfield( x( idx, : ), params, user );
        if( debug )
            fprintf( 1, 'RX:  jump over strip #%d\n', k );
        end
        assert( k < kmax, '(euler) strip iterations exceeded' );
        h = beta * h;
        t( idx + 1 ) = t( idx ) + h;
        x( idx + 1, : ) = x( idx, : ) + h * dx;
        g = guard( x( idx + 1, : ), params, user );
        k = k + 1;
    end
    
    if( debug )
        %fprintf( 1, '  :  j = %d t = %0.3f g = %0.3f s = %0.3f\n', params( idp.mode ), t(idx+1), g, s(idx,:) );
        js  = sprintf( '%d ', params( idp.mode ) );
        ts  = sprintf( '%0.1f ', t( idx ) );
        xs  = sprintf( '%0.1f ', x( idx, : ) );
        dxs = sprintf( '%0.1f ', dx);
        gs  = sprintf( '%0.1f ', g);
        disp( [ '  :  j = ', js, ' t = ', ts, ' g = ', gs ] );
        disp( [ '     x = ', xs ] );
        disp( [ '    dx = ', dxs ] );
    end
    
    % if state is on strip
    if( any( g < 0 ) )
        
        % spend time on the strip
        t( idx + 2 ) = t( idx + 1 ) + ( rx + min( g ) );
        x( idx + 2, : ) = x( idx + 1, : );
        
        if( debug )
            fprintf( 1, 'RX:  j = %d t = %0.3f \n', params( idp.mode ), t( idx + 2 ), g( 1 ) );
        end
        
        % append trj to trjs
        trjs( trj_idx ).t = t( 1:( idx + 2 ) );
        trjs( trj_idx ).x = x( 1:( idx + 2 ), : );
        trjs( trj_idx ).p = params;
        
        trj_idx = trj_idx + 1;
        
        % apply reset to modify trj
        [ x, params ] = reset( x( idx + 2, : ), params, user );
        t = t( idx + 2 );
        
        % re-initialize step size
        h = h_in;
        
        if( debug )
            g = guard( x( 1, : ), params, user );
            fprintf( 1, 'RX:  j = %d t = %0.3f g = %0.3f\n', params( idp.mode ), t( 1 ), g( 1 ) );
            js  = sprintf( '%d ', params( idp.mode ) );
            ts  = sprintf( '%0.1f ', t );
            xs  = sprintf( '%0.1f ', x );
            gs  = sprintf( '%0.1f ', g );
            disp( [ 'g :  j = ', js, ' t = ', ts, ' g = ', gs ] );
            disp( [ '     x = ', xs ] );
        end
        idx = 1;
    else
        idx = idx + 1;
    end
end % while

if( ( ~isnan( params( user.mdl.idp.mode ) ) ) )
    trjs( trj_idx ).t = t( 1:idx );
    trjs( trj_idx ).x = x( 1:idx, : );
    trjs( trj_idx ).p = params;
end

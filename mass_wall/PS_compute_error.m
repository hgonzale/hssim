% max_error = PS_compute_error( grnd_trjs, ps_trjs, user )
%
% Compute the error between trajectories using the \rho distance
%
% Inputs:
% grnd_trjs  - struct - ground truth trajectory data
% ps_trjs    - struct - PS computed trajectory
% user       - struct - holds all the simulation parameters
%
% Outputs:
% max_error  - scalar - corresponding to the maximum error
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function max_error = PS_compute_error( grnd_trjs, ps_trjs, user )

% initialize max_error
max_error = 0;

% setting up structures
t = ps_trjs.t;
x = ps_trjs.x;

% iterate through all trajectory data
for i = 1:length( t )
    % don't need to compute for times beyond ground truth values
    if ( t( i ) > user.time_arr( 2 ) )
        break;
    end
    index = 0;
    for j = 1:length( grnd_trjs )
        if ( grnd_trjs( j ).ti <= t( i ) && grnd_trjs( j ).tf >= t( i ) )
            index = j;
            break;
        end
    end
    
    error_i = abs( x( i, user.mdl.ids0.x ) -  grnd_trjs( index ).sol( t( i ) ) ); 
    max_error = max( [ max_error error_i ] );
end
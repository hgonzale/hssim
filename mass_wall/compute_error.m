% max_error = compute_error( grnd_trjs, trjs, user )
%
% Compute the error between trajectories using the \rho distance
%
% Inputs:
% grnd_trjs  - struct - ground truth trajectory data
% trjs       - struct - euler computed trajectory
% user       - struct - holds all the simulation parameters
% rx         - scalar - relaxation size
%
% Outputs:
% max_error  - scalar - corresponding to the maximum position error
% rho        - scalar - corresponding to the rho as in Equation (19) of our paper
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function [ max_error, rho ] = compute_error( grnd_trjs, trjs, user, rx )

% parameters
x_max = user.params( user.mdl.idp.x_max );
cr = user.params( user.mdl.idp.coefficient_restitution );

% initialize max_error
max_error = 0;
rho = 0;

% concatenate all of the trajectories to make analysis easier (too memory
% intensive)
% [ t, x, ~, ~ ] = stack( trjs );


% iterate through all trajectory data
for k = 1:length( trjs )
    for i = 1:( length( trjs( k ).t ) - 1 )
        % don't need to compute for times beyond ground truth values
        if ( trjs( k ).t( i + 1 ) > user.time_arr( 2 ) )
            break;
        end
        index = 0;
        for j = 1:length( grnd_trjs )
            if ( grnd_trjs( j ).ti <= trjs( k ).t( i ) && grnd_trjs( j ).tf >= trjs( k ).t( i ) )
                index = j;
                break;
            end
        end
        
        error_i = abs( trjs( k ).x( i, user.mdl.ids0.x ) -  grnd_trjs( index ).sol( trjs( k ).t( i ) ) );
        
        
%         error_i = abs( trjs( k ).x( i, user.mdl.ids0.x ) -  grnd_trjs( index ).sol( trjs( k ).t( i ) ) ) + ...
%             abs( trjs( k ).x( i, user.mdl.ids0.xdot ) -  grnd_trjs( index ).dsol( trjs( k ).t( i ) ) );
%         
%         % code snippet to determine if using the quotient distance across
%         % the strip reduces the distance
%         qerror_i = abs( trjs( k ).x( i, user.mdl.ids0.x ) - x_max - rx ) + abs( grnd_trjs( index ).sol( trjs( k ).t( i ) ) - x_max - rx );
%         
%         % using the quotient distance is worth it
%         if( qerror_i < error_i )
%             vdist = compute_vdiff( grnd_trjs( index ).dsol( trjs( k ).t( i ) ), trjs( k ).x( i, user.mdl.ids0.xdot ), cr );
%             error_i = min( [ error_i; qerror_i + vdist ] );
%         end
        max_error = max( [ max_error error_i ] );
        
        
        helper = @(y) compute_diff( y, grnd_trjs, trjs( k ).x( i, user.mdl.ids0.x ), trjs( k ).x( i + 1, user.mdl.ids0.x ), ...
            trjs( k ).x( i, user.mdl.ids0.xdot ), trjs( k ).x( i + 1, user.mdl.ids0.xdot ), trjs( k ).t( i ), trjs( k ).t( i + 1 ), x_max, rx, cr );
        
        blah = fminbnd( helper, trjs( k ).t( i ), trjs( k ).t( i + 1 ) );
        rho_i = -helper( blah ); % since we are trying to maximize the difference
        rho = max( [ rho rho_i ] );
    end
end

function diff = compute_diff( t, grnd_trjs, x0, x1, v1, v0, t0, t1, x_max, rx, cr )

index = 0;
for j = 1:length( grnd_trjs )
    if ( grnd_trjs( j ).ti <= t && grnd_trjs( j ).tf >= t )
        index = j;
        break;
    end
end

tx = ( x0 + ( x1 - x0 ) * ( t - t0 )/( t1 - t0  ) );
tv = ( v0 + ( v1 - v0 ) * ( t - t0 )/( t1 - t0  ) );

diff = abs( grnd_trjs( index ).sol( t ) - tx ) + ...
    abs( grnd_trjs( index ).dsol( t ) - tv ) ;

qdiff = abs( grnd_trjs( index ).sol( t ) - x_max - rx ) + abs( tx - x_max - rx );

if ( qdiff < diff )
    vdist = compute_vdiff( grnd_trjs( index ).dsol( t ), tv, cr );
    diff = min( [ diff; qdiff + vdist ] );
end

% the minus sign appears since we are trying to maximize
diff = -diff;

% compute distance across the quotient-ing of the strip
function vdiff = compute_vdiff( gv, tv, cr )

vdiff = abs( tv -  gv );
if ( abs( tv ) < abs( gv ) )
    qv1 = gv * -cr;
    v2 = tv;
elseif ( abs( gv ) < abs( tv ) )
    qv1 = tv * -cr;
    v2 = gv;
end
qqv1 = qv1 * -cr;
qdist = min( [ abs( qv1 - v2 ); abs( qqv1 - v2 ) ] );

while( qdist < vdiff && abs( qdist - vdiff ) > eps )
    vdiff = qdist;
    qv1 = qqv1 * -cr;
    qqv1 = qv1 * -cr;
    qdist = min( [ abs( qv1 - v2 ); abs( qqv1 - v2 ) ] );
end




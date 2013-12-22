% mw = forced_mass_wall_mdl( p_file, debug );
%
% Forced Mass Wall hybrid dynamical model
%
% Inputs:
% p_file - string - filename for the parameters in the model
% debug - bool - flag for printing debugging info
%
% Outputs:
% mdl - struct - simulation object
%   .U - continuous input u     = U( T, P, user )
%   .F - vector field     dX    = F( T, X, U, P, user )
%   .G - guard            g     = G( X, P, user )
%   .R - reset map        [x,p] = R( X, P, U, user )
%   .O - observations     o     = O( X, U, P, user )
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function mdl = forced_mass_wall_mdl( p_file , debug )

if( nargin == 0 )
    debug = 0;
end

mdl = read_config( p_file );
mdl.U = @U; mdl.F = @F;
mdl.G = @G; mdl.R = @R;
mdl.O = @O; mdl.debug = debug;

end % dpmdl



% input
function u = U( t, P, user )

amp = user.input_amp;
freq = user.input_freq;

u = amp * cos( freq * t );

end % input

% guard function
% g > 0 : guard inactive
% g = 0 : guard
function g = G( X, P, user )

switch( P( user.mdl.idp.mode ) )
    case 0
        % is position beyond the wall?
%         g = P( user.mdl.idp.x_max ) - X( user.mdl.ids0.x );
        g = max( [P( user.mdl.idp.x_max ) - X( user.mdl.ids0.x ) -X( user.mdl.ids0.xdot ) ] );
    otherwise
        fprintf( 1, '(G) unknown discrete mode\n' );
        g = NaN;
end % switch

end % guard function



% vector field
function dX = F( X, U, P, user )

% extract parameters
dr = P( user.mdl.idp.damping_ratio );
nf = P( user.mdl.idp.natural_freq );

switch( P( user.mdl.idp.mode ) )
    case 0
        % is position beyond the wall?
        dX = zeros( size( user.mdl.ids0.states ) );
        dX( user.mdl.ids0.x ) = X( user.mdl.ids0.xdot );
        dX( user.mdl.ids0.xdot ) = U - 2 * dr * X( user.mdl.ids0.xdot ) - nf^2 * X( user.mdl.ids0.x );
    otherwise
        error( '(F) unknown discrete mode' );
end % switch

end % vector field



% reset map
function [ x, p ] = R( X, P, user )

% extract parameters
cr = P( user.mdl.idp.coefficient_restitution );

% guard values
g = G( X, P, user );

switch( P( user.mdl.idp.mode ) )
    case 0
        if g < 0
            p = P; % the parameters always stay the same.
            x = zeros( size( user.mdl.ids0.states ) );
            x( user.mdl.ids0.x ) = X( user.mdl.ids0.x );
            x( user.mdl.ids0.xdot ) = -cr * X( user.mdl.ids0.xdot );
        else
            error( '(R) no applicable reset map' );
        end
    otherwise
        error( '(R) unknown discrete mode' );
end % switch

end % reset map



% observations
function o = O( X, U, P, user )

switch( P( user.mdl.idp.mode ) )
    case 0
        o = [ X( user.mdl.ids0.x ); X( user.mdl.ids0.xdot ); U ];
    otherwise
        error( '(O) unknown discrete mode' );
end % switch

end % observations

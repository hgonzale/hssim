% mw = nav_mdl( p_file, debug );
%
% Navigation hybrid dynamical model
%
% Inputs:
% p_file - string - filename for the parameters in the model
% debug - bool - flag for printing debugging info
%
% Outputs:
% mdl - struct - simulation object
%   .F - vector field     dX    = F( T, X, U, P, user )
%   .G - guard            g     = G( X, P, user )
%   .R - reset map        [x,p] = R( X, P, U, user )
%   .O - observations     o     = O( X, U, P, user )
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function mdl = nav_mdl( p_file , debug )

if( nargin == 0 )
    debug = 0;
end

mdl = read_config( p_file );
mdl.F = @F;
mdl.G = @G; mdl.R = @R;
mdl.O = @O; mdl.debug = debug;

end % dpmdl


% guard function
% g > 0 : guard inactive
% g = 0 : guard
function g = G( X, P, user )

% extract parameters
map = user.map;

% relevant params (the indices are the same for all states)
x = X( user.mdl.ids0.x );
y = X( user.mdl.ids0.y );
xdot = X( user.mdl.ids0.xdot );
ydot = X( user.mdl.ids0.ydot );
xmax = size( map, 2 );
ymax = size( map, 1 );

% one guard corresponding to each side of the cell
g = zeros( 4, 1 );

% have to do some calculations to get the real position instead of the
% matrix one
actual_locx = P( user.mdl.idp.locx );
actual_locy = size( user.map, 1 ) + 1 - P( user.mdl.idp.locy );

% strictly these guards are not 100% correct (due to the -x and -y terms),
% but as long as rx < 1 this will work just fine.
% west
g( 1 ) = max( [ ( x - ( actual_locx - 1 ) ) -x  xdot  ] ); 
% ( x <= ( actual_locx - 1 ) ) && ( x > 0 ) && ( xdot < 0 );
% north
g( 2 ) =  max( [ ( actual_locy - y ) ( y - ymax ) -ydot ] );
%( y >= actual_locy ) && ( y < ymax ) && ( ydot > 0 );
% east
g( 3 ) =  max( [ ( actual_locx - x ) ( x - xmax ) -xdot ] );
%( x >= actual_locx ) && ( x < xmax ) && ( xdot > 0 );
% south
g( 4 ) = max( [ ( y - ( actual_locy - 1 ) ) -y ydot ] );
%( y <= ( actual_locy - 1 ) ) && ( y > 0 ) && ( ydot < 0 );

end % guard function



% vector field
function dX = F( X, P, user )

% input params
A = user.A;

% relevant params (the indices are the same for all states)
v = X( user.mdl.ids0.xdot:user.mdl.ids0.ydot );

% desired velocity ( depends on the hybrid mode )
vd = [ sin( P( user.mdl.idp.mode ) * pi/4 ) cos( P( user.mdl.idp.mode ) * pi/4 ) ];

% same number of continuous states in all discrete modes
dX = zeros( size( user.mdl.ids0.states ) );
dX( user.mdl.ids0.x ) = X( user.mdl.ids0.xdot );
dX( user.mdl.ids0.y ) = X( user.mdl.ids0.ydot );
dX( user.mdl.ids0.xdot:user.mdl.ids0.ydot ) = ( A * ( v - vd )' )';

end % vector field


% reset map
function [ x, p ] = R( X, P, user )

% input params
map = user.map;

% output parameters
p = zeros( size( P ) ); 
x = X; % no state reset

% guard values
g = G( X, P, user );

% deciding which location in the map we switch to
if ( ( g( 1 ) <= 0 ) && ( g( 2 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) - 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) - 1;
elseif ( ( g( 2 ) <= 0 ) && ( g( 3 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) + 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) - 1;
elseif( ( g( 3 ) <= 0 ) && ( g( 4 ) <= 0 ) ) 
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) + 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) + 1;
elseif( ( g( 4 ) <= 0 ) && ( g( 1 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) - 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) + 1;
elseif( ( g( 1 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) - 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy );
elseif( ( g( 2 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx );
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) - 1;
elseif( ( g( 3 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx ) + 1;
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy );
elseif( ( g( 4 ) <= 0 ) )
    p( user.mdl.idp.locx ) = P( user.mdl.idp.locx );
    p( user.mdl.idp.locy ) = P( user.mdl.idp.locy ) + 1;
else 
    p = P;
end

% deciding which mode we are in
p( user.mdl.idp.mode ) = map( p( user.mdl.idp.locy ), p( user.mdl.idp.locx ) );

end % reset map


% observations
function o = O( X, P, user )

o = [ X( user.mdl.ids0.x ); X( user.mdl.ids0.y ); X( user.mdl.ids0.xdot ); X( user.mdl.ids0.ydot ) ];

end % observations

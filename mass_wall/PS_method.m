% ps_trjs = PS_method( mdl, user, h_in )
%
% Analytically compute trajectories
%
% Inputs:
% mdl  - struct - the hybrid model
% user - struct - holds all the simulation parameters
% h_in - scalar - step size (optional argument)
%
% Outputs:
% ps_trjs - list of trajectory structs
% ps_trjs - trajectory struct
%   .t - times
%   .x - just the position 
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function ps_trjs = PS_method( mdl, user, h_in )

% if optional arguments were not set
if ( nargin < 3 )
    h_in = user.step_size;
end

% simulation parameters
t_in = user.time_arr;
x_in = user.x0;
p_in = user.params;
input = mdl.U;
dr = user.params( user.mdl.idp.damping_ratio );
nf = user.params( user.mdl.idp.natural_freq );
cr = user.params( user.mdl.idp.coefficient_restitution );
x_max = user.params( user.mdl.idp.x_max );

% state indices
xidx = user.mdl.ids0.x;
xdidx = user.mdl.ids0.xdot;

% done the 'stupid-way' in order to make the times more comparable
ps_trjs.t = t_in( 1 );
ps_trjs.t = [ ps_trjs.t; t_in( 1 ) + h_in ];
ps_trjs.x = x_in( xidx );
ps_trjs.x = [ ps_trjs.x; x_in( xidx ) + x_in( xdidx ) * h_in + ...
    + h_in^2/2 * ( input( t_in( 1 ), p_in, user ) - 2 * dr * x_in( xdidx ) - nf^2 * x_in( xidx ) ) ];

while( ps_trjs.t( end ) <= t_in( 2 ) )
    ps_trjs.x = [ ps_trjs.x; two_step( ps_trjs.t( end ), ps_trjs.x( end - 1 ), ps_trjs.x( end ), h_in, user ) ];
    ps_trjs.t = [ ps_trjs.t; ps_trjs.t( end ) + h_in ];
end

% % setting up data structures 
% ps_trjs.t =  ( t_in( 1 ):h_in:t_in( 2 ) )';
% if ( ps_trjs.t( end ) < t_in( 2 ) )
%     ps_trjs.t = [ ps_trjs.t; t_in( 2 ) ];
% end
% ps_trjs.x = zeros( length( ps_trjs.t ), 1 );
% 
% % actually performing computation
% ps_trjs.x( 1 ) = x_in( xidx );
% ps_trjs.x( 2 ) = x_in( xidx ) + x_in( xdidx ) * h_in + ...
%     + h_in^2/2 * ( input( t_in( 1 ), p_in, user ) - 2 * dr * x_in( xdidx ) - nf^2 * x_in( xidx ) );
% for i = 2:( length( ps_trjs.t ) - 1 )
% %    helper = @(y) y + cr * ps_trjs.x( i - 1 ) - min( [ ( 1 + cr ) * x_max ...
% %        ( h_in^2 * input( ps_trjs.t( i ), p_in, user ) + ( 2 - h_in^2 * nf^2 ) * ps_trjs.x( i ) ...
% %        - ( ( 1 - cr ) - ( 1 + cr ) * dr * h_in ) * y )/( 1 + dr * h_in ) ] );
% 
% %    helper = @(y) two_step( y, ps_trjs.t( i ), ps_trjs.x( i - 1 ), ps_trjs.x( i ), h_in, user );
% %    ps_trjs.x( i + 1 ) = fzero( helper, 0 );
%     ps_trjs.x( i + 1 ) = two_step( ps_trjs.t( i ), ps_trjs.x( i - 1 ), ps_trjs.x( i ), h_in, user );
%    
% end


function out = two_step( t, y0, y1, h, user )

% simulation parameters
p_in = user.params;
input = user.mdl.U;
dr = user.params( user.mdl.idp.damping_ratio );
nf = user.params( user.mdl.idp.natural_freq );
cr = user.params( user.mdl.idp.coefficient_restitution );
x_max = user.params( user.mdl.idp.x_max );


% building terms of final expression
num = h^2 * input( t, p_in, user ) + ( 2 - h^2 * nf^2 ) * y1 - ( ( 1 - cr ) - ( 1 + cr ) * dr * h ) * y0;
denom = 1 + dr * h;

out = -cr * y0 + min( [ num/denom; ( 1 + cr ) * x_max ] );



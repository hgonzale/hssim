user = struct;

% simulation parameters
user.step_size = 1e-2;
user.rx = user.step_size/100;
user.debug = 0;
user.max_jumps = inf;
user.euler_beta = 0.9;
user.euler_kmax = 1000;

% model setup
user.mdl = nav_mdl( 'nav.conf', user.debug );

% model parameters
user.time_arr = [ 0, 20 ];

% input parameters
user.A = [ -0.8 -0.1; -0.2 -0.8 ];
user.map = [ 2 4 6 6 6; 2 4 7 7 4; 2 4 -1 3 4; 2 4 6 6 6; 2 8 0 0 0 ]; % if -1 then obstacle, if 8 then goal

% set of initial conditions
user.setx0 = [ 3.0 3.5; 3.0 3.5 ];
% user.setx0 = [ 0 1; 2 3 ];
% user.setv0 = [ -0.3 0.3; -0.3 0 ];
user.setv0 = [ 0.5 0.5; -0.5 0.5 ];
user.x0 = zeros( 1, 4 );

% how fine a discretization to choose
nx = 5;
nv = 5;

out = 0;
% actually performing verification
for i = linspace( user.setx0( 1, 1 ), user.setx0( 1, 2 ), nx )
    for j = linspace( user.setx0( 2, 1 ), user.setx0( 2, 2 ), nx )
        for k = linspace( user.setv0( 1, 1 ), user.setv0( 1, 2 ), nv )
            for l = linspace( user.setv0( 2, 1 ), user.setv0( 2, 2 ), nv )
                user.x0( user.mdl.ids0.x:user.mdl.ids0.y ) = [ i; j ];
                user.x0( user.mdl.ids0.xdot:user.mdl.ids0.ydot ) = [ k; l ];
                user.params( user.mdl.idp.locx ) = min( floor( user.x0( user.mdl.ids0.x ) + 1 ), size( user.map, 2 ) );
                user.params( user.mdl.idp.locy ) = size( user.map, 1 ) - ...
                    min( floor( user.x0( user.mdl.ids0.y ) + 1 ), size( user.map, 1 ) ) + 1;
                user.params( user.mdl.idp.mode ) = user.map( user.params( user.mdl.idp.locy ), ...
                    user.params( user.mdl.idp.locx ) );

                trjs = fwd_RK2( user.mdl, user );
                if ( trjs( end ).p( user.mdl.idp.mode ) == -1 )
                    out = -1; %% verification not satisfied
                    break;
                elseif( trjs( end ).p( user.mdl.idp.mode ) == 8 )
                    out = 1; %% its possible to reach the goal set!
                elseif( trjs( end ).t( end ) >= user.time_arr( end ) )
                    fprintf( 'Need longer orbits\n' );
                else
%                     fprintf( 'Cyclic Orbit\n' );
                end
%                 fprintf( ' done with i = %d, j = %d, k = %d, l = %d \n', i, j, k, l ); 
            end
            if ( out == -1 )
                break;
            end
        end
        if ( out == -1 )
            break;
        end
       fprintf( ' done with i = %d, j = %d \n', i, j ); 
    end
    if ( out == -1 )
        break;
    end
end

if ( out == -1 )
    fprintf( 'Verification Completely FAILED!\n');
elseif( out == 1 )
    fprintf( 'Verification Succeeded!\n');
else
    fprintf( 'Verification Partially Failed! \n');
end


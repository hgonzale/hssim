% plot_navbox( user )
%
% Generate navigation box
%
% Inputs:
% user  - struct - holds all the simulation parameters
% nTraj - scalar - number of sample trajectories to plot
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function plot_navbox( user, nTraj )

if ( nargin < 2 )
    nTraj = 5;
end

% extract parameters
map = user.map;

% relevant params
xmax = size( map, 2 );
ymax = size( map, 1 );

% plotting parameters
qarrowsize = 1.5;
qwidth = 0.2;
lwidth = 1.1;
trajwidth = 1.5;
fontsize = 32;
colors = { 'm', 'b', 'c', [ 255 128 0 ]/255, [ 255 102 178 ]/255, ...
    [ 128 128 128 ]/255, [ 102 102 255 ]/255, [ 51 0 102 ]/255 };

% discretization vf
Ax = { ( -1 + qwidth ):qwidth:-qwidth, ...
    ( -1 + qwidth - cos( pi/4 ) * qwidth/2 ):qwidth: (- qwidth - cos( pi/4 ) * qwidth/2 ), ...
    ( -1 + qwidth ):qwidth:( -2 * qwidth ), ...
    ( -1 + qwidth - cos( pi/4 ) * qwidth/2 ):qwidth: (- qwidth - cos( pi/4 ) * qwidth/2 ), ...
    ( -1 + qwidth ):qwidth:-qwidth, ...
    ( -1 + qwidth + cos( pi/4 ) * qwidth/2 ):qwidth: (- qwidth + cos( pi/4 ) * qwidth/2 ), ...
    ( -1 + 2 * qwidth ):qwidth:( -qwidth ), ...
    ( -1 + qwidth + cos( pi/4 ) * qwidth/2 ):qwidth: (- qwidth + cos( pi/4 ) * qwidth/2 )};
Ay = { ( -1 + qwidth ):qwidth:( -2 * qwidth ), ...
    ( -1 + qwidth - sin( pi/4 ) * qwidth/2 ):qwidth: (- qwidth - sin( pi/4 ) * qwidth/2 ), ...
    ( -1 + qwidth ):qwidth:-qwidth, ...
    ( -1 + qwidth + sin( pi/4 ) * qwidth/2 ):qwidth: (- qwidth + sin( pi/4 ) * qwidth/2 ), ...
    ( -1 + 2 * qwidth ):qwidth:( -qwidth ), ...
    ( -1 + qwidth + sin( pi/4 ) * qwidth/2 ):qwidth: (- qwidth + sin( pi/4 ) * qwidth/2 ), ...
    ( -1 + qwidth ):qwidth:-qwidth, ...
    ( -1 + qwidth - sin( pi/4 ) * qwidth/2 ):qwidth: (- qwidth - sin( pi/4 ) * qwidth/2 )};

% screen size
scrsz = get( 0, 'ScreenSize' );

% plotting things
figure('Position', [ 20 scrsz(4)-700 900 770 ] ); hold on;

% plot straight lines to delineate different boxes
for i = 0:xmax
    plot( [ i i ], [ 0 ymax ], 'k', 'linewidth', lwidth );
end
for j = 0:ymax
    plot( [ 0 xmax ], [ j j ], 'k', 'linewidth', lwidth );
end

% set tick marks
set( gca, 'YTick', 0:ymax, 'FontSize', fontsize );
set( gca, 'XTick', 0:xmax, 'FontSize', fontsize );

% actually plotting
for i = 1:xmax
    for j = 1:ymax
        locx = i;
        locy = ymax - j + 1;
        if( map( locy, locx ) == -1 )
            rectangle( 'Position', [ i - 1, j - 1, 1 1 ], 'FaceColor' , 'r', 'EdgeColor', 'k', 'LineWidth', lwidth );
            text( ( 2 * i - 1 )/2, ( 2 * j - 1 )/2, 'Obstacle', 'FontSize', fontsize, 'HorizontalAlignment', 'Center' );
            continue;
        elseif( map( locy, locx ) == 8 )
            rectangle( 'Position', [ i - 1, j - 1, 1 1 ], 'FaceColor' , 'g', 'EdgeColor', 'k', 'LineWidth', lwidth );
            text( ( 2 * i - 1 )/2, ( 2 * j - 1 )/2, 'Goal', 'FontSize', fontsize, 'HorizontalAlignment', 'Center' );
            continue;
        end
        %         [ X, Y ] = meshgrid( linspace( ( i - 1 ), i, nx ), linspace( ( j - 1 ), j, ny ) );
        [ X , Y ] = meshgrid( Ax{ map( locy, locx ) + 1 } + i, Ay{ map( locy, locx ) + 1 } + j );
        DX = ones( size( X ) ) .* sin( map( locy, locx ) * pi/4 ) .* qwidth;
        DY = ones( size( Y ) ) .* cos( map( locy, locx ) * pi/4 ) .* qwidth;
        h = quiver( X, Y, DX, DY, 'AutoScale','off', 'linewidth', lwidth, ...
            'MaxHeadSize', 100, 'Color', colors{ map( locy, locx ) + 1 } );
        adjust_quiver_arrowhead_size( h, qarrowsize )
        text( ( 2 * i - 1 )/2, ( 2 * j - 1 )/2, int2str( map( locy, locx ) ), 'FontSize', fontsize, 'HorizontalAlignment', 'Center' );
    end
end

% plot sample trajectories w/o randomness
% nTraj = 2;
% for i = linspace( user.setx0( 1, 1 ), user.setx0( 1, 2 ), nTraj )
%     for j = linspace( user.setx0( 2, 1 ), user.setx0( 2, 2 ), nTraj )
%         user.x0( user.mdl.ids0.x:user.mdl.ids0.y ) = [ i; j ];
%
%
%
%         % this is an artificial choice
%         user.x0( user.mdl.ids0.xdot:user.mdl.ids0.ydot ) = [ ( user.setv0( 1, 1 ) + user.setv0( 1, 1 ) )/2; ...
%             ( user.setv0( 2, 1 ) + user.setv0( 2, 1 ) )/2 ];
%
%         user.params( user.mdl.idp.locx ) = min( floor( user.x0( user.mdl.ids0.x ) + 1 ), size( user.map, 2 ) );
%         user.params( user.mdl.idp.locy ) = size( user.map, 1 ) - ...
%             min( floor( user.x0( user.mdl.ids0.y ) + 1 ), size( user.map, 1 ) ) + 1;
%         user.params( user.mdl.idp.mode ) = user.map( user.params( user.mdl.idp.locy ), ...
%             user.params( user.mdl.idp.locx ) );
%
%         trjs = fwd_RK2( user.mdl, user );
%         plot( trjs( 1 ).x( 1, user.mdl.ids0.x ), trjs( 1 ).x( 1, user.mdl.ids0.y ), 'k.', 'MarkerSize', 30 );
%         for k = 1:length( trjs )
%             plot( trjs( k ).x( :, user.mdl.ids0.x ), trjs( k ).x( :, user.mdl.ids0.y ), 'k:', 'linewidth', trajwidth );
%         end
%         plot( trjs( end ).x( end, user.mdl.ids0.x ), trjs( end ).x( end, user.mdl.ids0.y ), 'kx', 'MarkerSize', 30 );
%     end
% end

for d = 1:nTraj
    i = user.setx0( 1, 1 ) + ( user.setx0( 1, 2 ) - user.setx0( 1, 1 ) ) * rand( 1, 1 );
    j = user.setx0( 2, 1 ) + ( user.setx0( 2, 2 ) - user.setx0( 2, 1 ) ) * rand( 1 ,1 );
    k = user.setv0( 1, 1 ) + ( user.setv0( 1, 2 ) - user.setv0( 1, 1 ) ) * rand( 1, 1 );
    l = user.setv0( 2, 1 ) + ( user.setv0( 2, 2 ) - user.setv0( 2, 1 ) ) * rand( 1, 1 );
    user.x0( user.mdl.ids0.x:user.mdl.ids0.y ) = [ i; j ];
%     user.x0( user.mdl.ids0.x:user.mdl.ids0.y ) = [ 3; 3 ];
    user.v0( user.mdl.ids0.xdot:user.mdl.ids0.ydot ) = [ k; l ];
    
    % this is an artificial choice
    user.x0( user.mdl.ids0.xdot:user.mdl.ids0.ydot ) = [ ( user.setv0( 1, 1 ) + user.setv0( 1, 1 ) )/2; ...
        ( user.setv0( 2, 1 ) + user.setv0( 2, 1 ) )/2 ];
    
    user.params( user.mdl.idp.locx ) = min( floor( user.x0( user.mdl.ids0.x ) + 1 ), size( user.map, 2 ) );
    user.params( user.mdl.idp.locy ) = size( user.map, 1 ) - ...
        min( floor( user.x0( user.mdl.ids0.y ) + 1 ), size( user.map, 1 ) ) + 1;
    user.params( user.mdl.idp.mode ) = user.map( user.params( user.mdl.idp.locy ), ...
        user.params( user.mdl.idp.locx ) );
    
    trjs = fwd_RK2( user.mdl, user );
    plot( trjs( 1 ).x( 1, user.mdl.ids0.x ), trjs( 1 ).x( 1, user.mdl.ids0.y ), 'k.', 'MarkerSize', 30 );
    for k = 1:length( trjs )
        plot( trjs( k ).x( :, user.mdl.ids0.x ), trjs( k ).x( :, user.mdl.ids0.y ), 'k-', 'linewidth', trajwidth );
    end
    plot( trjs( end ).x( end, user.mdl.ids0.x ), trjs( end ).x( end, user.mdl.ids0.y ), 'kx', 'MarkerSize', 30, 'LineWidth',2,...
                'MarkerEdgeColor','k' );
end


% plot axis;
box on;
% axis( [ -1e-2 xmax + 1e-3 -1e-2 ymax + 1e-2] );





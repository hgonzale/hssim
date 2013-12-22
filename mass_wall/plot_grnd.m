% plot_grnd( grnd_trjs, user )
%
% Plot the ground truth trajectories
%
% Inputs:
% grnd_trjs  - struct - ground truth trajectory data
% user       - struct - holds all the simulation parameters
% gt         - scalar - 0 if zeno and 1 if normal
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function plot_grnd( grnd_trjs, user, gt )

% ground truth trajectories
figure; hold on;
for i = 1:length( grnd_trjs )
    plot_interval = linspace( grnd_trjs( i ).ti, grnd_trjs( i ).tf, 100 );
    plot_value = [];
    for j = plot_interval
        plot_value = [ plot_value; grnd_trjs( i ).sol( j ) ];
    end
    plot( plot_interval, plot_value, 'k', 'LineWidth', 2 );
end

if ( gt )
    % plotting a single period for fmw_normal
    plotT = 3 * 2 * pi/user.input_freq;
    plot( linspace( 0, plotT, 1000 ), user.params( user.mdl.idp.x_max ) * ones( 1000, 1 ), 'k:', 'LineWidth', 2 );
    %
    % axis scaling for fmw normal
    axis( [ 0 plotT -16 15 ] );
else
    % plotting a single period for fmw_zeno
    plotT = 2 * pi/user.input_freq;
    plot( linspace( 0, plotT, 1000 ), user.params( user.mdl.idp.x_max ) * ones( 1000, 1 ), 'k:', 'LineWidth', 2 );
    
    % axis scaling for fmw zeno
    axis( [ 0 plotT -0.86 -0.799 ] );
end

% label axes
xlabel( 'Time [s]', 'FontSize', 20 );
set( gca, 'fontsize', 20 );
ylabel( 'x(t)', 'FontSize', 20 );
set( gca, 'fontsize', 20 );







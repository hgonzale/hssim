% plot_comp
%
% Plot the \rho, max_error, and computation times
%
% by Sam Burden, Humberto Gonzalez, Ram Vasudevan, 2013

function plot_comp

load fmw_zeno_rx1to4;
PS_error_zeno = PS_max_error( :, 1 );
PS_comptime_zeno = PS_comp_time( :, 1 );
rho_zeno = rho;
maxerror_zeno = max_error;
comptime_zeno = comp_time;
h_zeno = h;
rx_zeno = rx;

load fmw_normal;
PS_error_normal = PS_max_error( :, 1 );
PS_comptime_normal = PS_comp_time( :, 1 );
rho_normal = rho;
maxerror_normal = max_error;
comptime_normal = comp_time;
h_normal = h;
rx_normal = rx;

% plotting pattern
pp = { 'm-+', 'c-o', 'r-^', 'g-v', 'b-x', 'k-s' };

% whichidx to plot
whichidx = 5:3:length( h_zeno );
idx = 1;
lw = 1.5;
ms = 15;

% screen size
scrsz = get( 0, 'ScreenSize' );

%plot rho
figure('Position', [ 20 scrsz(4)-700 900 770 ] ); hold on;
plot( h_normal( whichidx ), rho_normal( whichidx, idx ), pp{ 3 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_normal( whichidx ), rho_normal( whichidx, end ), pp{ 4 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_zeno( whichidx ), rho_zeno( whichidx, idx ), pp{ 1 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_zeno( whichidx ), rho_zeno( whichidx, end ), pp{ 2 }, 'LineWidth', lw, 'MarkerSize', ms );
set( gca, 'YScale', 'log', 'XScale', 'log' );
% make legend
hleg1 = legend( 'Ex. 1, \epsilon = 2e-2 * h', 'Ex. 1, \epsilon = 2e-5 * h', ...
    'Ex. 2, \epsilon = 2e-2 * h', 'Ex. 2, \epsilon = 2e-5 * h','Location', 'SouthEast' );
set( hleg1, 'FontSize', 30  );

% label axes
xlabel( 'h', 'FontSize', 30 );
set( gca, 'fontsize', 30 );
ylabel( '$$\rho^{\epsilon}_{[0,t_{max}]}$$', 'FontSize', 40, 'Interpreter', 'latex' );
set( gca, 'fontsize', 30 );


%plot max
figure('Position', [ 20 scrsz(4)-700 900 770 ] ); hold on;
plot( h_normal( whichidx ), maxerror_normal( whichidx, idx ), pp{ 3 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_normal( whichidx ), maxerror_normal( whichidx, end ), pp{ 4 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_normal( whichidx ), PS_error_normal( whichidx ), pp{ 6 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_zeno( whichidx ), maxerror_zeno( whichidx, idx ), pp{ 1 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_zeno( whichidx ), maxerror_zeno( whichidx, end ), pp{ 2 }, 'LineWidth', lw, 'MarkerSize', ms );
plot( h_zeno( whichidx ), PS_error_zeno( whichidx ), pp{ 5 }, 'LineWidth', lw, 'MarkerSize', ms );
set( gca, 'YScale', 'log', 'XScale', 'log' );
% axis( [ 1e-4 1e1 1e-5 1e2 ] )
% make legend
hleg1 = legend( 'Ex. 1, \epsilon = 2e-2 * h', 'Ex. 1, \epsilon = 2e-5 * h', 'Ex. 1, PS-Method', ...
    'Ex. 2, \epsilon = 2e-2 * h', 'Ex. 2, \epsilon = 2e-5 * h', 'Ex. 2, PS-Method','Location', ...
    'SouthEast' );
set( hleg1, 'FontSize', 30  );

% label axes
xlabel( 'h', 'FontSize', 30 );
set( gca, 'fontsize', 30 );
ylabel( '$$\hat{\rho}$$', 'FontSize', 40, 'Interpreter', 'latex' );
set( gca, 'fontsize', 30 );


%plot comptime
figure('Position', [ 20 scrsz(4)-700 900 770 ] ); hold on;
plot( h_normal( whichidx ), comptime_normal( whichidx, 1 ), pp{ 3 }, 'LineWidth', lw, 'MarkerSize', ms  );
plot( h_normal( whichidx ), comptime_normal( whichidx, end ), pp{ 4 }, 'LineWidth', lw, 'MarkerSize', ms  );
plot( h_normal( whichidx ), PS_comptime_normal( whichidx ), pp{ 6 }, 'LineWidth', lw, 'MarkerSize', ms  );
plot( h_zeno( whichidx ), comptime_zeno( whichidx, 1 ), pp{ 1 }, 'LineWidth', lw, 'MarkerSize', ms  );
plot( h_zeno( whichidx ), comptime_zeno( whichidx, end ), pp{ 2 }, 'LineWidth', lw, 'MarkerSize', ms  );
plot( h_zeno( whichidx ), PS_comptime_zeno( whichidx ), pp{ 5 }, 'LineWidth', lw, 'MarkerSize', ms  );
set( gca, 'YScale', 'log', 'XScale', 'log' );
% make legend
hleg1 = legend( 'Ex. 1, \epsilon = 2e-2 * h', 'Ex. 1, \epsilon = 2e-5 * h', 'Ex. 1, PS-Method', ...
    'Ex. 2, \epsilon = 2e-2 * h', 'Ex. 2, \epsilon = 2e-5 * h', 'Ex. 2, PS-Method','Location', ...
    'NorthEast' );
set( hleg1, 'FontSize', 30  );

% label axes
xlabel( 'h', 'FontSize', 30 );
set( gca, 'fontsize', 30 );
ylabel( 'Computation Time [s]', 'FontSize', 30 );
set( gca, 'fontsize', 30 );





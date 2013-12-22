% matlabpool;

user = struct;

% simulation parameters
user.step_size = 1e-4;
user.rx = user.step_size/100;
user.debug = 0;
user.max_jumps = inf;
user.euler_beta = 0.7;
user.euler_kmax = 1000;

% model setup
user.mdl = forced_mass_wall_mdl( 'forced_mass_wall.conf', user.debug );

% model parameters
user.time_arr = [ 0 50 * 2 *pi/2.5 ];
user.params( user.mdl.idp.mode ) = 0;
user.params( user.mdl.idp.damping_ratio ) = 0.05;
user.params( user.mdl.idp.natural_freq ) = 2.5;
user.params( user.mdl.idp.coefficient_restitution ) = 0.9;
user.params( user.mdl.idp.x_max ) = 14;

% input parameters
user.input_amp = 20;
user.input_freq = 2.5;

% initial condition
user.x0( user.mdl.ids0.x ) = 11.36263;  % since we start in mode 0
user.x0( user.mdl.ids0.xdot ) = 31.40358; % since we start in mode 0

% construct the ground truth data (only needs to be done once)
grnd_trjs = ground_truth( user.mdl, user );
fprintf( 1, 'done with ground truth data construction\n' );

% step-sizes and relaxation values to simulate
nh = 50; % number of time steps to sample
nrx = 20; % number of relaxation steps to sample
h = 5 * logspace( -1, -4, nh );
rx = logspace( -1, -4, nrx );
max_error = zeros( nh, nrx );
rho = zeros( nh, nrx );
comp_time = zeros( nh, nrx );
PS_max_error = zeros( nh );
PS_comp_time = zeros( nh );
for i = 1:nh
    for j = 1:nrx
        tic;
        trjs = fwd_RK2( user.mdl, user, h( i ), 0.2 * h( i ) * rx( j ) );
        comp_time( i, j ) = toc;
        [ max_error( i, j ), rho( i, j ) ] = compute_error( grnd_trjs, trjs, user, 0.2 * h( i ) * rx( j ) );
    end
    
    
    tic;
    trjs = PS_method( user.mdl, user, h( i ) );
    PS_comp_time( i ) = toc;
    PS_max_error( i ) = PS_compute_error( grnd_trjs, trjs, user );
end

save( 'fmw_normal.mat' );


clc
clear all
%% Symbols
syms px py pz vx vy vz roll pitch yaw wx wy wz bgx bgy bgz bax bay baz ax ay az dt real 
syms ngx ngy ngz nax nay naz nbgx nbgy nbgz nbax nbay nbaz real
syms npx npy npz
%% State Variables
position = [px; py; pz];          % position
orientation = [roll; pitch; yaw]; % orientation
angular_vel = [wx; wy; wz];       % angular velocity
linear_vel = [vx; vy; vz];        % Linear velocity
linear_acc = [ax; ay; az];        % Linear Acceleration
bg = [bgx; bgy; bgz];             % bias gyroscope
ba = [bax; bay; baz];             % bias accerelometer
ng = [ngx; ngy; ngz];             % gaussian noise gyroscope
na = [nax; nay; naz];             % gaussian noise accelerometer
Qg = 0.01;                        % nbg ~ N(0,Qg)
Qba = 0.01;                       % nba ~ N(0,Qba)
nbg = [nbgx; nbgy; nbgz];         % bias gaussian gyroscope
nba = [nbax; nbay; nbaz];         % bias gaussian accerometer
g = [0; 0; -9.81];                % gravity

%% Rotation Matrix and angular velocity in body frame

Rx = [1, 0,         0;
      0, cos(roll), -sin(roll); 
      0, sin(roll), cos(roll)];

Ry = [cos(pitch), 0, sin(pitch); 
      0,          1, 0; 
     -sin(pitch), 0, cos(pitch)];

Rz = [ cos(yaw), -sin(yaw), 0;
       sin(yaw),  cos(yaw), 0;
       0,         0,        1];


%% Unit vectors defining the coordinate axes
x_hat = [1; 0; 0];
y_hat = [0; 1; 0];
z_hat = [0; 0; 1];

% Angular velocity in the body frame
G = simplify(horzcat( z_hat, Rz*y_hat, Rz*Ry*x_hat));
% Rotation ZYX
R = simplify(Rz * Ry * Rx);
G_inverse = simplify(inv(G));

%% State Model
x1 = position;
x2 = orientation;
x3 = linear_vel;
x4 = bg;
x5 = ba;
x = vertcat(x1,x2,x3,x4,x5);
u = vertcat(linear_acc, angular_vel);
n = vertcat(ng, na, nba, nbg);

%% Process Model
x_dot = [               x3;
        G_inverse*(angular_vel - x4 -ng);
        g + R*(linear_acc - x5 - na);
                       nbg;
                       nba              ];

%% Jacobians
A = simplify(jacobian(x_dot,x));
% B = simplify(jacobian(x_dot,u));
U = simplify(jacobian(x_dot,n));

% save('U_jac.mat', "U")
% save('A_jac.mat', "A")
% save('G_inv.mat', "G_inverse")
% save('Rzyx.mat', "R")
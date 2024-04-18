 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time
%% Declaring the state variables 
% px = uPrev(1,1);
% py = uPrev(2,1);
% pz = uPrev(3,1);
roll = uPrev(4,1);
pitch = uPrev(5,1);
yaw = uPrev(6,1);
% vx = uPrev(7,1);
% vy = uPrev(8,1);
% vz = uPrev(9,1);
wx = angVel(1,1);
wy = angVel(2,1);
wz = angVel(3,1);
ax = acc(1,1);
ay = acc(2,1);
az = acc(3,1);
g = [0 ; 0; -9.81];

x1 = uPrev(1:3,1);
x2 = uPrev(4:6,1);
x3 = uPrev(7:9,1);
x4 = uPrev(10:12,1);
x5 = uPrev(13:15,1);
am = acc;
wm = angVel;

%% Rotation Matrix for ZYX rotation
R = rotz(yaw)*roty(pitch)*rotx(roll);

%% Angular velocity in the body frame
G = [0, -sin(yaw), cos(pitch)*cos(yaw);
     0,  cos(yaw), cos(pitch)*sin(yaw);
     1,         0,         -sin(pitch)];
G_inverse = [(cos(yaw)*sin(pitch))/cos(pitch), (sin(pitch)*sin(yaw))/cos(pitch), 1;
                       -sin(yaw),                         cos(yaw), 0;
             cos(yaw)/cos(pitch),              sin(yaw)/cos(pitch), 0];

% Bias Assumed from tuning 
ngx = 0.01;
ngy = 0.01;
ngz = 0.01;
nax = 0.01;
nay = 0.01;
naz = 0.01;
bgx = 0.01;
bgy = 0.01;
bgz = 0.01;
bax = 0.01;
bay = 0.01;
baz = 0.01; 


nbg = [bgx; bgy; bgz];             % bias gyroscope
nba = [bax; bay; baz];             % bias accerelometer
ng = [ngx; ngy; ngz];              % gaussian noise gyroscope
na = [nax; nay; naz];              % gaussian noise accelerometer



x_dot = [          x3;
        G_inverse*(wm - x4 -ng);
        g + R*(am - x5 - na);
                   nbg;
                   nba      ];

%% At
At = [0, 0, 0,                                                                                                                                         0,                                                                                                                                      0,                                                                                                                                                                                    0, 1, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 1, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 1,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                  -(bgx*cos(yaw) + ngx*cos(yaw) + bgy*sin(yaw) - wx*cos(yaw) + ngy*sin(yaw) - wy*sin(yaw))/cos(pitch)^2,                                                                                (sin(pitch)*sin(yaw)*(bgx + ngx - wx))/cos(pitch) - (cos(yaw)*sin(pitch)*(bgy + ngy - wy))/cos(pitch), 0, 0, 0, -(cos(yaw)*sin(pitch))/cos(pitch), -(sin(pitch)*sin(yaw))/cos(pitch), -1,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                cos(yaw)*(bgx + ngx - wx) + sin(yaw)*(bgy + ngy - wy), 0, 0, 0,                          sin(yaw),                         -cos(yaw),  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                            - (cos(yaw)*sin(pitch)*(bgx + ngx - wx))/cos(pitch)^2 - (sin(pitch)*sin(yaw)*(bgy + ngy - wy))/cos(pitch)^2,                                                                                                      (sin(yaw)*(bgx + ngx - wx))/cos(pitch) - (cos(yaw)*(bgy + ngy - wy))/cos(pitch), 0, 0, 0,              -cos(yaw)/cos(pitch),              -sin(yaw)/cos(pitch),  0,                    0,                                                    0,                                                    0;
0, 0, 0, - (sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))*(bay - ay + nay) - (cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll))*(baz - az + naz), cos(yaw)*sin(pitch)*(bax - ax + nax) - cos(pitch)*cos(roll)*cos(yaw)*(baz - az + naz) - cos(pitch)*cos(yaw)*sin(roll)*(bay - ay + nay), (cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))*(bay - ay + nay) - (cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw))*(baz - az + naz) + cos(pitch)*sin(yaw)*(bax - ax + nax), 0, 0, 0,                                 0,                                 0,  0, -cos(pitch)*cos(yaw),   cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll), - sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch);
0, 0, 0,   (cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw))*(bay - ay + nay) + (cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))*(baz - az + naz), sin(pitch)*sin(yaw)*(bax - ax + nax) - cos(pitch)*cos(roll)*sin(yaw)*(baz - az + naz) - cos(pitch)*sin(roll)*sin(yaw)*(bay - ay + nay), (cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll))*(bay - ay + nay) - (sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))*(baz - az + naz) - cos(pitch)*cos(yaw)*(bax - ax + nax), 0, 0, 0,                                 0,                                 0,  0, -cos(pitch)*sin(yaw), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw),   cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw);
0, 0, 0,                                                                   cos(pitch)*sin(roll)*(baz - az + naz) - cos(pitch)*cos(roll)*(bay - ay + nay),                            cos(pitch)*(bax - ax + nax) + cos(roll)*sin(pitch)*(baz - az + naz) + sin(pitch)*sin(roll)*(bay - ay + nay),                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,           sin(pitch),                                -cos(pitch)*sin(roll),                                -cos(pitch)*cos(roll);
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0;
0, 0, 0,                                                                                                                                               0,                                                                                                                                      0,                                                                                                                                                                                    0, 0, 0, 0,                                 0,                                 0,  0,                    0,                                                    0,                                                    0];
Ft = eye(15) + At*dt;

%% Ut

Ut = [                          0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
-(cos(yaw)*sin(pitch))/cos(pitch), -(sin(pitch)*sin(yaw))/cos(pitch), -1,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
                         sin(yaw),                         -cos(yaw),  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
             -cos(yaw)/cos(pitch),              -sin(yaw)/cos(pitch),  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0, -cos(pitch)*cos(yaw),   cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll), - sin(roll)*sin(yaw) - cos(roll)*cos(yaw)*sin(pitch), 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0, -cos(pitch)*sin(yaw), - cos(roll)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw),   cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw), 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0,           sin(pitch),                                -cos(pitch)*sin(roll),                                -cos(pitch)*cos(roll), 0, 0, 0, 0, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 1, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 1, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 0, 0, 0, 1;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 1, 0, 0, 0, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 1, 0, 0, 0, 0;
                                0,                                 0,  0,                    0,                                                    0,                                                    0, 0, 0, 1, 0, 0, 0];

Qd = diag([ng(1,1) ng(2,1) ng(3,1) na(1,1) na(2,1) na(3,1) nbg(1,1) nbg(2,1) nbg(3,1) nbg(1,1) nbg(2,1) nbg(3,1)]) * dt;
uEst = uPrev + x_dot*dt;
covarEst = Ft*covarPrev*Ft' + Ut*Qd*Ut';

end

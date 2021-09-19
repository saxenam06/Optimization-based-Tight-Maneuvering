function [Fyf,Fyr] = Tire_model_param(v,vy,wz,delta_f)

%%%% Vehicle Parameters
mass = 1600;  % in kgs
lf=4.47/2; lr=4.47/2;  % front, rear axle to c.g. location in m
w=1.82; % Track width
H=0.51; % C.G. Height
Iz=2100; %Vehicle inertia about Z axis in kg-m^2
L=lf+lr; % Wheel Base
Cf=57000; % N/rad Front tire linear cornering stiffness
Cr=36700; %NN/rad Rear Tire Linear Cornering Stiffness


beta = atan(((lr./(lf+lr))*tan(delta_f)));
vp_x=v*cos(beta);
% 
% sf=delta_f-beta-(lf*wz/vp_x);
% sr=-beta+(lr*wz/vp_x);

% 
sf=delta_f-atan((vy+(lf*wz))/vp_x);
sr=-atan((vy-(lr*wz))/vp_x);


Fyf=2*Cf*sf;
Fyr=2*Cr*sr;

end

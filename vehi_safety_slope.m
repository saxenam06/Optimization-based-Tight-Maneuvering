function [Fzr_left, Fzr_right, Fzf_left, Fzf_right] = vehi_safety_slope(v,vy,wz,a,Fyf,Fyr)

mass = 1600;  % in kgs
lf=4.47/2; lr=4.47/2;  % front, rear axle to c.g. location in m
w=1.82; % Track width
H=0.51; % C.G. Height
Iz=2100; %Vehicle inertia about Z axis in kg-m^2
L=lf+lr; % Wheel Base
Cf=57000; % N/rad Front tire linear cornering stiffness
Cr=36700; %NN/rad Rear Tire Linear Cornering Stiffness

Mu=0.142*mass;
Muf=Mu/2;
Mur=Mu/2;
Ms=mass-Mu;

g=9.81;

Fzf0= (Ms*lr/L  +  Muf)*g;
Fzr0= (Ms*lf/L  +  Mur)*g;

uzx=806;%mass*H/L; %%

Fzf=Fzf0-uzx*(a-(vy*wz));
Fzr=Fzr0+uzx*(a-(vy*wz));

uzyf=675;%mass*H*lr/(w*L); %%
uzyr=1076;%mass*H*lf/(w*L); %%

dFzyf=uzyf*((Fyf+Fyr)/mass);
dFzyr=uzyr*((Fyf+Fyr)/mass);

Fzf_left=0.5*Fzf - dFzyf;
Fzf_right = 0.5*Fzf + dFzyf;
Fzr_left =0.5*Fzr - dFzyr;
Fzr_right  = 0.5*Fzr + dFzyr;

end

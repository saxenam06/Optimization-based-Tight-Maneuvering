clc
clear all
close all

addpath('C:\CASADI')
import casadi.*

%%%Planning of Prediction Horizon, time step
N=30;  % No of steps in future for which the states of the Ego Vehicle are predicted
dt=0.05;   % Time Step
Lanewidth=3.5;

%%%%%Define Linear Obstacles: road boundary coordinates
ref_obs{1}=[linspace(0,120,240)' linspace(-3.5,-3.5,240)'];
ref_obs{2}=[linspace(0,120,240)' linspace(3.5,3.5,240)'];

%%%Define Reference Lane Waypoints for Vehicle to follow
ref_traj= ((ref_obs{1}+ref_obs{2})/2)';
ref_traj(2,:)=ref_traj(2,:)-1.75;
ref_traj=[ref_traj repmat(ref_traj(:,end),1,2000)];
traj_states=2;  % Only two states (X,Y) are given for the planned trajectory
%%% N+1 points of the Lane are given as reference to be achieved by the respective N+1 steps in the predicted
%%% trajectory of the ego vehicle 
n_traj_states=(N+1)*traj_states; % Total no of parameters to be given for the reference trajectory

%%%% Define no of obstacles and total obstacle states
obs_states=2;
no_of_obstacles=3;
%%%%%%%%Initial Position of Moving obstacles
CG(1:no_of_obstacles,1:2)=[25.5 41.0 70.5; -1.75 -1.75+3.5 -1.75]';


%%%%%%%STep: Start Defining Model Parameter and State Space model
%%%% Vehicle Parameters
mass = 1600;  % in kgs
lf=4.47/2; lr=4.47/2;  % front, rear axle to c.g. location in m
w=1.82; % Track width
H=0.51; % C.G. Height
Iz=2100; %Vehicle inertia about Z axis in kg-m^2
L=lf+lr; % Wheel Base
Cf=57000; % N/rad Front tire linear cornering stiffness
Cr=36700; %NN/rad Rear Tire Linear Cornering Stiffness
%%%%%Ego Vehicle Dynamic model or System model used for making Predictions
Pos_X=SX.sym('s1');
Pos_Y=SX.sym('s2');
theta=SX.sym('theta');
v=SX.sym('v');
vy=SX.sym('vy');
wz=SX.sym('wz');
a=SX.sym('acce');
delta_f=SX.sym('del_f');
states=[Pos_X;Pos_Y;theta;v;vy;wz;a;delta_f];
nvehs=length(states);
mvobs_states=SX.sym('mvobs',no_of_obstacles*obs_states,1);
states=[states;mvobs_states];
n_states=length(states);
n_mvobs_velo=no_of_obstacles;
n_mvobs_theta=no_of_obstacles;
n_ref_states=0;  % not used here reference goal location is set explicitly in the waypoints itself

%%%% Define P which will contain all the parameters that will be used by
%%% MPC: Initial position of ego vehicle, Reference trajectory states of Lane to be followed,
%%%% Target States of the Goal (not used here),
%%%% Velocity of all the obstacles and Angular position of all the obstacles
P=SX.sym('P',n_states+n_traj_states+n_ref_states+n_mvobs_velo+n_mvobs_theta);

%%% define control iputs
jerk=SX.sym('Jk'); 
str_rate=SX.sym('str_rate');
T=SX.sym('Tp');
controls=[jerk;str_rate;T];
n_controls=length(controls);

%%% Define State space model of Ego vehicle
beta = atan(((lr./(lf+lr))*tan(delta_f)));
xdot= v*cos(theta+beta);
ydot= v*sin(theta+beta);
thetadot=wz+(v*sin(beta)/lr);
vdot=a;
[Fyf, Fyr]=Tire_model_param(v,vy,wz,delta_f);
vydot= (Fyf+Fyr)/mass -(v*cos(beta)*wz);
wzdot= (Fyf*lf -Fyr*lr)/Iz;
adot=jerk;
delta_dot=str_rate;

%%% Define State space model of Moving Obstacles
obs_vel=P(n_states+n_traj_states+n_ref_states+1:n_states+n_traj_states+n_ref_states+n_mvobs_velo);
obs_theta=P(n_states+n_traj_states+n_ref_states+n_mvobs_velo+1:n_states+n_traj_states+n_ref_states+n_mvobs_velo+n_mvobs_theta);
mvobsd=[];
for i=1:no_of_obstacles
    mvobsdx=obs_vel(i)*cos(obs_theta(i));
    mvobsdy=obs_vel(i)*sin(obs_theta(i));
    mvobsd=[mvobsd mvobsdx mvobsdy];
end

%%%Complete system model with states of the Eqo vehicle and the States of
%%%the moving obstacles
sys_model=[xdot;ydot;thetadot;vdot;vydot;wzdot;adot;delta_dot;mvobsd'];
model_fn=Function('f',{states,controls,obs_vel,obs_theta},{sys_model});

%%%Bounding box around the vehicle as a function of the Vehicle Position and yaw angle
%%%%Modeled as a polytopic set Ax<=b
R = [cos(theta), -sin(theta); 
     sin(theta), cos(theta)];
b = [L/2 w/2 L/2 w/2]';
A = [R'; -R']; 
b = b + A*[Pos_X;Pos_Y]; 
box=[A, b];
Bound_box=Function('b',{[Pos_X;Pos_Y],theta},{box});
dmin=0.3; % minimum distance to be maintained between vehciles


%%%%% Step: Start formulating MPC problem using only symbolic variables%%%
%%%%%Define control variables used for complete prediction horizon
U=SX.sym('U',n_controls,N);
%%Define state variables used for complete prediction horizon
pred_state_matrix=SX.sym('X_pred',n_states,(N+1));

%%%%  Target States of the Goal (not used here) Goal location defined
%%%%  explicitly in the Reference trajectory see later
n_ref_states=0;  % x-Tar, y_tar

%%%% No. of vehicle Safety Wheel Lift-off constraints
n_veh_safty_states=4;
%%%% No. of lagrange multipliers for inequality Constraints
nlambda=4;
U_lmvo=SX.sym('Ulmvo',nlambda*no_of_obstacles,N+1); % Lagrange multiplier for Vehicle Polytopic set
U_lmov=SX.sym('Ulmov',nlambda*no_of_obstacles,N+1);  % Lagrange multiplier for Obstacle Polytopic set
%%%% No. of lagrange variables for the dual problem
ns=2;
U_svo=SX.sym('Usvo',ns*no_of_obstacles,N+1);

%%%% Define Objective function and Constraint Vector
obj=0; %Initialize objective function 
g=[]; %Constraint Vector

%%%%Integral gains
Q0=[0.1 0; 0 0.1]; % Weight for reference trajectory tracking
Q_delta=0.005;   %Weight for the Steering Angle
R=zeros(3,3); R(1,1)=0.005; R(2,2)=0.0005;  % Weight for control actions Jerk Rate and steering rate
R(3,3)=0.05; % Weight for time step
w_Fz=0.005; % Weight for Wheel lift off constraint

%%%% Vehicle Safety Constraint terms
Fz_thr=1000; % Threshhold vertical load for Vehicle Safety Constraint in N
Fz_off=100;
a_term=Fz_thr+3*Fz_off; % Smoothing
b_term=Fz_off;

Lidar_dist_limit=15; % Required minimum view of 15m at all the speeds
n_pred_len=1; %% Only one preview constraint of Lidar Limit of 15m

%%%Define Objective functions and Constraints for MPC
states_curr=pred_state_matrix(:,1);  
%%%% Constraint of initial state with initial position of ego vehicle
g=[g;states_curr-P(1:n_states)];

%%%% Define objective function
for k=1:N
    states_curr=pred_state_matrix(:,k);
    control_curr=U(:,k);
    reftraj_curr=[P(n_states+(2*k-1)); P(n_states+(2*k))];
    
    %%%Minimise distance between each of the N+1 points of the Lane given as reference to be achieved by the respective N+1 steps in the predicted
    %%% trajectory of the ego vehicle 
    states_term=(states_curr(1:2) - reftraj_curr(1:2));
    obj=obj + (states_term'*Q0*states_term)+(N*control_curr(3));
    
    %%%%Calculate Tire lateral forces as a function of Longitudinal &
    %%%%Lateral speed, Yaw rate & Steering Angle
    [Fyf,Fyr]=Tire_model_param(states_curr(4),states_curr(5),states_curr(6),states_curr(8));
    
    %%%%Calculate Tire Vertical forces as a function of Longitudinal &
    %%%%Lateral speed, Yaw rate, Longitudinal Acceleration & Lateral forces   
    
    [Fzr_left,Fzr_right,Fzf_left,Fzf_right]= vehi_safety_slope(states_curr(4),states_curr(5),states_curr(6),states_curr(7),Fyf,Fyr);
    
    %% Calculate cost term on Wheel lift off If greater than 1300N, -ve cost, If less than 1300 positive cost
    Fzr_L_term=tanh(-(Fzr_left-a_term)/b_term);
    Fzr_R_term=tanh(-(Fzr_right-a_term)/b_term);
    Fzf_L_term=tanh(-(Fzf_left-a_term)/b_term);
    Fzf_R_term=tanh(-(Fzf_right-a_term)/b_term);
    ARC_term2=Fzr_L_term+Fzr_R_term+Fzf_L_term+Fzf_R_term;
    %%% Other cost terms on control, Wheel lift off Cost term and cost on
    %%% large steering angle 
    obj=obj+control_curr'*R*control_curr + w_Fz*ARC_term2 + Q_delta*states_curr(8)^2;
    
    states_next=pred_state_matrix(:,k+1);
    %%%%Next step state prediction and constraint
    model_fn_value=model_fn(states_curr, control_curr, obs_vel, obs_theta);
    states_next_euler= states_curr + (control_curr(3)*model_fn_value);
    %%%%Next step state prediction and constraint
    g=[g;(states_next-states_next_euler)];
    %%%%Constraint for Wheel lift off to be greater than 1000N
    g=[g;-(Fzr_left-Fz_thr)];
    g=[g;-(Fzr_right-Fz_thr)];
    g=[g;-(Fzf_left-Fz_thr)];
    g=[g;-(Fzf_right-Fz_thr)];
    %%%%Constraint for minimum Lidar preview to be at 15m during longitudinal speed variation by adjusting the time step 
    pred_dist=states_curr(4)*control_curr(3)*N;
    g=[g;-(pred_dist-Lidar_dist_limit)];
     %%%%Collision avoidacne Constraints assure minimum distance between the polytopic sets of Vehicle and Obstacle to be
     %%%greater than dmin for their entire prediction horizon, Defined as dual problem 
    for O=1:no_of_obstacles
        Bv=Bound_box([states_curr(1);states_curr(2)],states_curr(3));
        Av{k}=Bv(:,1:2);
        bv{k}=Bv(:,3);
        Bo=Bound_box([states_curr(nvehs+((2*O)-1));states_curr(nvehs+(2*O))],obs_theta(O));
        Ao{k}=Bo(:,1:2);
        bo{k}=Bo(:,3);
        g=[g;bv{k}'*U_lmvo(nlambda*O-3:nlambda*O,k)+bo{k}'*U_lmov(nlambda*O-3:nlambda*O,k)];
        g=[g;Av{k}'*U_lmvo(nlambda*O-3:nlambda*O,k)+ U_svo(ns*O-1:ns*O,k)];
        g=[g;Ao{k}'*U_lmov(nlambda*O-3:nlambda*O,k)- U_svo(ns*O-1:ns*O,k)];
        g=[g; ((U_svo(ns*O-1,k)^2 + U_svo(ns*O,k)^2)) ];
    end
end
  
reftraj_next=[P(n_states+(2*k+1));P(n_states+(2*k+2))];
%%%Minimise distance between last points of the Lane given as reference to be achieved by the respective last step in the predicted
%%% trajectory of the ego vehicle 
s_f=[states_next(1:2)-reftraj_next(1:2)];
obj=obj+s_f'*Q0*s_f;
%%
%%%%Collision avoidacne Constraints assure minimum distance between the polytopic sets of Vehicle and Obstacle to be
     %%%greater than dmin for the last step in the predicion horizon, Defined as dual problem 
for O=1:no_of_obstacles
        Bv=Bound_box([states_next(1);states_next(2)],states_next(3));
        Av{N+1}=Bv(:,1:2);
        bv{N+1}=Bv(:,3);
        Bo=Bound_box([states_next(nvehs+((2*O)-1));states_next(nvehs+(2*O))],obs_theta(O));
        Ao{N+1}=Bo(:,1:2);
        bo{N+1}=Bo(:,3); 
        g=[g;bv{N+1}'*U_lmvo(nlambda*O-3:nlambda*O,N+1)+bo{N+1}'*U_lmov(nlambda*O-3:nlambda*O,N+1)];
        g=[g;Av{N+1}'*U_lmvo(nlambda*O-3:nlambda*O,N+1)+ U_svo(ns*O-1:ns*O,N+1)];
        g=[g;Ao{N+1}'*U_lmov(nlambda*O-3:nlambda*O,N+1)- U_svo(ns*O-1:ns*O,N+1)];
        g=[g; ((U_svo(ns*O-1,N+1)^2 + U_svo(ns*O,N+1)^2)) ];
end
%%
%%%%Function for longitudinal acceleration constraint as a function of
%%%%longitudinal speed similar to the thesis
c1 = -1.28e-4; c2 =8.59e-3; c3 =-0.2257; c4 =3.0828; c5=-1.38e-4; c6 =6.85e-3; c7 =-0.1204; c8 =-3.5589;
ax_max = c1*v^3+ c2*v^2+ c3*v+c4;
ax_min = c5*v^3+ c6*v^2+ c7*v+c8;
constraint_model = [ ax_max ; ax_min ];
con_fn = Function ('f' ,{ states },{ constraint_model });

%%%%% Optimiation variables as all Ego Vehicle states and the Control
%%%%% variables, lagrange multipliers of the dual problem
OPT_variables = [ reshape( pred_state_matrix,(n_states)*(N+1),1); reshape(U,n_controls*N,1);...
    reshape(U_lmvo,nlambda*no_of_obstacles*(N+1),1);reshape(U_lmov,nlambda*no_of_obstacles*(N+1),1);reshape(U_svo,ns*no_of_obstacles*(N+1),1)];
probNLP = struct ('f',obj , 'x', OPT_variables , 'g', g, 'p', P) ;


%%%Define Solver settings and assign it to a structure
opts = struct ;
opts.ipopt.max_iter = 200; 
opts.ipopt.print_level = 0; %0,3
opts.print_time = 0; %0, 1
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6; %
S = nlpsol ('solver','ipopt',probNLP, opts ) ;

%%%%%%%%%%MPC Problem formulation over%%%%%%%%%%%

%%%%%%%%% Input values into the MPC problem %%%%%%%%%%%%%
%%%%%%% Define Bounds on constraints
args = struct ;
args.lbg (1: n_states ) = 0;
args.ubg (1: n_states ) = 0;


%%%Define array which can be used to later refer the sequence of State
%%%Equality constraints, Vehicle Safety constraint & Lidar distance limit
%%% constraints and  lagrange multiplers of the collision avidance constraints 
factor = n_states ;
pred_st_con_len = zeros(n_states,N) ;
pred_safty_con_len = zeros(n_veh_safty_states,N) ;
pred_dist_con_len = zeros(n_pred_len,N) ;
ndistmin=1;
ndistspatialv=2;
ndistspatialo=2;
ndistdual=1;

pred_min=zeros(ndistmin,N+1);
pred_distspatialv=zeros(ndistspatialv,N+1);
pred_distspatialo=zeros(ndistspatialo,N+1);
pred_dual=zeros(ndistdual,N+1);

for k = 1:N
pred_st_con_len(:,k)=(factor+1):( factor+n_states ) ;
pred_safty_con_len(:,k) = ((factor+n_states)+1):(factor+(n_states+n_veh_safty_states));
pred_dist_con_len(:,k)=(factor+(n_states+n_veh_safty_states)+1):(factor+(n_states+n_veh_safty_states+n_pred_len));
for O=1:no_of_obstacles
    pred_distmin(O,k)=(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin));
    pred_distspatialv(ndistspatialv*O-1:ndistspatialv*O,k)=(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin)+1):(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin+ndistspatialv));
    pred_distspatialo(ndistspatialo*O-1:ndistspatialo*O,k)=(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin+ndistspatialv)+1):(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin+ndistspatialv+ndistspatialo));
    pred_dual(O,k)=(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin+ndistspatialv+ndistspatialo)+1):(factor+(n_states+n_veh_safty_states+n_pred_len+ndistmin+ndistspatialv+ndistspatialo+ndistdual));
    factor = factor+(ndistmin+ndistspatialv+ndistspatialo+ndistdual);
end
factor=factor+(n_states+n_veh_safty_states+n_pred_len );
end

%% Redo for N+1 step 
for O=1:no_of_obstacles
    pred_distmin(O,k+1)=(factor+ndistmin);
    pred_distspatialv(ndistspatialv*O-1:ndistspatialv*O,k+1)=(factor+ndistmin+1):(factor+ndistmin+ndistspatialv);
    pred_distspatialo(ndistspatialo*O-1:ndistspatialo*O,k+1)=(factor+ndistmin+ndistspatialv+1):(factor+ndistmin+ndistspatialv+ndistspatialo);
    pred_dual(O,k+1)=factor+ndistmin+ndistspatialv+ndistspatialo+ndistdual;
    factor = factor+(ndistmin+ndistspatialv+ndistspatialo+ndistdual);
end

%%%%Defining bounds on State equality cinstraints
args.lbg(reshape(pred_st_con_len ,1,(n_states *N))) = 0;
args.ubg(reshape(pred_st_con_len ,1,(n_states *N))) = 0;
%%%%Defining bounds on Wheel lift off constraint
args.lbg(reshape(pred_safty_con_len ,1,(n_veh_safty_states*N))) = -inf ;
args.ubg(reshape(pred_safty_con_len ,1,(n_veh_safty_states*N))) = -5;
%%%%Defining bounds on Lidar Distance limit
args.lbg(reshape(pred_dist_con_len ,1,( n_pred_len*N))) = -1;
args.ubg(reshape(pred_dist_con_len ,1,( n_pred_len*N))) = 1;
%% Define upper and lower bounds for collision avoidance constraints reformulated as a Lagrange dual problem
%%% to assure minimum distance between the polytopic sets of Vehicle and Obstacle 
args.lbg(reshape(pred_distmin ,1,( ndistmin*O*(N+1)))) = -inf;
args.ubg(reshape(pred_distmin ,1,( ndistmin*O*(N+1)))) = -dmin;
args.lbg(reshape(pred_distspatialv ,1,( ndistspatialv*O*(N+1)))) = 0;
args.ubg(reshape(pred_distspatialv ,1,( ndistspatialv*O*(N+1)))) = 0;
args.lbg(reshape(pred_distspatialo ,1,( ndistspatialo*O*(N+1)))) = 0;
args.ubg(reshape(pred_distspatialo ,1,( ndistspatialo*O*(N+1)))) = 0;
args.lbg(reshape(pred_dual ,1,( ndistdual*O*(N+1)))) = 0;
args.ubg(reshape(pred_dual ,1,( ndistdual*O*(N+1)))) = 1;

%%Define bounds for the Optimization variables
args.lbw(1: n_states : n_states *(N +1) ,1) = -1;
args.ubw(1: n_states : n_states *(N +1) ,1) = 120;
args.lbw(2: n_states : n_states *(N +1) ,1) = -3.5+w/2;
args.ubw(2: n_states : n_states *(N +1) ,1) = 3.5-w/2;
args.lbw(3: n_states : n_states *(N +1) ,1) = -2*pi;
args.ubw(3: n_states : n_states *(N +1) ,1) = 2*pi;
args.lbw(4: n_states : n_states *(N +1) ,1) = 0.001;
args.ubw(4: n_states : n_states *(N +1) ,1) = 15.5;
args.lbw(5: n_states : n_states *(N +1) ,1) = -inf;
args.ubw(5: n_states : n_states *(N +1) ,1) = inf ;
args.lbw(6: n_states : n_states *(N +1) ,1) = -1;
args.ubw(6: n_states : n_states *(N +1) ,1) = 1 ;
args.lbw(8: n_states : n_states *(N +1) ,1) = -pi/4;
args.ubw(8: n_states : n_states *(N +1) ,1) = pi/4;
for O=1:no_of_obstacles
    args.lbw(nvehs+((2*O)-1): n_states : n_states*(N +1) ,1) = -1010;
    args.ubw(nvehs+((2*O)-1): n_states : n_states*(N +1) ,1) = 1010;
    args.lbw(nvehs+((2*O)): n_states : n_states*(N +1) ,1) = -1010;
    args.ubw(nvehs+((2*O)): n_states : n_states*(N +1) ,1) = 1010;
end
    
J_max = 5; J_min = - J_max ;
str_rate_max =50*pi/180; str_rate_min =-str_rate_max ;
%%%%% Bounds on Control inputs
args.lbw(n_states*(N +1)+1: n_controls : n_states*(N +1)+n_controls*N ,1) = J_min ;
args.ubw(n_states*(N +1)+1: n_controls : n_states*(N +1)+n_controls*N ,1) = J_max ;
args.lbw(n_states*(N +1)+2: n_controls : n_states*(N +1)+n_controls*N ,1) = str_rate_min ;
args.ubw(n_states*(N +1)+2: n_controls : n_states*(N +1)+n_controls*N ,1) = str_rate_max ;
args.lbw(n_states*(N +1)+3: n_controls : n_states*(N +1)+n_controls*N ,1) = 0.05;
args.ubw(n_states*(N +1)+3: n_controls : n_states*(N +1)+n_controls*N ,1) = 1;
%% Bounds on the lagrange mulitpliers of the dual problem
args.lbw(n_states*(N +1)+n_controls*N+1 : n_states*(N +1)+n_controls*N+nlambda*no_of_obstacles*(N+1),1) = 0;
args.ubw(n_states*(N +1)+n_controls*N+1 : n_states*(N +1)+n_controls*N+nlambda*no_of_obstacles*(N+1),1) = inf;
args.lbw(n_states*(N +1)+n_controls*N+nlambda*no_of_obstacles*(N+1)+1 : n_states*(N +1)+n_controls*N+(nlambda+nlambda)*no_of_obstacles*(N+1),1) = 0;
args.ubw(n_states*(N +1)+n_controls*N+nlambda*no_of_obstacles*(N+1)+1 : n_states*(N +1)+n_controls*N+(nlambda+nlambda)*no_of_obstacles*(N+1),1) = inf;
args.lbw(n_states*(N +1)+n_controls*N+(nlambda+nlambda)*no_of_obstacles*(N+1)+1 : n_states*(N +1)+n_controls*N+(nlambda+nlambda+ns)*no_of_obstacles*(N+1),1) = -1;
args.ubw(n_states*(N +1)+n_controls*N+(nlambda+nlambda)*no_of_obstacles*(N+1)+1 : n_states*(N +1)+n_controls*N+(nlambda+nlambda+ns)*no_of_obstacles*(N+1),1) = 1;

%%%% Initial States
x0 = [0.0; -1.75; 0; 5.5; 0.0; 0.0; 0.0; 0.0; reshape(CG', no_of_obstacles*2,1)];
%Initialise Dynamic constraint on longitudinal acceleration as a function
%of states
con_fn_value=con_fn(x0);
args.lbw(7: n_states : n_states *(N +1) ,1) = full(con_fn_value(2));
args.ubw(7: n_states : n_states *(N +1) ,1) = full(con_fn_value(1));
%%%%%Define Goal Location
x_ref = [120; -1.75];

%%%%% Define Variable to store Control actions determined by Optimizer at
%%%%% each time step
u0_int=[0.0;0.0;0.05];
ux(:,1)=u0_int;

%%%%% Define Variable to store evolution of Vehicle after the determined
%%%%% control actions at each time step
xx(:,1)=x0;
%%%Define Variable required for
%%%initialisation of Optimization Variables at every time step
u0=repmat(u0_int,1,N)';
pred_state_mat_init=repmat(x0,1,N+1)';
%%%Define Variable required for
%%%initialisation of lagrange multipliers at every time step
ulmvo_init=zeros(nlambda*no_of_obstacles,1);
ulmov_init=zeros(nlambda*no_of_obstacles,1);
usvo_init=zeros(ns*no_of_obstacles,1);
ulmvo0=repmat(ulmvo_init,1,N+1)';
ulmov0=repmat(ulmov_init,1,N+1)';
usvo0=repmat(usvo_init,1,N+1)';

ulmvoS=reshape(ulmvo0,nlambda*no_of_obstacles*(N+1),1);
ulmovS=reshape(ulmov0,nlambda*no_of_obstacles*(N+1),1);
usvoS=reshape(usvo0,ns*no_of_obstacles*(N+1),1);
ulm0=[ulmvoS;ulmovS;usvoS];

%%%%%Define variable to store predicted trajectory along the entire
%%%%%prediction horizon at each time step
xx1=[];

%%%%%Define variable to store predicted control trajectory along the entire
%%%%%prediction horizon at each time step
u_cl=[];

veh_sfty_N(:,1)=[3;3;3;3]; % arbritary values above 1KN
veh_sfty_N1= [];

%%%%Start MPC
%%%For Counting iterations
mpciter=0;
%%%%%%Define Initial time
t0=0;
t(1)=t0;

main_loop=tic;
%%%Define variable for storing points of the Lane  given as reference to be 
%%%achieved by the respective N+1 steps in the predicted
%%% trajectory of the ego vehicle 

ref_traj_cls=zeros(2,N+1);
ref_traj_cls_X=zeros(2,N+1);  
ref_traj_cls_Y=zeros(2,N+1);
%%%%%Define moving obstacle states values as per the giveen scenario
mv_obs_velo_ms=1.0*[1.5;3.5;1.0];
mv_obs_theta_rad=0.*ones(no_of_obstacles,1);

while norm((x0(1:2)-x_ref(1:2)),2)>2.e-1 %&& (mpciter < 10*sim_time/dt)
%%%Define points of the Lane that will be given as reference to be 
%%%achieved by the respective N+1 states in the predicted
%%% trajectory of the ego vehicle 


for idx=1:size(pred_state_mat_init,1)
    val_x=pred_state_mat_init(idx,1);
    [matchvalx,matchidx]=min(abs(ref_traj(1,:) - val_x));
    ref_traj_cls_X(1,idx)=ref_traj(1,matchidx);
    ref_traj_cls_X(2,idx)=ref_traj(2,matchidx);
    
    val_y=pred_state_mat_init(idx,2);
    [matchvaly,matchidy]=min(abs(ref_traj(2,:) - val_y));
    ref_traj_cls_Y(2,idx)=ref_traj(2,matchidy);
    ref_traj_cls_Y(1,idx)=ref_traj(1,matchidy);
    
    dist_X=sqrt((ref_traj_cls_X(1,idx)-x_ref(1))^2 + (ref_traj_cls_X(2,idx)-x_ref(2))^2);
    dist_Y=sqrt((ref_traj_cls_Y(1,idx)-x_ref(1))^2 + (ref_traj_cls_Y(2,idx)-x_ref(2))^2);
    
    if dist_X<dist_Y
        ref_traj_cls(:,idx)=ref_traj_cls_X(:,idx);
    else
        ref_traj_cls(:,idx)=ref_traj_cls_Y(:,idx);
    end   
end


%%% Once the distance of the ego vehicle from the goal location becomes
%%% lesser than the minimum lidar distance limit
%%% we will change values of the constraint bounds as follows
%%%% Switch to constant time step from varying time step i.e. no constraint on the
%%%% preview distance as it may now reduce with the reduction of the speed of the ego
%%%% vehicle, Minimum longitudinal speed of ego vehicle to 0.0001,
%%%%  No bounds on longitudinal acceleration and Jerk rate
%%% Reference tajectory is nothing but now a constant goal location 
if  norm((x0(1:2)-x_ref(1:2)),2)<=Lidar_dist_limit
args.ubw(n_states*(N+1)+3: n_controls : n_states*(N +1)+n_controls*N ,1) = 0.05;
args.lbw(4: n_states : n_states *(N +1) ,1) = 0.0001;
args.lbw(7: n_states : n_states *(N +1) ,1) = -inf;
args.ubw(7: n_states : n_states *(N +1) ,1) = inf;
args.lbg(reshape(pred_dist_con_len ,1,( n_pred_len *N))) = -inf;
args.ubg(reshape(pred_dist_con_len ,1,( n_pred_len *N))) = inf;
args.lbw(n_states *(N+1)+1: n_controls : n_states*(N +1)+n_controls*N ,1) = -inf ;
args.ubw(n_states *(N+1)+1: n_controls : n_states*(N +1)+n_controls*N ,1) = inf ;
ref_traj_cls=repmat(ref_traj(:,end),1,N+1);
end
%%%Give the above defined parameter values to the MPC 
%%%% Initial position of ego vehicle, Reference trajectory states of Lane to be followed,
%%%% Velocity of all the obstacles, Angular position of all the obstacles
args.p = [x0; reshape(ref_traj_cls,n_traj_states,1); mv_obs_velo_ms;mv_obs_theta_rad];
% Initialize optimization varibales in the first iteration 
%%%In the next iteration initialize with the
% trajectory calculated from MPC solution
args.w_init = [reshape(pred_state_mat_init', n_states*(N+1) ,1); reshape(u0',n_controls*N ,1);ulm0]; 
sol = S('x0',args.w_init, 'lbx', args.lbw , 'ubx', args.ubw,...
'lbg',args.lbg, 'ubg', args.ubg, 'p', args.p);
% Get solution of entire trajectory of controls and store for each MPC iteration
u = reshape(full(sol.x(n_states*(N+1)+1:n_states*(N+1)+n_controls*N))',n_controls ,N)'; 
u_cl (:,1:n_controls , mpciter +1) = u;
% Get solution of entire trajectory of Ego vehicle states and store for each MPC iteration
xx1 (:,1: n_states , mpciter +1) = reshape(full(sol.x(1:n_states*(N+1)))',n_states ,N+1)';
%%%Store constraints resulting values
% veh_sfty = (-1*( reshape ( full ( sol.g( reshape (pred_safty_con_len ,1,(n_veh_safty_states*N))))',n_veh_safty_states,N)') + Fz_thr )/1000;
% veh_lid = (1*( reshape ( full ( sol.g( reshape (pred_dist_con_len ,1,(n_pred_len*N))))',n_pred_len,N)' ));
% veh_sfty_N1 (:,1: n_veh_safty_states , mpciter +1) = veh_sfty ;
% veh_dist_lid(:,mpciter+1) = veh_lid;
% obs_dist_constr(:,1:no_of_obstacles)=reshape(1*full(sol.g(n_states*(N+1)+(n_veh_safty_states*N)+(n_pred_len*N)+1: ...
%     n_states*(N+1)+(n_veh_safty_states*N)+(n_pred_len*N)+no_of_obstacles*(N+1)))',no_of_obstacles,N+1)';
% obs_dist_lin(:,1)=reshape(1*full(sol.g(n_states*(N+1)+(n_veh_safty_states*N)+(n_pred_len*N)+(no_of_obstacles*(N+1)*N/mv_fact)+1:...
%     n_states*(N+1)+(n_veh_safty_states*N)+(n_pred_len*N)+(no_of_obstacles*(N+1)*N/mv_fact)+N+1))',1,N+1)';
%%%% Update time as per the the calculated time step by MPC
t(mpciter+1)=t0;
%%% Update initial position and the control trajectory calculated from MPC solution for
%%%%next mpc iteration
[t0, x0, u0]= shift_mv_obs(u(1,3), t0,x0,u, mv_obs_velo_ms, mv_obs_theta_rad, model_fn);
%%% recalculate acceleraation bounds which depend on current states
con_fn_value=con_fn(x0);
args.lbw(7: n_states : n_states*(N +1) ,1) = full(con_fn_value(2));
args.ubw(7: n_states : n_states*(N +1) ,1) = full(con_fn_value(1));
%%%%% Store the current calculated states    
xx(:,mpciter+2)=x0;
%%%%% Store the current calculated control actions    
ux(:,mpciter+2)=u(1,:);

%%%% Store the current solution of constraints
%%%%% Wheel liftoffconstraint
% veh_sfty_N(:,mpciter+2)=veh_sfty(1,:);
% %%%%% Lidar limit
% veh_dist(:,mpciter+1)=veh_lid(1,:);
% %%%%% MOving Obstacle abvoidance
% Mobs_dist(:,mpciter+1)=obs_dist_constr(1,:);
% %%%%% MOving Obstacle abvoidance
% Lobs_dist(:,mpciter+1)=obs_dist_lin(1,:);
%%%%Initialise the states trajectory calculated from MPC solution for
%%%%next mpc iteration
pred_state_mat_init = reshape ( full (sol.x(1:n_states*(N+1)))', n_states, N+1)';
pred_state_mat_init = [pred_state_mat_init(2:end,:);pred_state_mat_init(end,:)];
mpciter;
mpciter=mpciter+1
temp(mpciter)=norm((x0(1:2)-x_ref(1:2)),2);

end
main_loop_time=toc(main_loop);
ss_error= norm ((x0(1:2) - x_ref(1:2)),2)
average_mpc_time=main_loop_time/(mpciter+1)
Draw_MPC_Goal_stop_Mv_obs_param(t,xx,xx1,ux,x_ref,N,ref_obs,ref_traj,CG,lf,lr,w)
        
        
        
        
        
        

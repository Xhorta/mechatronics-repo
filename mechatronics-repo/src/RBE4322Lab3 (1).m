clc
clear

%% Assumptions
%link Material is 
% All links are rigid
% we can make the assumption that the scale is perfectly symmetrical
% Spring mass is negligible
% mass of the L-shape pin is included in M_rack
% spring is completely attached to the V Lever
% Bottom plate is assumed to be ground
% all levers are rectangular for MMoI
% L-piece that moves rack is rigidly attached to V Lever
% initial velocity is zero
% Using K = d^4*G/(8*D^3*n)
% the use of an equivalent force at each link for the calculation of Deq
% M_washer is the sum of the plastic washer and metal washer
%% Mass of each link
%converted mass from g to kg
g= 9.81;  %m/s^2
M_top= 0.5831;
M_lever = 0.02121;
M_Vlever = 0.1028;
M_rack = 0.009;
M_washer=.00383;
M_Lpiece= .001;

%% Distance
% Converted from cm to meter 
L1 = 0.132;
L2 = 0.105;
L3 = 0.245;
L4 = 0.009;
L5= 0.017;
L6= L1/2; 
L7 = L3*.25;


%Pinions Radius
R_pinion = 0.003;

% Washer radius
R_washer= .0159;

%% Dimensions of Springs

d_Spring1 = 0.001;
d_Spring3 = .0025;
d_Spring4 = .0003;

D_Spring1 = .008;
D_Spring3 = .018;
D_Spring4 = .006;

L_Spring1 = .01;
L_Spring3 = .016;
L_Spring4 = .0166;

%% MMoI
%Parallel axis theorem 
% Washer rotational inertia about the pinion axis
J_washer = 1/2*(M_washer*(R_washer^2));

%% MMoI (Kg*m^2)
MMoI_Vlever = 1/3*(M_Vlever*(L3^2));
MMoI_lever = 1/3*(M_lever*(L1^2));
MMoI_pinion = 3.15e-8;
MMoI_PaperWheel=2.012e-5;
MMoI_Lpiece = (1/12)*(M_Lpiece)*(L5^2+L4^2)+(M_Lpiece)*((sqrt(L4^2+L5^2)/2)^2); %Derived by treating L Piece as a rectangular plate
%Assumes that L-Piece is a perfect rectangle

%% ________________Spring information______________________________________
% stiffness calculations in (N/m)

%% No. Coils of Springs
N1 = 9;
N2 = 3.75;
N3 = 40;

%% Dimensions of Springs
d_k1 = 0.001;
d_k2 = .0025;
d_k3 = .0003;

D_k1 = .008;
D_k2 = .018;
D_k3 = .006;

L_k1 = .01;
L_k2 = .016;
L_k3 = .0166;

%% Spring stiffness calculations using k = d^4G/(8D^3n)
G_phosphorBronze = 43.1e9;
G_carbonSteel = 79.3e9;
G_musicWire = 81.7e9;

%% Define Springs stiffness for levers k1_S and k1_V (in N/m)
% because the top and bottom plates are grounded along with the rack pins
% and K1_Oringal connects the top and bottom plate, we can assume the K1_Original is acting on the levers as well
% we assume that K1_Original is split into two perfectly identical springs (k1/2),(k1/2)
% are placed the same distance along there respected lever pin
% spring setup of the external spring forces acting on our small levers(k1_S) and V-lever (k1_V)

%% Assuming springs k1_S and k1_V material is Phosphor Bronze (ASTM B159)
%% K1_Oringal => k1/2, k1/2
% Having a shear modulus of G=43.1 GPa, we multiply G*1000 = 43100 MPa, converting gigapascal to megapascal

%K1_Oringal = 1707.28/2; % N/m
K1_Oringal=((d_k1)^4*G_phosphorBronze)/(8*(D_k1)^3*N1);

%spring acting on small lever
k1_S = K1_Oringal/2; 

%spring acting on V-lever
k1_V = K1_Oringal/2; 

%% spring for Mass-plate
% Assuming the spring ,k2, material is made out of carbon steel 79.3 GPa
% Having a shear modulus of G=79.3 GPa, we multiply G*1000 = 79300 MPa, converting gigapascal to megapascal

k2=((d_k2)^4*G_carbonSteel)/(8*(D_k2)^3*N2);
%k2 10101.9 N/m.


%% spring for Rack 
% Assuming the spring k3 material is made out of music wire (ASTM A228) 
% Having a shear modulus of G=81.7 GPa, we multiply G*1000 = 81700 MPa, converting gigapascal to megapascal
k3=((d_k3)^4*G_musicWire)/(8*(D_k3)^3*N3);
%k3 14.61 % N/m.

%% Meq Calculations
% Mx''+ Dx' + Kx = 0 

Meq1 = MMoI_lever/(L1^2);
Jeq1 = (Meq1*(L2^2)) + MMoI_Vlever;
Meq2 = (Jeq1/(L3^2))*2 + M_top;
Jeq2 = Meq2*(L4^2)+ MMoI_Lpiece;
Meq3 = (Jeq2 /(L5^2))+ M_rack;


Meq = (Meq3*(R_pinion^2)) + J_washer + MMoI_pinion + MMoI_PaperWheel; % pinion+washer+paper wheel


%% Keq calculations
Keq1_R = k1_S*(L6^2);
Keq2_T = Keq1_R/(L1^2);
Keq3_R = (Keq2_T*(L2^2))+ (k1_V*(L7^2));
Keq4_T = (Keq3_R/(L3^2))*(2);
Keq5_T = Keq4_T + k2;
Keq6_R = (L4^2)*Keq5_T;
Keq7_T= (Keq6_R/(L5^2)) ;
Keq8_T= Keq7_T+k3;

%pinion
Keq = Keq8_T*R_pinion^2;  

%% Deq calculations
%Xavier's Body weight in kg
Xavier = 74.827;
% Convert body weight from kg to newtons
W= Xavier*g;

% Assuming the steels coefficient of friction u=0.6 at lever pivots and
% u=.4 at top plate and rack
% Assuming friction force is greatest at initial contact
% the use of an equivalent force at each link for the calculation of Deq

% Assuming that normal force in each link equal to W/4

% D1 = coef_friction * Force_normal * 4 * 1/pi
% assuming input frequency and amplitude of oscillation to be 1


D_lever= 4*(1/pi)*0.6*(W/4)*L1^2;   % small lever pivot 
D_Vlever= 4*(1/pi)*0.6*(W/4)*L1^2;   % V-lever pivot 
D_top= 4*(1/pi)*0.4*45;             % top plate contact 
D_rack= 4*(1/pi)*0.4*19;             % rack rail

Deq1_T= D_lever/(L1^2);                  % small lever -> lever tip 
Deq2_R = (Deq1_T*(L2^2)) + D_Vlever;        % -> V-lever rotation 
Deq3_T= (Deq2_R/(L3^2))*2 + D_top;         % -> V-lever ground / top plate, symmetrical 
Deq4_R= Deq3_T*(L4^2);                     % -> Lpin  
Deq5_T= (Deq4_R/(L5^2)) + D_rack;          % -> rack 

Deq = Deq5_T*(R_pinion^2);                % at pinion 

%% Feq calculations
%body weight acts on the top plate, so it is applied going down on M_top

Feq1 = W;                               

% @ L-pin 
Feq2_R= Feq1*L4;   

% @ rack
Feq3_T =Feq2_R/L5;                         

%%Torque @ pinion
Feq= Feq3_T*R_pinion;                  

%% Plot
% x'=diff(x,t)
syms x(t)

%first equation + switching the coefficients to rotational units since we are lumping to a rotational component
Jeq = Meq;
Beq = Deq;
K_Teq = Keq;
Teq = Feq;

eqn1 = Jeq*diff(x,t,2) + Beq*diff(x,t) + K_Teq*x == Teq;

%defining variable
Dx = diff(x,t);

%defining initial conditions
initialCondition = [x(0)==0,Dx(0)==0];

% Solve
sol = dsolve(eqn1,initialCondition);

%% Display lumped parameters for report table
fprintf('\nLumped parameters at pinion:\n');
fprintf('Jeq = %.4e kg*m^2\n', Jeq);
fprintf('Beq = %.4e N*m*s/rad\n', Beq);
fprintf('Keq = %.4e N*m/rad\n', K_Teq);
fprintf('Teq = %.4e N*m\n', Teq);

%natural frequency 
oscilation_freq= sqrt(K_Teq/Jeq);

% Damping ratio 
D_ratio = Beq/(2*sqrt(K_Teq*Jeq));

%Steady state angle for needle position 
ss_theta = Teq/K_Teq;

%% J·θ'' + B·θ' + K·θ = T
fprintf('oscilation_freq= %.2f rad/s, D_ratio= %.3f, ss_theta= %.2f rad (%.0f deg)\n', ...
    oscilation_freq, D_ratio, ss_theta, rad2deg(ss_theta));


%% Plots

% Angular Displacement Vs Time Graph
figure;
fplot(sol, [0 1]);
title("Angular Displacement vs Time");
xlabel("Time (s)");
ylabel("Angular Displacement (rad)");
legend('\theta(t)')     

%Angular Velocity Vs Time Graph
figure;
fplot(diff(sol, t), [0 1]);
title("Angular Velocity vs Time");
xlabel("Time (s)");
ylabel("Angular Velocity (rad/s)");
legend('\omega(t)')      


%Angular Acceleration Vs Time Graph
figure;
fplot(diff(diff(sol, t)),[0 1]);
title 'Angular Acceleration vs time';
xlabel 'time (s)';
ylabel 'Angular Acceleration (rad/s^2)';



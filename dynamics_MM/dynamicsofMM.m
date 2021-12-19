clc
clear all
close all
%% Introduction
% this code is for modelling the coupling dynamics of a mobile robot manipulator
% Zhengxue Zhou
% email: zhouzx@mpe.au.dk
%% Mobile Robot
%Symbolic Definition
%m_c is the mass of the base/platform with wheels
%m_t is the mass of each wheel
%r is the radius of each wheel
%I_mc the MoI of each wheel about the wheel diameter
%I_c the MoI of vhechile without wheels and the rotors of the motors about
%a vetical axis through the middle point between two wheels
%psi the orientation angle of the mobile robot measured from the x-axis
%q_r is the angular of right wheel
%e_x distance from mobile robot frame to the base of arm along X
%e_y distance from mobile robot frame to the base of arm along Y
%e_z distance from mobile robot frame to the base of arm along Z
%L The half distance between the actuated wheels of the mobile robot
syms t I_c m_c m_t I_mc r L dq_x dq_y psi dpsi q_r q_l dq_r dq_l
m_m = m_c+2*m_t;
I = I_c+2*m_t*L^2+2*I_mc;
psi=r/2/L*(q_r-q_l);
S_m=[r/2*cos(psi),r/2*cos(psi);
     r/2*sin(psi),r/2*sin(psi);
     r/2/L,-r/2/L;
     1,0;
     0,1];

T_b=1/2*m_m*(dq_x^2+dq_y^2)+(r/2)^2*I*(dq_r-dq_l)^2;
%% parameters for calculating Jacobian of the manipulator
syms q_x q_y x_c y_c e_x e_y e_z;

p_c_m=[x_c;y_c;0];%CM of the mobile robot.
p_0_x=[q_x;0;0]+p_c_m;%the first P joint along x
p_0_y=[q_x;q_y;0]+p_c_m;%the second P joint along y
R_0_c=[cos(psi),-sin(psi),p_0_x(1,1);sin(psi),cos(psi),p_0_y(2,1);0,0,1];
p_0_c=p_0_y+R_0_c*p_c_m;%the third R joint along z

p_0_b=[e_x;e_y;e_z]; %position of the base of arm under mobile base coordinate
T_0_b=[1,0,0,e_x;0,1,0,e_y;0,0,1,e_z;0,0,0,1];%from first joint to mobile base, assume there no rotation between  the two frames

R_b_w=[cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];%rotation matrix from base of arm to world
p_b_w=p_0_b+[q_x;q_y;0];%inertial position of the base under world coordinate
T_b_w=[cos(psi),-sin(psi),0,p_b_w(1,1);sin(psi),cos(psi),0,p_b_w(2,1);0,0,1,0;0,0,0,1];%from mobile base to word, translation fristly then rotation

%% Defination of the generalized coordinate
syms q_1 q_2 q_3 q_4 q_5 q_6 ...
    dq_1 dq_2 dq_3 dq_4 dq_5 dq_6;
q=[q_x;q_y;q_r;q_l;q_1;q_2;q_3;q_4;q_5;q_6];
dq=[dq_x;dq_y;dq_r;dq_l;dq_1;dq_2;dq_3;dq_4;dq_5;dq_6];
%% Defination of the DH and dynamic Parameters of UR5e
syms In_1x In_1y In_1z In_2x In_2y In_2z ...
     In_3x In_3y In_3z In_4x In_4y In_4z ...
     In_5x In_5y In_5z In_6x In_6y In_6z ...
     alpha_0 alpha_1 alpha_2 alpha_3 alpha_4 alpha_5 ...
     a_0 a_1 a_2 a_3 a_4 a_5 d_1 d_2 d_3 d_4 d_5 d_6 ...
     p_cx1 p_cy1 p_cz1 p_cx2 p_cy2 p_cz2 p_cx3 p_cy3 p_cz3 ...
     p_cx4 p_cy4 p_cz4 p_cx5 p_cy5 p_cz5 p_cx6 p_cy6 p_cz6 ...
     m_1 m_2 m_3 m_4 m_5 m_6 g;
In_1 = [In_1x,0,0;0,In_1y,0;0,0,In_1z];%MoI of CoM of each joint,coordinate of MoI is the Center of rotation. 
In_2 = [In_2x,0,0;0,In_2y,0;0,0,In_2z];
In_3 = [In_3x,0,0;0,In_3y,0;0,0,In_3z];
In_4 = [In_4x,0,0;0,In_4y,0;0,0,In_4z];
In_5 = [In_5x,0,0;0,In_5y,0;0,0,In_5z];
In_6 = [In_6x,0,0;0,In_6y,0;0,0,In_6z];
%% General transformation between link1 to 6
%ROTATION MATRICEs 
R_1=[cos(q_1) -sin(q_1) 0;
     sin(q_1)*cos(alpha_0) cos(q_1)*cos(alpha_0) -sin(alpha_0);
     sin(q_1)*sin(alpha_0) cos(q_1)*sin(alpha_0)  cos(alpha_0)];
R_2=[cos(q_2) -sin(q_2) 0;
     sin(q_2)*cos(alpha_1) cos(q_1)*cos(alpha_1) -sin(alpha_1);
     sin(q_2)*sin(alpha_1) cos(q_2)*sin(alpha_1)  cos(alpha_1)];
R_3=[cos(q_3) -sin(q_3) 0;
     sin(q_3)*cos(alpha_2) cos(q_3)*cos(alpha_2) -sin(alpha_2);
     sin(q_3)*sin(alpha_2) cos(q_3)*sin(alpha_2)  cos(alpha_2)];
R_4=[cos(q_4) -sin(q_4) 0;
     sin(q_4)*cos(alpha_3) cos(q_4)*cos(alpha_3) -sin(alpha_3);
     sin(q_4)*sin(alpha_3) cos(q_4)*sin(alpha_3)  cos(alpha_3)];
R_5=[cos(q_5) -sin(q_5) 0;
     sin(q_5)*cos(alpha_4) cos(q_5)*cos(alpha_4) -sin(alpha_4);
     sin(q_5)*sin(alpha_4) cos(q_5)*sin(alpha_4)  cos(alpha_4)];
R_6=[cos(q_6) -sin(q_6) 0;
     sin(q_6)*cos(alpha_5) cos(q_6)*cos(alpha_5) -sin(alpha_5);
     sin(q_6)*sin(alpha_5) cos(q_6)*sin(alpha_5)  cos(alpha_5)];
% POSITION VECTORs
p_1=[a_0;-sin(alpha_0)*d_1;cos(alpha_0)*d_1];
p_2=[a_1;-sin(alpha_1)*d_2;cos(alpha_1)*d_2];
p_3=[a_2;-sin(alpha_2)*d_3;cos(alpha_2)*d_3];
p_4=[a_3;-sin(alpha_3)*d_4;cos(alpha_3)*d_4];
p_5=[a_4;-sin(alpha_4)*d_5;cos(alpha_4)*d_5];
p_6=[a_5;-sin(alpha_5)*d_6;cos(alpha_5)*d_6];
%% TRANSLATION MATRICES AND FORWARD KINEMATICS
T_1 = [R_1,p_1;zeros(1,3),1];
T_2 = [R_2,p_2;zeros(1,3),1];
T_3 = [R_3,p_3;zeros(1,3),1];
T_4 = [R_4,p_4;zeros(1,3),1];
T_5 = [R_5,p_5;zeros(1,3),1];
T_6 = [R_6,p_6;zeros(1,3),1];
T = T_1*T_2*T_3*T_4*T_5*T_6;%based on base of arm frame
T_6_w=T_b_w*T_0_b*T_1*T_2*T_3*T_4*T_5*T_6;
%% ROTATION MATRICEs FROM world
R_1=R_b_w*R_1;
R_20=R_1*R_2;
R_30=R_20*R_3;
R_40=R_30*R_4;
R_50=R_40*R_5;
R_60=R_50*R_6;
%% COMs' POSITION VECTORs based on base arm frame
p_c1=p_1+R_1*[p_cx1;p_cy1;p_cz1];
p_c2=p_1+R_1*(p_2+R_2*[p_cx2;p_cy2;p_cz2]);
p_c3=p_1+R_1*(p_2+R_2*(p_3+R_3*[p_cx3;p_cy3;p_cz3]));
p_c4=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*[p_cx4;p_cy4;p_cz4])));
p_c5=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*(p_5+R_5*([p_cx5;p_cy5;p_cz5])))));
p_c6=p_1+R_1*(p_2+R_2*(p_3+R_3*(p_4+R_4*(p_5+R_5*(p_6+R_6*[p_cx6;p_cy6;p_cz6])))));
% from base of arm to world
p_c1=p_b_w+R_b_w*(p_0_b+p_c1);
p_c2=p_b_w+R_b_w*(p_0_b+p_c2);
p_c3=p_b_w+R_b_w*(p_0_b+p_c3);
p_c4=p_b_w+R_b_w*(p_0_b+p_c4);
p_c5=p_b_w+R_b_w*(p_0_b+p_c5);
p_c6=p_b_w+R_b_w*(p_0_b+p_c6);
%% LINEAR and ANGULAR PART of JACOBIANs
ee = [T_6_w(1,4);T_6_w(2,4);T_6_w(3,4)];%position of the EE
J_v = jacobian(ee,q);
J_o=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),R_50(:,3),R_60(:,3)];
Jacob = [J_v;J_o];
%% pseudo-inertia matrices linear part
J_vx=jacobian(p_0_x,q);% Jacobian of center of mass,
J_vy=jacobian(p_0_y,q);
J_v1=jacobian(p_c1,q);
J_v2=jacobian(p_c2,q);
J_v3=jacobian(p_c3,q);
J_v4=jacobian(p_c4,q);
J_v5=jacobian(p_c5,q);
J_v6=jacobian(p_c6,q);

%% pseudo-inertia matrices angular part
J_oz=[zeros(3,2),r/2/L,-r/2/L,zeros(3,6)];
J_o1=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),zeros(3,5)];
J_o2=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),zeros(3,4)];
J_o3=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),R_30(:,3),zeros(3,3)];
J_o4=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),zeros(3,2)];
J_o5=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),R_50(:,3),zeros(3,1)];
J_o6=[zeros(3,2),r/2/L,-r/2/L,R_1(:,3),R_20(:,3),R_30(:,3),R_40(:,3),R_50(:,3),R_60(:,3)];
%% elements of the Kinetic energy
JJ_v1=J_v1.'*m_1*eye(3)*J_v1;
JJ_o1=J_o1.'*R_1*In_1*R_1.'*J_o1;
JJ_v2=J_v2.'*m_2*eye(3)*J_v2;
JJ_o2=J_o2.'*R_20*In_2*R_20.'*J_o2;
JJ_v3=J_v3.'*m_3*eye(3)*J_v3;
JJ_o3=J_o3.'*R_30*In_3*R_30.'*J_o3;
JJ_v4=J_v4.'*m_4*eye(3)*J_v4;
JJ_o4=J_o4.'*R_40*In_4*R_40.'*J_o4;
JJ_v5=J_v5.'*m_5*eye(3)*J_v5;
JJ_o5=J_o5.'*R_50*In_5*R_50.'*J_o5;
JJ_v6=J_v6.'*m_6*eye(3)*J_v6;
JJ_o6=J_o6.'*R_60*In_6*R_60.'*J_o6;
%% Kinetic energy
M_a = JJ_v1+JJ_o1+JJ_v2+JJ_o2+JJ_v3+JJ_o3+JJ_v4+JJ_o4+JJ_v5+JJ_o5+JJ_v6+JJ_o6;
K_a = 0.5*dq.'*JJ_v1*dq+dq.'*JJ_o1*dq+ ...,
 dq.'*JJ_v2*dq+dq.'*JJ_o2*dq+ ...,
 dq.'*JJ_v3*dq+dq.'*JJ_o3*dq+ ...,
 dq.'*JJ_v4*dq+dq.'*JJ_o4*dq+ ...,
 dq.'*JJ_v5*dq+dq.'*JJ_o5*dq+ ...,
 dq.'*JJ_v6*dq+dq.'*JJ_o6*dq;
%% Kinetic energy
K=K_a+T_b;
%% ROBOT's INERTIA (MASS) MATRIX 

M_d=zeros(10,10);
M_d=[diff(K,dq_x).';diff(K,dq_y).';diff(K,dq_r).';diff(K,dq_l).';diff(K,dq_1).';diff(K,dq_2).';diff(K,dq_3).';diff(K,dq_4).';diff(K,dq_5).';diff(K,dq_6).'];
M=[diff(M_d(1,1),dq_x),diff(M_d(1,1),dq_y),diff(M_d(1,1),dq_r),diff(M_d(1,1),dq_l),diff(M_d(1,1),dq_1),diff(M_d(1,1),dq_2),diff(M_d(1,1),dq_3),diff(M_d(1,1),dq_4),diff(M_d(1,1),dq_5),diff(M_d(1,1),dq_6); ...
   diff(M_d(2,1),dq_x),diff(M_d(2,1),dq_y),diff(M_d(2,1),dq_r),diff(M_d(2,1),dq_l),diff(M_d(2,1),dq_1),diff(M_d(2,1),dq_2),diff(M_d(2,1),dq_3),diff(M_d(2,1),dq_4),diff(M_d(2,1),dq_5),diff(M_d(2,1),dq_6); ...
   diff(M_d(3,1),dq_x),diff(M_d(3,1),dq_y),diff(M_d(3,1),dq_r),diff(M_d(3,1),dq_l),diff(M_d(3,1),dq_1),diff(M_d(3,1),dq_2),diff(M_d(3,1),dq_3),diff(M_d(3,1),dq_4),diff(M_d(3,1),dq_5),diff(M_d(3,1),dq_6); ...
   diff(M_d(4,1),dq_x),diff(M_d(4,1),dq_y),diff(M_d(4,1),dq_r),diff(M_d(4,1),dq_l),diff(M_d(4,1),dq_1),diff(M_d(4,1),dq_2),diff(M_d(4,1),dq_3),diff(M_d(4,1),dq_4),diff(M_d(4,1),dq_5),diff(M_d(4,1),dq_6); ...
   diff(M_d(5,1),dq_x),diff(M_d(5,1),dq_y),diff(M_d(5,1),dq_r),diff(M_d(5,1),dq_l),diff(M_d(5,1),dq_1),diff(M_d(5,1),dq_2),diff(M_d(5,1),dq_3),diff(M_d(5,1),dq_4),diff(M_d(5,1),dq_5),diff(M_d(5,1),dq_6); ...
   diff(M_d(6,1),dq_x),diff(M_d(6,1),dq_y),diff(M_d(6,1),dq_r),diff(M_d(6,1),dq_l),diff(M_d(6,1),dq_1),diff(M_d(6,1),dq_2),diff(M_d(6,1),dq_3),diff(M_d(6,1),dq_4),diff(M_d(6,1),dq_5),diff(M_d(6,1),dq_6); ...
   diff(M_d(7,1),dq_x),diff(M_d(7,1),dq_y),diff(M_d(7,1),dq_r),diff(M_d(7,1),dq_l),diff(M_d(7,1),dq_1),diff(M_d(7,1),dq_2),diff(M_d(7,1),dq_3),diff(M_d(7,1),dq_4),diff(M_d(7,1),dq_5),diff(M_d(7,1),dq_6); ...
   diff(M_d(8,1),dq_x),diff(M_d(8,1),dq_y),diff(M_d(8,1),dq_r),diff(M_d(8,1),dq_l),diff(M_d(8,1),dq_1),diff(M_d(8,1),dq_2),diff(M_d(8,1),dq_3),diff(M_d(8,1),dq_4),diff(M_d(8,1),dq_5),diff(M_d(8,1),dq_6); ...
   diff(M_d(9,1),dq_x),diff(M_d(9,1),dq_y),diff(M_d(9,1),dq_r),diff(M_d(9,1),dq_l),diff(M_d(9,1),dq_1),diff(M_d(9,1),dq_2),diff(M_d(9,1),dq_3),diff(M_d(9,1),dq_4),diff(M_d(9,1),dq_5),diff(M_d(9,1),dq_6); ...
   diff(M_d(10,1),dq_x),diff(M_d(10,1),dq_y),diff(M_d(10,1),dq_r),diff(M_d(10,1),dq_l),diff(M_d(10,1),dq_1),diff(M_d(10,1),dq_2),diff(M_d(10,1),dq_3),diff(M_d(10,1),dq_4),diff(M_d(10,1),dq_5),diff(M_d(10,1),dq_6)];
M_n=SSS'*M*SSS;
%% ROBOT's CORIOLIS and CENTRIFUGAL MATRIX

for k=1:10
   for s=1:10
      C(k,s)=.5*((diff(M(k,s),q_x)+diff(M(k,1),q(s,1))-diff(M(1,s),q(k,1)))*dq_x...
                +(diff(M(k,s),q_y)+diff(M(k,2),q(s,1))-diff(M(2,s),q(k,1)))*dq_y...
                +(diff(M(k,s),q_r)+diff(M(k,3),q(s,1))-diff(M(3,s),q(k,1)))*dq_r...
                +(diff(M(k,s),q_l)+diff(M(k,4),q(s,1))-diff(M(4,s),q(k,1)))*dq_l...
                +(diff(M(k,s),q_1)+diff(M(k,5),q(s,1))-diff(M(5,s),q(k,1)))*dq_1...
                +(diff(M(k,s),q_2)+diff(M(k,6),q(s,1))-diff(M(6,s),q(k,1)))*dq_2...
                +(diff(M(k,s),q_3)+diff(M(k,7),q(s,1))-diff(M(7,s),q(k,1)))*dq_3...
                +(diff(M(k,s),q_4)+diff(M(k,8),q(s,1))-diff(M(8,s),q(k,1)))*dq_4...
                +(diff(M(k,s),q_5)+diff(M(k,9),q(s,1))-diff(M(5,s),q(k,1)))*dq_5...
                +(diff(M(k,s),q_6)+diff(M(k,10),q(s,1))-diff(M(10,s),q(k,1)))*dq_6);
   end
end

fid = fopen('C16.txt', 'w');
fwrite(fid, char(vvv), 'char');
fclose(fid);
%% POTENTIAL ENERGIES and GRAVITY VECTOR
P1=m_1*[0,0,g]*p_c1;
P2=m_2*[0,0,g]*p_c2;
P3=m_3*[0,0,g]*p_c3;
P4=m_4*[0,0,g]*p_c4;
P5=m_5*[0,0,g]*p_c5;
P6=m_6*[0,0,g]*p_c6;
P=P2+P3+P4+P5+P6;
g_x=diff(P,q_x);
g_y=diff(P,q_y);
g_r=diff(P,q_r);
g_l=diff(P,q_l);
g_1=diff(P,q_1);
g_2=diff(P,q_2);
g_3=diff(P,q_3);
g_4=diff(P,q_4);
g_5=diff(P,q_5);
g_6=diff(P,q_6);
G=[g_x;g_y;g_r;g_l;g_1;g_2;g_3;g_4;g_5;g_6];

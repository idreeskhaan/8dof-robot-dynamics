clc
close all

syms  theta d a alpha;
syms q1 q2 q3 q4 q5 q6 q7 q8; 
syms L1 L2 L3 L4 L5 L6 L7 L8;
%%syms d4    %to create symbolic variables

b=100;

T =[ cos(theta)   -cos(alpha)*sin(theta)   sin(alpha)*sin(theta)   a*cos(theta)
    sin(theta)    cos(alpha)*cos(theta)   -cos(theta)*sin(alpha)   a*sin(theta)
    0             sin(alpha)               cos(alpha)              d
    0             0                        0                       1           ];

%% Transformation Matrices for All Links

%% Link 1

T1=subs(T,{theta, d, a, alpha,},{q1,0,0,pi/2});
%% Link 2
% 

syms b;
T2=subs(T,{theta, d, a, alpha,},{q2,0,b,-pi/2});
%% Link 3

T3=subs(T,{theta, d, a, alpha,},{q3,0,0,pi/2});
%% Link 4

T4=subs(T,{theta, d, a, alpha,},{q4,0,b,-pi/2});
%% Link 5

T5=subs(T,{theta, d, a, alpha,},{q5,0,0,pi/2});
%% Link 6

T6=subs(T,{theta, d, a, alpha,},{q6,0,b,-pi/2});
%% Link 7

T7=subs(T,{theta, d, a, alpha,},{q7,0,0,pi/2});
%% Link 8

T8=subs(T,{theta, d, a, alpha,},{q8,0,b,-pi/2});
%%


% Transformation Matrices for each link wrt Base Frame {0}
T01=T1;
T02= simplify(T1*T2);
T03= simplify(T02*T3);
T04= simplify(T03*T4);
T05= simplify(T04*T5);
T06= simplify(T05*T6);
T07= simplify(T06*T7);
T08= simplify(T07*T8);

%% Position of Center of masses of each link
PC1= T01*[0;0;0;1];
PC2= T02*[-b/2;0;0;1];
PC3= T03*[0;0;0;1];
PC4= T04*[-b/2;0;0;1];
PC5= T05*[0;0;0;1];
PC6= T06*[-b/2;0;0;1];
PC7= T07*[0;0;0;1];
PC8= T08*[-b/2;0;0;1];

%define the gravity direction based on the WORLD FRAME {0}
gx=0;
gy=-9.81;
gz=0;

%Potential Energies for each link
syms m1 m2 m3 m4 m5 m6 m7 m8;
PE1= -m1*[gx gy gz]*PC1(1:3);
PE2= -m2*[gx gy gz]*PC2(1:3);
PE3= -m3*[gx gy gz]*PC3(1:3);
PE4= -m4*[gx gy gz]*PC4(1:3);
PE5= -m5*[gx gy gz]*PC5(1:3);
PE6= -m6*[gx gy gz]*PC6(1:3);
PE7= -m7*[gx gy gz]*PC7(1:3);
PE8= -m8*[gx gy gz]*PC8(1:3);
PE= PE1+PE2+PE3+PE4+PE5+PE6+PE7+PE8; %Total PE

%Velocities[Linear] of center of masses of each link [Chain Rule]
syms q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q7_dot q8_dot;
VC1= diff(PC1,q1)*q1_dot + diff(PC1,q2)*q2_dot + diff(PC1,q3)*q3_dot + diff(PC1,q4)*q4_dot + diff(PC1,q5)*q5_dot + diff(PC1,q6)*q6_dot + diff(PC1,q7)*q7_dot + diff(PC1,q8)*q8_dot; %fun of q1
VC2= diff(PC2,q1)*q1_dot + diff(PC2,q2)*q2_dot + diff(PC2,q3)*q3_dot + diff(PC2,q4)*q4_dot + diff(PC2,q5)*q5_dot + diff(PC2,q6)*q6_dot + diff(PC2,q7)*q7_dot + diff(PC2,q8)*q8_dot; %fun of q1 q2
VC3= diff(PC3,q1)*q1_dot + diff(PC3,q2)*q2_dot + diff(PC3,q3)*q3_dot + diff(PC3,q4)*q4_dot + diff(PC3,q5)*q5_dot + diff(PC3,q6)*q6_dot + diff(PC3,q7)*q7_dot + diff(PC3,q8)*q8_dot; %fun of q1 q2
VC4= diff(PC4,q1)*q1_dot + diff(PC4,q2)*q2_dot + diff(PC4,q3)*q3_dot + diff(PC4,q4)*q4_dot + diff(PC4,q5)*q5_dot + diff(PC4,q6)*q6_dot + diff(PC4,q7)*q7_dot + diff(PC4,q8)*q8_dot; %fun of q1 q2 q3 q4
VC5= diff(PC5,q1)*q1_dot + diff(PC5,q2)*q2_dot + diff(PC5,q3)*q3_dot + diff(PC5,q4)*q4_dot + diff(PC5,q5)*q5_dot + diff(PC5,q6)*q6_dot + diff(PC5,q7)*q7_dot + diff(PC5,q8)*q8_dot; %fun of q1 q2 q3 q4
VC6= diff(PC6,q1)*q1_dot + diff(PC6,q2)*q2_dot + diff(PC6,q3)*q3_dot + diff(PC6,q4)*q4_dot + diff(PC6,q5)*q5_dot + diff(PC6,q6)*q6_dot + diff(PC6,q7)*q7_dot + diff(PC6,q8)*q8_dot; %fun of q1 q2 q3 q4 q5 q6
VC7= diff(PC7,q1)*q1_dot + diff(PC7,q2)*q2_dot + diff(PC7,q3)*q3_dot + diff(PC7,q4)*q4_dot + diff(PC7,q5)*q5_dot + diff(PC7,q6)*q6_dot + diff(PC7,q7)*q7_dot + diff(PC7,q8)*q8_dot; %fun of q1 q2 q3 q4 q6 q6
VC8= diff(PC8,q1)*q1_dot + diff(PC8,q2)*q2_dot + diff(PC8,q3)*q3_dot + diff(PC8,q4)*q4_dot + diff(PC8,q5)*q5_dot + diff(PC8,q6)*q6_dot + diff(PC8,q7)*q7_dot + diff(PC8,q8)*q8_dot; %fun of q1 q2 q3 q4 q6 q6 q7 q8


%Angular Velocity Propagation
%First calculate the Rotation Matrices out of Transf Matrices
R21= transpose(T2(1:3,1:3));
R32= transpose(T3(1:3,1:3));
R43= transpose(T4(1:3,1:3));
R54= transpose(T5(1:3,1:3));
R65= transpose(T6(1:3,1:3));
R76= transpose(T7(1:3,1:3));
R87= transpose(T8(1:3,1:3));

%Now Apply iterative formula for Angular Velocities
w11=[0; 0; q1_dot];
w22= R21*w11+ [0; 0; q2_dot];
w33= R32*w22+ [0; 0; q3_dot];
w44= R43*w33+ [0; 0; q4_dot];
w55= R54*w44+ [0; 0; q5_dot];
w66= R65*w55+ [0; 0; q6_dot];
w77= R76*w66+ [0; 0; q7_dot];
w88= R87*w77+ [0; 0; q8_dot];

%Define Inertia Matrices
syms Ixx1 Iyy1 Izz1;
syms Ixx2 Iyy2 Izz2;
syms Ixx3 Iyy3 Izz3;
syms Ixx4 Iyy4 Izz4;
syms Ixx5 Iyy5 Izz5;
syms Ixx6 Iyy6 Izz6;
syms Ixx7 Iyy7 Izz7;
syms Ixx8 Iyy8 Izz8;

I1= [Ixx1, 0, 0;
    0, Iyy1, 0;
    0, 0, Izz1];
I2= [Ixx2, 0, 0;
    0, Iyy2, 0;
    0, 0, Izz2];
I3= [Ixx3, 0, 0;
    0, Iyy3, 0;
    0, 0, Izz3];
I4= [Ixx4, 0, 0;
    0, Iyy4, 0;
    0, 0, Izz4];
I5= [Ixx5, 0, 0;
    0, Iyy5, 0;
    0, 0, Izz5];
I6= [Ixx6, 0, 0;
    0, Iyy6, 0;
    0, 0, Izz6];
I7= [Ixx7, 0, 0;
    0, Iyy7, 0;
    0, 0, Izz7];
I8= [Ixx8, 0, 0;
    0, Iyy8, 0;
    0, 0, Izz8];

%Kinetic Energies
KE1= 0.5*m1*transpose(VC1(1:3))*VC1(1:3) + 0.5*transpose(w11)*I1*w11;
KE2= 0.5*m2*transpose(VC2(1:3))*VC2(1:3) + 0.5*transpose(w22)*I2*w22;
KE3= 0.5*m3*transpose(VC3(1:3))*VC3(1:3) + 0.5*transpose(w33)*I3*w33;
KE4= 0.5*m4*transpose(VC4(1:3))*VC4(1:3) + 0.5*transpose(w44)*I4*w44;
KE5= 0.5*m5*transpose(VC5(1:3))*VC5(1:3) + 0.5*transpose(w55)*I5*w55;
KE6= 0.5*m6*transpose(VC6(1:3))*VC6(1:3) + 0.5*transpose(w66)*I6*w66;
KE7= 0.5*m7*transpose(VC7(1:3))*VC7(1:3) + 0.5*transpose(w77)*I7*w77;
KE8= 0.5*m8*transpose(VC8(1:3))*VC8(1:3) + 0.5*transpose(w88)*I8*w88;
KE= KE1+KE2+KE3+KE4+KE5+KE6+KE7+KE8; %total KE


%Lagrangian Dynamics
%Define angular accelerations for each link
syms q1_2dot q2_2dot q3_2dot q4_2dot q5_2dot q6_2dot q7_2dot q8_2dot; 

%For tau1
temp1_q1= diff(KE, q1_dot);
temp2_q1 = diff(temp1_q1, q1).*q1_dot + diff(temp1_q1, q2).*q2_dot + diff(temp1_q1, q3).*q3_dot + diff(temp1_q1, q4).*q4_dot + diff(temp1_q1, q5).*q5_dot + diff(temp1_q1, q6).*q6_dot + diff(temp1_q1, q7).*q7_dot +diff(temp1_q1, q8).*q8_dot + diff(temp1_q1, q1_dot).*q1_2dot + diff(temp1_q1, q2_dot).*q2_2dot + diff(temp1_q1, q3_dot).*q3_2dot + diff(temp1_q1, q4_dot).*q4_2dot + diff(temp1_q1, q5_dot).*q5_2dot + diff(temp1_q1, q6_dot).*q6_2dot+ diff(temp1_q1, q7_dot).*q7_2dot+ diff(temp1_q1, q8_dot).*q8_2dot;
temp3_q1 = - diff(KE, q1);
temp4_q1 = diff(PE, q1) ;
tau1= temp2_q1 +temp3_q1 +temp4_q1;

%For tau2
temp1_q2= diff(KE, q2_dot);
temp2_q2 = diff(temp1_q2, q1).*q1_dot + diff(temp1_q2, q2).*q2_dot + diff(temp1_q2, q3).*q3_dot + diff(temp1_q2, q4).*q4_dot + diff(temp1_q2, q5).*q5_dot + diff(temp1_q2, q6).*q6_dot + diff(temp1_q2, q7).*q7_dot +diff(temp1_q2, q8).*q8_dot + diff(temp1_q2, q1_dot).*q1_2dot + diff(temp1_q2, q2_dot).*q2_2dot + diff(temp1_q2, q3_dot).*q3_2dot + diff(temp1_q2, q4_dot).*q4_2dot + diff(temp1_q2, q5_dot).*q5_2dot + diff(temp1_q2, q6_dot).*q6_2dot+ diff(temp1_q2, q7_dot).*q7_2dot+ diff(temp1_q2, q8_dot).*q8_2dot;
temp3_q2 = - diff(KE, q2);
temp4_q2 = diff(PE, q2) ;
tau2= temp2_q2 +temp3_q2 +temp4_q2;

%For tau3
temp1_q3= diff(KE, q3_dot);
temp2_q3 = diff(temp1_q3, q1).*q1_dot + diff(temp1_q3, q2).*q2_dot + diff(temp1_q3, q3).*q3_dot + diff(temp1_q3, q4).*q4_dot + diff(temp1_q3, q5).*q5_dot + diff(temp1_q3, q6).*q6_dot + diff(temp1_q3, q7).*q7_dot +diff(temp1_q3, q8).*q8_dot + diff(temp1_q3, q1_dot).*q1_2dot + diff(temp1_q3, q2_dot).*q2_2dot + diff(temp1_q3, q3_dot).*q3_2dot + diff(temp1_q3, q4_dot).*q4_2dot + diff(temp1_q3, q5_dot).*q5_2dot + diff(temp1_q3, q6_dot).*q6_2dot+ diff(temp1_q3, q7_dot).*q7_2dot+ diff(temp1_q3, q8_dot).*q8_2dot;
temp3_q3 = - diff(KE, q3);
temp4_q3 = diff(PE, q3) ;
tau3= temp2_q3 +temp3_q3 +temp4_q3;

%For tau4
temp1_q4= diff(KE, q4_dot);
temp2_q4 = diff(temp1_q4, q1).*q1_dot + diff(temp1_q4, q2).*q2_dot + diff(temp1_q4, q3).*q3_dot + diff(temp1_q4, q4).*q4_dot + diff(temp1_q4, q5).*q5_dot + diff(temp1_q4, q6).*q6_dot + diff(temp1_q4, q7).*q7_dot +diff(temp1_q4, q8).*q8_dot + diff(temp1_q4, q1_dot).*q1_2dot + diff(temp1_q4, q2_dot).*q2_2dot + diff(temp1_q4, q3_dot).*q3_2dot + diff(temp1_q4, q4_dot).*q4_2dot + diff(temp1_q4, q5_dot).*q5_2dot + diff(temp1_q4, q6_dot).*q6_2dot+ diff(temp1_q4, q7_dot).*q7_2dot+ diff(temp1_q4, q8_dot).*q8_2dot;
temp3_q4 = - diff(KE, q4);
temp4_q4 = diff(PE, q4) ;
tau4= temp2_q4 +temp3_q4 +temp4_q4;

%For tau5
temp1_q5= diff(KE, q5_dot);
temp2_q5 = diff(temp1_q5, q1).*q1_dot + diff(temp1_q5, q2).*q2_dot + diff(temp1_q5, q3).*q3_dot + diff(temp1_q5, q4).*q4_dot + diff(temp1_q5, q5).*q5_dot + diff(temp1_q5, q6).*q6_dot + diff(temp1_q5, q7).*q7_dot +diff(temp1_q5, q8).*q8_dot + diff(temp1_q5, q1_dot).*q1_2dot + diff(temp1_q5, q2_dot).*q2_2dot + diff(temp1_q5, q3_dot).*q3_2dot + diff(temp1_q5, q4_dot).*q4_2dot + diff(temp1_q5, q5_dot).*q5_2dot + diff(temp1_q5, q6_dot).*q6_2dot+ diff(temp1_q5, q7_dot).*q7_2dot+ diff(temp1_q5, q8_dot).*q8_2dot;
temp3_q5 = - diff(KE, q5);
temp4_q5 = diff(PE, q5) ;
tau5= temp2_q5 +temp3_q5 +temp4_q5;

%For tau6
temp1_q6= diff(KE, q6_dot);
temp2_q6 = diff(temp1_q6, q1).*q1_dot + diff(temp1_q6, q2).*q2_dot + diff(temp1_q6, q3).*q3_dot + diff(temp1_q6, q4).*q4_dot + diff(temp1_q6, q5).*q5_dot + diff(temp1_q6, q6).*q6_dot + diff(temp1_q6, q7).*q7_dot +diff(temp1_q6, q8).*q8_dot + diff(temp1_q6, q1_dot).*q1_2dot + diff(temp1_q6, q2_dot).*q2_2dot + diff(temp1_q6, q3_dot).*q3_2dot + diff(temp1_q6, q4_dot).*q4_2dot + diff(temp1_q6, q5_dot).*q5_2dot + diff(temp1_q6, q6_dot).*q6_2dot+ diff(temp1_q6, q7_dot).*q7_2dot+ diff(temp1_q6, q8_dot).*q8_2dot;
temp3_q6 = - diff(KE, q6);
temp4_q6 = diff(PE, q6) ;
tau6= temp2_q6 +temp3_q6 +temp4_q6;

%For tau7
temp1_q7= diff(KE, q7_dot);
temp2_q7 = diff(temp1_q7, q1).*q1_dot + diff(temp1_q7, q2).*q2_dot + diff(temp1_q7, q3).*q3_dot + diff(temp1_q7, q4).*q4_dot + diff(temp1_q7, q5).*q5_dot + diff(temp1_q7, q6).*q6_dot + diff(temp1_q7, q7).*q7_dot +diff(temp1_q7, q8).*q8_dot + diff(temp1_q7, q1_dot).*q1_2dot + diff(temp1_q7, q2_dot).*q2_2dot + diff(temp1_q7, q3_dot).*q3_2dot + diff(temp1_q7, q4_dot).*q4_2dot + diff(temp1_q7, q5_dot).*q5_2dot + diff(temp1_q7, q6_dot).*q6_2dot+ diff(temp1_q7, q7_dot).*q7_2dot+ diff(temp1_q7, q8_dot).*q8_2dot;
temp3_q7 = - diff(KE, q7);
temp4_q7 = diff(PE, q7) ;
tau7= temp2_q7 +temp3_q7 +temp4_q7;

%For tau8
temp1_q8= diff(KE, q8_dot);
temp2_q8 = diff(temp1_q8, q1).*q1_dot + diff(temp1_q8, q2).*q2_dot + diff(temp1_q8, q3).*q3_dot + diff(temp1_q8, q4).*q4_dot + diff(temp1_q8, q5).*q5_dot + diff(temp1_q8, q6).*q6_dot + diff(temp1_q8, q7).*q7_dot +diff(temp1_q8, q8).*q8_dot + diff(temp1_q8, q1_dot).*q1_2dot + diff(temp1_q8, q2_dot).*q2_2dot + diff(temp1_q8, q3_dot).*q3_2dot + diff(temp1_q8, q4_dot).*q4_2dot + diff(temp1_q8, q5_dot).*q5_2dot + diff(temp1_q8, q6_dot).*q6_2dot+ diff(temp1_q8, q7_dot).*q7_2dot+ diff(temp1_q8, q8_dot).*q8_2dot;
temp3_q8 = - diff(KE, q8);
temp4_q8 = diff(PE, q8) ;
tau8= temp2_q8 +temp3_q8 +temp4_q8;

%Now substitute some constant values
q1=0;q2=0;q3=0;q4=0;q5=0;q6=0;q7=0;q8=0;
q1_dot=0;q2_dot=0;q3_dot=0;q4_dot=0;q5_dot=0;q6_dot=0;q7_dot=0;q8_dot=0;
q1_2dot=0;q2_2dot=0;q3_2dot=0;q4_2dot=0;q5_2dot=0;q6_2dot=0;q7_2dot=0;q8_2dot=0;


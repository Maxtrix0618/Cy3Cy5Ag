% 【单染料-银膜耦合 耦合参数计算】
% 输入纯银等离激元峰位、纯染料的左右峰位和极化子色散三峰位（除中）共5个波长，计算输出2个耦合参数

clear; clc;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!

L=[512.5 555];          % 纯染料峰波长 [左 右]
La=[575	480	525	610];       % 纯银峰位和极化子色散的三个峰位波长 [纯银 左 中 右]

E=He./L;
S=He./La(1);
Ea=He./La(2:4);

syms D1 D2;
H=[
     S  , D1 , D2 ;
     D1 ,E(1),  0 ;
     D2 ,  0 ,E(2)
];
I=eye(3);

Eq1 = sym(det(H-(Ea(1)*I))==0);
Eq2 = sym(det(H-(Ea(3)*I))==0);

Sol = solve(Eq1,Eq2);

disp('Matlab求解结果')
disp('      D1        D2      ');
disp([double(Sol.D1),double(Sol.D2)]);
disp('参考（绝对值）')
disp('      D1        D2      ');
disp([double(abs(Sol.D1(1))),double(abs(Sol.D2(1)))]);



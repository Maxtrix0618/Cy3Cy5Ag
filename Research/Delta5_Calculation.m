% 【双染料-银膜耦合 耦合参数计算】
% 输入纯银等离激元峰位、两种纯染料的左右峰位和极化子色散五峰位（除中）共9个波长，计算输出4个耦合参数

clear; clc;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!

L=[525 565 605 660];            % 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
La=[749	508	545	603	640	751];	% 纯银峰位和极化子色散的五个峰位波长 [纯银 左2 左1 中 右1 右2]

E=He./L;
S=He./La(1);
Ea=He./La(2:6);

syms D1 D2 D3 D4;
H=[
     S  , D1 , D2 , D3 , D4 ;
     D1 ,E(1),  0 ,  0 ,  0 ;
     D2 ,  0 ,E(2),  0 ,  0 ;
     D3 ,  0 ,  0 ,E(3),  0 ;
     D4 ,  0 ,  0 ,  0 ,E(4);
];
I = eye(5);

Eq1 = sym(det(H-(Ea(1)*I))==0);
Eq2 = sym(det(H-(Ea(2)*I))==0);
Eq3 = sym(det(H-(Ea(4)*I))==0);
Eq4 = sym(det(H-(Ea(5)*I))==0);

Sol = solve(Eq1,Eq2,Eq3,Eq4);

disp('Matlab求解结果')
disp('      D1        D2        D3        D4      ');
disp([double(Sol.D1),double(Sol.D2),double(Sol.D3),double(Sol.D4)]);
disp('参考（绝对值）')
disp('      D1        D2        D3        D4      ');
disp([double(abs(Sol.D1(1))),double(abs(Sol.D2(1))),double(abs(Sol.D3(1))),double(abs(Sol.D4(1)))]);



% 【双染料-银膜耦合 平均耦合参数计算】

clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
Im = eye(5);

fid=fopen('Exp5.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[6,inf]);
n=length(Exp(1,:));

DA = zeros(1,4);

for I=1:n

    L=zeros(1,10);
    L(1)=Exp(1,I);      % 纯银峰位波长
    L(2)=525;           % Cy3左峰波长
    L(3)=565;           % Cy3右峰波长
    L(4)=605;           % Cy5左峰波长
    L(5)=660;           % Cy5右峰波长

    L(6)=Exp(2,I);     % 极化子色散左2峰
    L(7)=Exp(3,I);     % 极化子色散左1峰
    L(8)=Exp(4,I);     % 极化子色散中峰（不使用）
    L(9)=Exp(5,I);     % 极化子色散右1峰
    L(10)=Exp(6,I);    % 极化子色散右2峰

    E=He./L;

    syms D1 D2 D3 D4;
    H=[
        E(1), D1 , D2 , D3 , D4 ;
         D1 ,E(2),  0 ,  0 ,  0 ;
         D2 ,  0 ,E(3),  0 ,  0 ;
         D3 ,  0 ,  0 ,E(4),  0 ;
         D4 ,  0 ,  0 ,  0 ,E(5);
    ];

    Eq1 = sym(det(H-(E(6)*Im))==0);
    Eq2 = sym(det(H-(E(7)*Im))==0);
    Eq3 = sym(det(H-(E(8)*Im))==0);
    Eq4 = sym(det(H-(E(9)*Im))==0);

    Sol = solve(Eq1,Eq2,Eq3,Eq4);

    disp(['组',num2str(I)]);
    disp('      D1        D2        D3        D4      ');
    disp([double(Sol.D1), double(Sol.D2), double(Sol.D3), double(Sol.D4)]);
    
    DA(1) = DA(1) + abs(Sol.D1(1));
    DA(2) = DA(2) + abs(Sol.D2(1));
    DA(3) = DA(3) + abs(Sol.D3(1));
    DA(4) = DA(4) + abs(Sol.D4(1));

end

for i=1:4
    DA(i) = DA(i)./n;
end
disp('绝对值均值');
disp('      D1        D2        D3        D4      ');
disp([DA(1),DA(2),DA(3),DA(4)]);





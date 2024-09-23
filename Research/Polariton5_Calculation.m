% 【双染料-银膜耦合 极化子色散峰计算】
% 输入两种纯染料的两峰波长和耦合参数，计算输出以各波长为等离激元峰的银膜同这些染料的耦合所产生的5个极化子色散峰峰位波长

clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda） 

L=[525 565 605 660];	% 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
E=He./L;                % 各波长对应能量

D=[0.16139      0.26151      0.15529      0.35152];   % 耦合参数Delta_1~4

p = 1;          % 步长(nm)
Lp=450:p:850;   % 横轴：波长450~850(nm)
n=(400/p)+1;    % 离散点组数目
Eg=zeros(n,5);  % n行5列的2D数组，记录第行数个波长的5个矩阵特征值


% 遍历范围内所有波长
for i=1:n
    S=He/Lp(i);  % 当前波长能量（作为纯银膜峰位波长）
    A=[
    S   ,D(1),D(2),D(3),D(4);
    D(1),E(1),  0 ,  0 ,  0 ;
    D(2),  0 ,E(2),  0 ,  0 ;
    D(3),  0 ,  0 ,E(3),  0 ;
    D(4),  0 ,  0 ,  0 ,E(4);
    ];
    Eg(i,:) = eig(A);   % 求A的（5个）特征值（能量）
end

Lac = He ./ Eg;    % 将特征能量转换到特征波长（数组）

% 打印
disp('     NSF       Le2       Le1       Mid       Ri1       Ri2      ');
disp([Lp(:),Lac(:,5),Lac(:,4),Lac(:,3),Lac(:,2),Lac(:,1)]);

% 作图
figure
axis([450,850 450,850]);
darkGreen = [4 157 107]/255; lw = 1.5;
for i=1:4
    plot([450,850],[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on;            % 染料4峰参考线（Cy3左右、Cy5左右）
end
plot(Lp,Lp, 'b-. ');                                                                        % 纯银膜参考线
plot(Lp,Lac(:,1),'r- ',Lp,Lac(:,2),'r- ',Lp,Lac(:,3),'r- ',Lp,Lac(:,4),'r- ',Lp,Lac(:,5),'r- ');	% 极化子色散峰（计算结果）


fid=fopen('Exp5.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[6,inf]);

for i=2:6
    plot(Exp(1,:),Exp(i,:),'k+ '); hold on;   % 离散点（实验数据）
end

title('Cy3-Cy5@Ag 五峰实验数据与计算')
xlabel('Bare plasmon peak(nm)');
ylabel('Polariton peak(nm)')


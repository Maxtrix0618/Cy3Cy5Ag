% 【双染料-银膜耦合 染料交叉项调试(仅C24)-变染料浓度比】
% 考虑染料峰C24相互作用，生成理论值曲线，两条对照

clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda） 

L=[523 561 621 665];	% 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
E=He./L;                % 各波长对应能量

S=He/610;

KD3=0.6;        % Cy3变化调比
KD5=0.8;        % Cy5变化调比
D=[0 0 0 0];    % 耦合参数Delta_1~4

p=0.01;         % 步长
Lp=0:p:1;       % 横轴：u[Cy3]/(u[Cy3]+u[Cy5])
n=1/p+1;        % 离散点组数目
CV=[0, 0.5];    % 染料交叉项取值
Eg=zeros(2,n,5);    % 2*n*5列的2D数组，记录在一个染料交叉值下第行数个Dye比值的5个矩阵特征值

% 染料交叉项
C12=0;
C13=0;
C14=0;
C23=0;
C24=0;
C34=0;

% 染料交叉项取值遍历
for ci=1:2
    C24=CV(ci);
    for xi=1:n
        x=xi-1;
        D(1)=KD3*sqrt(x*p);
        D(2)=KD3*sqrt(x*p);
        D(3)=KD5*sqrt(1-x*p);
        D(4)=KD5*sqrt(1-x*p);
        A=[
            S   ,D(1),D(2),D(3),D(4);
            D(1),E(1),C12 ,C13 ,C14 ;
            D(2),C12 ,E(2),C23 ,C24 ;
            D(3),C13 ,C23 ,E(3),C34 ;
            D(4),C14 ,C24 ,C34 ,E(4);
        ];
        Eg(ci,xi,:) = sort(eig(A),'descend');   % 求A的（5个）特征值（能量）
    end
end
Lac = He ./ Eg;     % 将特征能量转换到特征波长（数组）


% 作图
figure
set(gcf, 'Position', [0, 0, 700, 600]); 
Set = struct('color',{{'k','r','b','g'}},'lineStyle', {{'--', '-'}});
for lc=1:2
    for li=1:3
        plot(Lp,Lac(lc,:,li),Set.lineStyle{lc},'color',Set.color{li}); hold on;
    end
    plot(Lp,Lac(lc,:,4),Set.lineStyle{lc},'color',[4 157 107]/255); hold on;
    plot(Lp,Lac(lc,:,5),Set.lineStyle{lc},'color',[192 0 237]/255); hold on;
end
xlabel('u[Cy3] / (u[Cy3]+u[Cy5])');
ylabel('Polariton peak(nm)');
title(['Cy3-Cy5@Ag 染料交叉项:C24']);


% 输出数据
data=zeros(n,11);
data(:,1)=Lp;
for lc=1:2
    for li=1:5
%         data(:,(li-1)*2+lc+1)=Lac(lc,:,li);
        data(:,(lc-1)*5+li+1)=Lac(lc,:,li);
    end
end
save('DyeDL_C24.txt', 'data', '-ascii');

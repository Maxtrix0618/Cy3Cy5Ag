% 【单染料-银膜耦合 手动调整耦合参数D】
% 对同一批次的每个样品分别手动输入D

clear; clc;

I=eye(3);
darkGreen = [4 157 107]/255; lw = 1.5;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
L=[512.5 555];	% 纯染料峰波长 [左 右]
E=He./L;

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';
N=length(Exp(:,1));

DS=zeros(N,2);
DS(1,:)=[0.2,0.2];
DS(2,:)=[0.2,0.2];
DS(3,:)=[0.2,0.2];
DS(4,:)=[0.2,0.2];
DS(5,:)=[0.2,0.2];
DS(6,:)=[0.2,0.2];
DS(7,:)=[0.2,0.2];

CalS=zeros(N,3);
for EI=1:N
    S=He/Exp(EI,1);
    D=DS(EI,:);
    A=[
        S   ,D(1),D(2);
        D(1),E(1),  0 ;
        D(2),  0 ,E(2);
      ];
    CalS(EI,:) = He ./ sort(eig(A),'descend');
end


% 作图
figure
Range=[450,850];
Lp=Range(1):1:Range(2);
axis([Range(1),Range(2) Range(1),Range(2)]);
for i=1:2
    plot(Range,[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
end
plot(Lp,Lp, 'b-. ');    hold on;

% 实验数据
for i=1:3
    plot(Exp(:,1),CalS(:,i),'r^ '); hold on;
    plot(Exp(:,1),Exp(:,i+1),'ko '); hold on;
end

title(['三峰计算与实验数据'])
xlabel('Bare plasmon peak(nm)');
ylabel('Polariton peak(nm)');





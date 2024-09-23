% 【双染料-银膜耦合 自动调整耦合参数D】
% 对同一批次的每个样品分别自动调整计算D以求令极化峰计算值拟合实验数据

clear; clc;

I=eye(5);
darkGreen = [4 157 107]/255; lw = 1.5;

He=1239.841;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）
L=[525 565 605 660];	% 纯染料峰波长 [左 右]
E=He./L;

fid=fopen('Exp5.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[6,inf]);
Exp=Exp';
N=length(Exp(:,1));

DS=0.1.*ones(N,4);
% DS(7,:)=0.1.*ones(1,4);
% DS(9,:)=0.1.*ones(1,4);

p=0.0001;       % 拟合调整步长
d=1;            % 拟合调整方向
DIJS=[2,3,4,1]; % 四元倒反矢


CalS=zeros(N,5);
for EI=1:N
    S=He/Exp(EI,1);
    D=DS(EI,:);
    disp(['>>开始组',num2str(EI)]);
    ErrLast=10000;  % 上次迭代误差，这里随意给个大数
    gts2=0;         % 连续顺行次数标：连续顺行2次才说明未打转
    spin=0;         % 拟合值打转次数标：调头一次+1，顺行一次归0，不小于2说明拟合程序达到底点
    spinBig=0;      % 大打转次数标：在两个耦合系数之间轮转，原理同上
    IOfD=2;         % 现在调整的D指标（1,2,3,4）
    while spinBig<4
        A=[
            S   ,D(1),D(2),D(3),D(4);
            D(1),E(1),  0 ,  0 ,  0 ;
            D(2),  0 ,E(2),  0 ,  0 ;
            D(3),  0 ,  0 ,E(3),  0 ;
            D(4),  0 ,  0 ,  0 ,E(4);
          ];
        Cal=He ./ sort(eig(A),'descend');
        ErrCurr = (Cal(1)-Exp(EI,2))^2+(Cal(2)-Exp(EI,3))^2+(Cal(3)-Exp(EI,4))^2+(Cal(4)-Exp(EI,5))^2+(Cal(5)-Exp(EI,6))^2;
        disp(['Err',num2str(ErrCurr)]);

        if ErrCurr<ErrLast
            % 顺行
            disp('->')
            if gts2>=2
                gts2=0;
                spin=0;
                spinBig=0;
            else
                gts2=gts2+1;
            end
        else
            % 调头
            disp('O')
            gts2=0;
            spin=spin+1;
            d=-d;
        end
        if spin<4
            D(IOfD)=D(IOfD)+(d*p);
        else
            spin=0;
            spinBig=spinBig+1;
            IOfD=DIJS(IOfD);
            disp(['>组',num2str(EI),'现在调整：D',num2str(IOfD)]);
        end
        ErrLast=ErrCurr;
    end
    disp(['!组',num2str(EI),'调整完成：[D1 D2 D3 D4]=[',num2str(D(1)),' ',num2str(D(2)),' ',num2str(D(3)),' ',num2str(D(4)),']']);
    DS(EI,:)=D;
    CalS(EI,:) = Cal;
end

disp('NSF    D1    D2    D3    D4');
for i=1:N
    disp([num2str(Exp(i,1)),' ',num2str(DS(i,1)),' ',num2str(DS(i,2)),' ',num2str(DS(i,3)),' ',num2str(DS(i,4))]);
end

% 作图
figure
Range=[450,850];
Lp=Range(1):1:Range(2);
axis([Range(1),Range(2) Range(1),Range(2)]);
for i=1:4
    plot(Range,[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
end
plot(Lp,Lp, 'b-. ');    hold on;

% 实验数据
for i=1:5
    plot(Exp(:,1),CalS(:,i),'r^ '); hold on;
    plot(Exp(:,1),Exp(:,i+1),'ko '); hold on;
end

title('Cy3-Cy5-Ag')
xlabel('Bare plasmon peak(nm)');
ylabel('Polariton peak(nm)');





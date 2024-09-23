% 【自动调整简并耦合参数D 单染料-银膜耦合】

clear; clc;

I=eye(3);
darkGreen = [4 157 107]/255; lw = 1.5;

He=1239.841;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）
L=[498 538];	% 纯染料峰波长 [左 右]
E=He./L;

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';
N=length(Exp(:,1));

DS=0.2.*ones(N);

p=0.0001;       % 拟合调整步长
d=1;            % 拟合调整方向

CalS=zeros(N,3);
for EI=1:N
    S=He/Exp(EI,1);
    D=DS(EI);
    disp(['>>开始组',num2str(EI)]);
    ErrLast=10000;  % 上次迭代误差，这里随意给个大数
    gts2=0;         % 连续顺行次数标：连续顺行2次才说明未打转
    spin=0;         % 拟合值打转次数标：调头一次+1，顺行一次归0，不小于2说明拟合程序达到底点
    DO=1;
    while DO==1
        A=[
            S  , D  , D  ;
            D  ,E(1), 0  ;
            D  , 0  ,E(2);
          ];
        Cal=He ./ sort(eig(A),'descend');
        ErrCurr = (Cal(1)-Exp(EI,2))^2+(Cal(2)-Exp(EI,3))^2+(Cal(3)-Exp(EI,4))^2;
        disp(['Err ',num2str(ErrCurr)]);

        if ErrCurr<ErrLast
            % 顺行
            disp('->')
            disp(gts2);
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
            D=D+(d*p);
        else
            DO=0;
        end
        ErrLast=ErrCurr;
    end
    disp(['!组',num2str(EI),'调整完成：D=',num2str(D)]);
    DS(EI)=D;
    CalS(EI,:) = Cal;
end

disp('NSF    D');
for i=1:N
    disp([num2str(Exp(i,1)),' ',num2str(DS(i,1))]);
end

% 作图
figure
Range=[450,750];
Lp=Range(1):1:Range(2);
axis([Range(1) Range(2), Range(1) Range(2)]);
for i=1:2
    plot(Range,[L(i),L(i)], ':','color',darkGreen,'linewidth',lw);  hold on
end
plot(Lp,Lp, 'b-. ');    hold on;

% 实验数据
for i=1:3
    plot(Exp(:,1),CalS(:,i),'r^ '); hold on;
    plot(Exp(:,1),Exp(:,i+1),'ko '); hold on;
end

title('RhB-Ag')
xlabel('Bare plasmon peak(nm)');
ylabel('Polariton peak(nm)');





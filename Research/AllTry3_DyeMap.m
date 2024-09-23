% 【单染料-染料微调地图（耦合系数手动输入：应为同批次实验数据均值）】

clear; clc;

He=1243.125;	% 波长和能量换算常数 （nm->eV, E=hc/lamda）  :H=1239.841 ?!
D=[0.28164  0.10271];

fid=fopen('Exp3.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[4,inf]);
Exp=Exp';
n=length(Exp(:,1));
I=eye(3);

O=[512.5 555];	% 初始纯染料峰波长 [左 右]
p=5;  R=40;    % 微调步长与单向步数
N=2*R+1;
Lo=[O(1)-p*R, O(1)+p*R];
Ro=[O(2)-p*R, O(2)+p*R];
LeRange = Lo(1) : p : Lo(2);
RiRange = Ro(1) : p : Ro(2);

ErrMap=zeros(N,N,3);

for Le=LeRange
    for Ri=RiRange
        E=He./[Le,Ri];

        Eg=zeros(n,4);
        for i=1:n
            S=He/Exp(i,1);
            A=[
            S   ,D(1),D(2);
            D(1),E(1),  0 ;
            D(2),  0 ,E(2);
            ];
            Eg(i,2:4) = sort(eig(A),'descend');
        end
        Lac = He ./ Eg;
        Lac(:,1)=Exp(:,1);

        % 给出三条计算峰线与实验数据的标准差
        le = (Le-Lo(1))/p +1;
        ri = (Ri-Ro(1))/p +1;
        for line=1:3
            for i=1:n
                ErrMap(le,ri,line) = ErrMap(le,ri,line) + (Lac(i,line+1)-Exp(i,line+1))^2;
            end
        end
        
%         if le>ri
%             ErrMap(le,ri,1)=0; ErrMap(le,ri,2)=0; ErrMap(le,ri,3)=0;
%         end
        
    end
end

ErrMap = ErrMap./n;
ErrMap = sqrt(ErrMap);

% 创建坐标网格
x=zeros(N,N);
y=zeros(N,N);
for i=1:N
    for j=1:N
        x(i,j)=p*i+Lo(1);
        y(i,j)=p*j+Ro(1);
    end
end

Mid=floor((N+1)/2);
for line=1:3
    mesh(x,y,ErrMap(:,:,line)); hold on;
    disp([num2str(line),'  ',num2str(ErrMap(Mid,Mid,line))])
end
plot3([O(1),O(1)],[O(2),O(2)],[0,40],'-.','Color',[0.4,0.4,0.4]);

title('微调染料峰：极化峰实验与计算的标准差')
xlabel('染料左峰');
ylabel('染料右峰');




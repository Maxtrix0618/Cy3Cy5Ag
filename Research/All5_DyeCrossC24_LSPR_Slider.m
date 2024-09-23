function All5_DyeCrossC24_LSPR_Slider
clear; clc;

He=1243.125; % 波长和能量换算常数 （nm->eV, E=hc/lamda） 

L=[523 561 621 665];	% 纯染料峰波长 [Cy3左 Cy3右 Cy5左 Cy5右]
E=He./L;                % 各波长对应能量

S=He/610;

KD1=@(VK1) VK1; VK1=0;	% D1变调比
KD3=@(VK3) VK3; VK3=0;	% D3变调比
CV2=@(VC) VC; VC=0.5;   % C24取值
U=@(Vu) Vu; Vu=0.5;     % u3/uT浓度比值

D=[0 0 0 0];    % 耦合参数Delta_1~4

p = 1;          % 步长(nm)
Lp=450:p:850;   % 横轴：波长450~850(nm)
n=(400/p)+1;    % 离散点组数目


% 染料交叉项
C12=0; C13=0; C14=0; C23=0 ;C24=0; C34=0;

fid=fopen('Exp5.txt','r');       % 从文件读取实验数据
Exp=fscanf(fid,'%f',[6,inf]);
Exp=Exp';

F = figure(1);
set(gcf, 'Position', [800, 200, 700, 700]); 
LineSet = struct('lineStyle', {{'--', '-'}});
ColorSet = [[0 0 0];[255 0 0];[0 0 255];[24 157 137];[192 0 237]]/255;

Hu = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.14 0.75 0.03],...
    'min' , 0,...
    'max' , 1,...
    'value' , Vu,...
    'callback', @(s, e) @sliderCallback);
H1 = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.10 0.75 0.03],...
    'min' , 0,...
    'max' , 1,...
    'value' , VK1,...
    'callback', @(s, e) @sliderCallback);
H3 = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.06 0.75 0.03],...
    'min' , 0,...
    'max' , 1,...
    'value' , VK3,...
    'callback', @(s, e) @sliderCallback);
HC = uicontrol(...
    'parent' , F,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.2 0.02 0.75 0.03],...
    'min' , 0,...
    'max' , 2,...
    'value' , VC,...
    'callback', @(s, e) @sliderCallback);

addlistener(Hu, 'Value', 'PostSet', @sliderCallback);
addlistener(H1, 'Value', 'PostSet', @sliderCallback);
addlistener(H3, 'Value', 'PostSet', @sliderCallback);
addlistener(HC, 'Value', 'PostSet', @sliderCallback);

Tu = annotation('textbox',[0.05 0.14 0.8 0.04],'String',['u3/uT= ', num2str(sprintf('%.2f',Vu))],'EdgeColor','none');
T1 = annotation('textbox',[0.05 0.10 0.8 0.04],'String',['KD1= ', num2str(sprintf('%.2f',VK1))],'EdgeColor','none');
T3 = annotation('textbox',[0.05 0.06 0.8 0.04],'String',['KD3= ', num2str(sprintf('%.2f',VK3))],'EdgeColor','none');
TC = annotation('textbox',[0.05 0.02 0.8 0.04],'String',['C24= ', num2str(sprintf('%.3f',VC))],'EdgeColor','none');

P = [];
for lcc=1:2
    for lii=1:5
        P(lcc,lii) = plot(0:0,0:0); hold on;
    end
end
set(gca, 'position', [0.1 0.25 0.85 0.7]);

function sliderCallback(~,~)
    for lc=1:2
        for li=1:5
            delete(P(lc,li));
        end
    end
    delete(Tu); delete(T1); delete(T3); delete(TC);
    
    U = get(Hu, 'value');
    KD1 = get(H1, 'value');
    KD3 = get(H3, 'value');
    CV2 = get(HC, 'value');
    CVs=[0 CV2];
    
    Tu = annotation('textbox',[0.05 0.14 0.8 0.04],'String',['u3/uT= ', num2str(sprintf('%.3f',U))],'EdgeColor','none');
    T1 = annotation('textbox',[0.05 0.10 0.8 0.04],'String',['KD1= ', num2str(sprintf('%.2f',KD1))],'EdgeColor','none');
    T3 = annotation('textbox',[0.05 0.06 0.8 0.04],'String',['KD3= ', num2str(sprintf('%.2f',KD3))],'EdgeColor','none');
    TC = annotation('textbox',[0.05 0.02 0.8 0.04],'String',['C24= ', num2str(sprintf('%.3f',CV2))],'EdgeColor','none');
    
    Eg=zeros(2,n,5);
    for ci=1:2
        C24=CVs(ci);
        for i=1:n
            S=He/Lp(i);
            D(1)=KD1*sqrt(U*p);
            D(2)=KD1*sqrt(U*p);
            D(3)=KD3*sqrt(1-U*p);
            D(4)=KD3*sqrt(1-U*p);
            A=[
                S   ,D(1),D(2),D(3),D(4);
                D(1),E(1),C12 ,C13 ,C14 ;
                D(2),C12 ,E(2),C23 ,C24 ;
                D(3),C13 ,C23 ,E(3),C34 ;
                D(4),C14 ,C24 ,C34 ,E(4);
            ];
            Eg(ci,i,:) = sort(eig(A),'descend');   % 求A的（5个）特征值（能量）
        end
    end

    for i=1:2
        for j=1:n
            for k=1:5
                if (imag(Eg(i,j,k))~=0) % 可能会出现虚数特征值，此处除去之
                    Eg(i,j,k) = 0;
                end
            end
        end
    end
    Lac = He ./ Eg;     % 将特征能量转换到特征波长（数组）
    
    % 实验数据
    for i=2:6
        plot(Exp(:,1),Exp(:,i),'+','color',ColorSet(i-1,:)); hold on;
    end
    
    for lc=1:2
        for li=1:5
            P(lc,li) = plot(Lp,Lac(lc,:,li),LineSet.lineStyle{lc},'color',ColorSet(li,:)); hold on;
        end
    end
    xlabel('Bare plasmon peak(nm)');
    ylabel('Polariton peak(nm)');
    title('Cy3-Cy5@Ag [C24]');
    
end

end


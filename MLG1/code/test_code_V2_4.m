% 2019/11/10 第一篇论文的代码原始来源
% 修改最后的GF层，原来是使用sum(abs(G_f))的，继续使用这个，一二次导数联立叠加
% v2.4 分块计算,但是需要自己定义每个块


clc; clear all; close all;tic;

% 图像的属性

% imageAddress='.\HM_L32.0_V142.5_H60\';  % fake【需要减1】
% %【1-5】
% imageAddress='.\20191020T182622\';  % 81  ---球速快
% imageAddress='.\20191029_ball_720_f18\';  % 132 ---球速正常
% imageAddress='.\20191103_1_ball_720_25\';  % 160
% imageAddress='.\20191103_1_car_720_25\';  % 190
% imageAddress='.\20191103_2_carcamera_720_25_v1\';  % 247
% %【6-10】
% imageAddress='.\20191103_2_carcamera_720_25_v2\';  % 82 【需要减1】----车速不一样
% imageAddress='.\20191103_2_carmoving_720_25\';  % 163
% imageAddress='.\20191105_ball_crosswise_far_720_25\';  % 150 横向
% imageAddress='.\20191105_ball_crosswise_med_720_25\';  % 122 横向
% imageAddress='.\20191105_ball_crosswise_near_720_25\';  % 114 横向
% %【11-14】
% imageAddress='.\20191110_car-ball_720_25\';  % 253，有背板，小球目标
% imageAddress='.\20191110_car-box_720_25\';  % 299，有背板，盒子目标
% imageAddress='.\20191110_car-box-NoBG_720_25\'; % 316, 无背板，盒子目标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imageAddress='.\NEW\PIC_T\V30_50\';     % tanslating
% imageAddress='.\PIC_FOR_PAPER\PIC_VEHICLE\V14_w_A\';     % approaching
% imageAddress='.\NEW\V7_B_L\';     % leave
imageAddress='.\tests images\4_leaving sample\plus\';


% %帧的总数 fake=60; 1103_2_carcamera_v2=82
% numFrames = 82;
numFrames = length(dir([imageAddress '*.jpg']));

Hor=720;  %视频帧的宽度 fake=270 real_182622=720
Ver=720; %视频帧的高度 fake=270 real_182622=720




% 将文件夹里所有的图像放在一个三维矩阵里,或者直接使用下面的代码读取已保存的文件
seriesImages=zeros(Ver,Hor,numFrames);     %声明一个三维矩阵
for k=1:numFrames
    %     FileName = strcat((imageAddress),'f',num2str(k-1),'.png');  % fake
    %     image = imread(FileName,'png');     % fake
    
    FileName = strcat((imageAddress),'f',num2str(k),'.jpg');  % real
    image = imread(FileName,'jpg');     % real
    
    seriesImages(:,:,k)=image;
    
end


% % 分块
j = -4;   % 范围 【-7~8】    需要针对【1-16标号】减【8】    
imageMatrix=imageDepart_v3(seriesImages,j);




% LGMD参数设置部分
Ver = size(imageMatrix,1);  %图片高度
Hor = size(imageMatrix,2);  %水平宽度

frame = size(imageMatrix,3);         %帧数


% The persistence of the luminance change
np = 1;
u = 1;

% local inhibition weight
wi = [0.125,0.25,0.125;0.25,0,0.25;0.125,0.25,0.125];   %侧抑制层，与论文中一致
% global inhibition weight
WI = 0.3;   %S层抑制权重，与论文中一致

% convolution mask (from S-->G) %G层均值滤波

we = ones(3,3).*1/9;



% Scale (from S-->G)
Cw = 4;     % G层参数，算w的
delta_c = 0.01;      % G层参数，算w的

% decay coefficient range (0,1) （G）
Cde = 0.5;      %衰减系数
% dacay threshold  （G）
Tde = 15;       %衰减阈值 15

% Flag of using fixed Ts, 1 for use and 0 for no
Flag_FixedTs = 0;
% Flag of using fixed Tffi, 1 for use and 0 for no
Flag_FixedTffi = 0;

if Flag_FixedTffi == 1
    % fixed FFI threshold
    Tffi = 80;
elseif Flag_FixedTffi == 0
    % Parameters in adapatable FFI threshold
    Tffi = 10;
    AFffi = 0.02;
    Tfo = 30;
end

if Flag_FixedTs == 0
    % fixed spiking threshold (No FFM cell)
    Ts = 0.7;
else
    % Parameters in FFM cell - No evidence in locust
    Ts = 0;
    AFlt = 4;
    AFmp = 1;
    Tmp = 0.9;
    AFl = 1;
    Tlto = 0;
    T_u = 230;
    T_l = 180;
    DetaTlt = 0.03;
    Ld = uint16(0);
end

% Creat the memory of data
% The current Image
Lf = uint8(zeros(Ver,Hor));
% The Image one frame delay
Lf_1 = uint8(zeros(Ver,Hor));
% The P layer
Pf = uint8(zeros(Ver,Hor));
Pf_1 = uint8(zeros(Ver,Hor));

p = zeros(1,10);
% Luminace persistence matrix, length defined by np and 10 is enough
% Actually, np is 1 in the paper
for i=1:10
    p(i) = 1/(1+exp(u*i));
end

% The E/I layer
Ef = uint8(zeros(Ver,Hor));
If = zeros(Ver,Hor);
% The S layer,note S may be negative
Sf = zeros(Ver,Hor);
% The G layer - No evidence in locust
Gf = zeros(Ver,Hor);    % 阈值判断后的Gf
G_f=zeros(Ver,Hor);     % 阈值判断前的Gf
Cef = zeros(Ver,Hor);   % Cef值
w = 0;
% LGMD cell,note that I set it as a 1*20 vector to stored it for abserving
Kf = zeros(1,frame);
kf = ones(1,frame)*0.5;


% integral kf
Gf_number=zeros(1,frame);   % 阈值滤波后的Gf矩阵中非0元素的个数

first_Intergral = zeros(1,frame);   % 一次导
second_Integral = zeros(1,frame);   % 二次导

Df = zeros(1,frame);                  % 等于 当前Gf_number减去上一个Gf_number
df = ones(1,frame)*0.5;            %

% kkf=zeros(1,frame);
ntimes = 1;






% The FFI cell
FFI = zeros(1,frame);
% The spiking
Spike = zeros(1,frame);
Spike2 = zeros(1,frame);

figure('color','w');
set(gcf,'outerposition',get(0,'screensize'));


for k = 1:frame
    % Make some of the data used in the last step be zero
    If = zeros(Ver,Hor);
    % tansmit the data to the Lf-1,Pf-1 from the Lf, Pf-1
    Lf_1 = double(Lf);
    Pf_1 = Pf;
    
    Lf=double(imageMatrix(:,:,k));
    
    % LGMD algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Pf = abs(Lf - Lf_1) + sum(p(1:np));     %这里有问题，少乘了Pf-1(x,y)
    Pf = abs(Lf - Lf_1) ;       %修改1：直接去掉过去几帧对当前帧影响部分 --output2
    %     Pf = abs(Lf - Lf_1) + sum(p(1:np).*Pf_1);     %修改2：乘以Pf-1(x,y)
    
    
    
    Ef = Pf;        %激励层等于P层
    
    If=If+conv2(double(Pf_1),double(wi),'same');    %侧抑制层
    
    % Note that I ignored the boundary,
    % i.e.,If(1,:)=If(:,1)=If(:,Hro)=If(Ver,:)=0
    
    Sf = double(Ef) - If.*WI;   %WI是常数=0.3      -- Sf output3
    
    
    Cef=conv2(double(Sf),double(we),'same');    %均值滤波       --output4
    w = delta_c + max(max(abs(Cef)))/Cw;        % delta_c是极小实数=0.01，Cw是常数=4
    G_f = Sf.*Cef/w;         %   -output5
    
    Gf=G_f;             % 数值传递
    
    for i = 1:Ver
        for j = 1:Hor
            if Gf(i,j)*Cde < Tde         % Cde=0.5，Tde=15
                
                Gf(i,j) = 0;                % --output6
            end
        end
    end
    
    
    Kf(1,k) = sum(sum(abs(Gf)));     % 小kf
    kf(1,k) = 1/(1+exp(-Kf(1,k)*3.87/(Ver*Hor)));    % sig函数         --output8
    
    % test code
    Gf_number(1,k) = Kf(1,k);
    if k>2
        first_Intergral(1,k) = (Gf_number(k)-Gf_number(k-2))/0.04;        % 一阶导 25帧率
        second_Integral(1,k) = (Gf_number(k-2)+Gf_number(k)-2*Gf_number(k-1))/0.0016;  % 二阶导
        
        if k == 38
            ntimes = 0.0001;
        end
        
        if first_Intergral(1,k) > 0 && second_Integral(1,k) > 0
            ntimes = ntimes + 0.3;
        elseif first_Intergral(1,k) > 0 && second_Integral(1,k) < 0
            ntimes = ntimes + 0.2;
        else
            ntimes = ntimes - 0.3;
        end
        
        if ntimes <= 0
            ntimes = 0.0001;
        end
        
        
        Df(k) = abs(Gf_number(k)-Gf_number(k-1));
        
        df(1,k) = 1/(1+exp((-Gf_number(1,k)*ntimes*3.87)/(Ver*Hor)));

        
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   FFI，转向抑制机制
    if Flag_FixedTffi == 0
        % adaptable Tffi
        Tffi = Tfo + AFffi*Tffi;
    end
       
    nnn = (sum(sum(abs(Pf_1)))*3.87 )/ (Ver*Hor)
    
    if (sum(sum(abs(Pf_1)))*3.87 )/ (Ver*Hor) > Tffi        
        FFI(1,k) = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   FFM
    if Flag_FixedTs == 1
        % adaptable Ts using FFM cell
        Ld = (sum(max(Lf)) + sum(max(Lf')))/(Ver + Hor);
        if Ld > T_u
            Tlt = Tlto + AFl*DetaTlt;
        else if Ld < T_l
                Tlt = Tlto - AFl*DetaTlt;
            else
                Tlt = 0;
            end
        end
        Ts = AFlt*Tlt + AFmp*Tmp;   % 在机器人上AFlt大于0，在电脑上应该等于0；  AFmp=1  Tmp=0.86
    end
    
    %   Spike!!
    if kf(1,k) >= Ts && FFI(1,k)==0 && k>3
        Spike(1,k) = 1;
    end
    
    if k>4 && (Spike(1,k-3)+Spike(1,k-2)+Spike(1,k-1)+Spike(1,k) >= 5)
        Spike2(1,k) = 1;
    end
    
    
    % 实时输出图像，并且随着一帧帧更新
    subplot(341);imshow(mat2gray(imageMatrix(:,:,k)));title('原图');
    
    subplot(342);imshow(Pf);title('P');
    
    subplot(343);imshow(Sf);title('S');
    
    subplot(344);imshow(Cef);title('G->Cef');
    
    subplot(345);imshow(G_f);title('G->G_f');
    
    subplot(346);imshow(Gf);title('Gf');
    
    subplot(347);
    plot(Kf);
    xlim([1,frame]);        % 必须是5 开始，否则前几帧数值太大
    title('Kf：Gf层非0元素绝对值之和');
    
    
    subplot(348);
    plot(kf);hold on;
    
    %     set(gca,'XTick',1:1:frame);
    xlim([1 frame]);
    ylim([0.4 1]);
    hold on; plot([0 frame],[0.86 0.86],'r--');
    hold on; plot([0 frame],[0.7 0.7],'g--');
    title('kf-原始最终输出');
    
    subplot(3,4,9);
    plot(first_Intergral);
    xlim([10,frame]);        % 必须是5 开始，否则前几帧数值太大
    title('kf一次导');
    
    subplot(3,4,10);
    plot(second_Integral);
    xlim([10,frame]);        % 必须是5 开始，否则前几帧数值太大
    title('二次导');
    
    subplot(3,4,11);
    plot(Spike2);
    xlim([5,frame]);        % 必须是5 开始，否则前几帧数值太大
    title('Spike2');
    
    subplot(3,4,12);
    plot(df,'m');hold on;
    plot(kf,'c');
    %             axis([63 frame 0.4 1]);
    %         set(gca,'XTick',65:1:frame);
%     xlim([4,frame]);
    ylim([0.4 1]);
    hold on; plot([0 frame],[0.86 0.86],'r--');
    hold on; plot([0 frame],[0.7 0.7],'g--');
    title(ntimes);
    
    
    suptitle(num2str(k));
    hold on;
    
    
end




toc;    % 计时结束



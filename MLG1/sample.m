% LGMD's Matlab implementation
% Xuelong Sun  01/09/2016  CIL in UoL

clc;clear all;

% Parameters setting in the LGMD model
  % Size of Image
Ver = 240;  %图片高度240
Hor = 432;  %水平宽度432
  % The total frame to test
frame = 60;
  % The persistence of the luminance change 
np = 1;
u = 1;
  % local inhibition weight
wi = [0.125,0.25,0.125;0.25,0,0.25;0.125,0.25,0.125];
  % global inhibition weight
WI = 0.3;
  % convolution mask (from S-->G)
we = ones(3,3).*1/9;
  % Scale (from S-->G)
Cw = 4;
delta_c = 0.01;
  % decay coefficient range (0,1)
Cde = 0.5;
  % dacay threshold
Tde = 15;
  % Flag of using fixed Ts, 1 for use and 0 for no 
Flag_FixedTs = 0;
  % Flag of using fixed Tffi, 1 for use and 0 for no 
Flag_FixedTffi = 1;

if Flag_FixedTffi == 1
  % fixed FFI threshold
  Tffi = 80;
else
  % Parameters in adapatable FFI threshold
  Tffi = 0;
  AFffi = 0.02;
  Tfo = 10;
end

if Flag_FixedTs == 1
  % fixed spiking threshold (No FFM cell)
  Ts = 0.8;
else
  % Parameters in FFM cell - No evidence in locust
  Ts = 0;
  AFlt = 4;
  AFmp = 1;
  Tmp = 0.86;
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
for i=1:10
    p(i) = 1/(1+exp(u*i));
end
  % The E/I layer
Ef = uint8(zeros(Ver,Hor));
If = zeros(Ver,Hor);
  % The S layer,note S may be negative
Sf = zeros(Ver,Hor);
  % The G layer - No evidence in locust
Gf = zeros(Ver,Hor);
Cef = zeros(Ver,Hor);
w = 0;
  % LGMD cell,note that I set it as a 1*20 vector to stored it for abserving 
Kf = zeros(1,frame);
  % The FFI cell
FFI = zeros(1,frame);
  % The spiking
Spike = zeros(1,frame);

   
for k = 1:frame
    % Make some of the data used in the last step be zero
    If = zeros(Ver,Hor);
    Cef = zeros(Ver,Hor);
    % tansmit the data to the Lf-1,Pf-1 from the Lf, Pf-1
    Lf_1 = Lf;
    Pf_1 = Pf;

     
     
    % load the current Image
    if k >= 100
        FileName = strcat(('D:\Code\Matlab\LGMD-xuelong\ImageTest\00'),num2str(k),'.jpg');
    else if k>=10
            FileName = strcat(('D:\Code\Matlab\LGMD-xuelong\ImageTest\00'),num2str(k),'.jpg');
        else
            FileName = strcat(('D:\Code\Matlab\LGMD-xuelong\ImageTest\00'),num2str(k),'.jpg');
        end
    end 
    Lf = imread(FileName,'jpg');
    
    % my code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [m,n]=size(Lf);
% 
% newGrayPic=Lf;
% 
% sobelNum=0;
% 
% sobelThreshold=0.7;
% 
% for j=2:m-1
% 
%     for k=2:n-1
% 
%         sobelNum=abs(Lf(j-1,k+1)+2*Lf(j,k+1)+Lf(j+1,k+1)-Lf(j-1,k-1)-2*Lf(j,k-1)-Lf(j+1,k-1))+abs(Lf(j-1,k-1)+2*Lf(j-1,k)+Lf(j-1,k+1)-Lf(j+1,k-1)-2*Lf(j+1,k)-Lf(j+1,k+1));
% 
%         if(sobelNum > sobelThreshold)
% 
%             newGrayPic(j,k)=255;
% 
%         else
% 
%             newGrayPic(j,k)=0;
% 
%         end
% 
%     end
% 
% end
%     Lf=newGrayPic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LGMD algorithm
    Pf = abs(Lf - Lf_1) + sum(p(1:np));
    Ef = Pf;
    for i = 2:Ver-1
        for j = 2:Hor-1
            for s = -1:1
                for t = -1:1
                    if s ~= 0 || t ~= 0
                        If(i,j) = If(i,j) + Pf_1(i+s,j+t).*wi(s+2,t+2);
                    end
                end
            end
        end
    end
    % Note that I ignored the boundary,
    % i.e.,If(1,:)=If(:,1)=If(:,Hro)=If(Ver,:)=0
    Sf = double(Ef) - If.*WI;
    for i = 2:Ver-1
        for j = 2:Hor-1
            for s = -1:1
                for t = -1:1
                    Cef(i,j) = Cef(i,j) + Sf(i+s,j+t).*wi(s+2,t+2);  
                end
            end
        end
    end
    w = delta_c + max(max(abs(Cef)))/Cw;
    Gf = Sf.*Cef/w;
    for i = 1:Ver
       for j = 1:Hor
           if Gf(i,j)*Cde < Tde
               Gf(i,j) = 0;
           end
       end
    end
    Kf(1,k) = sum(sum(abs(Gf)));
    Kf(1,k) = 1/(1+exp(-Kf(1,k)/(Ver*Hor)));
    
    if ~Flag_FixedTffi
        % adaptable Tffi
        Tffi = Tfo + AFffi*Tffi;
    end
          
    if sum(sum(Pf_1))/(Ver*Hor) > Tffi
        FFI(1,k) = 1;
    end
    
    if ~Flag_FixedTs
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
        Ts = AFlt*Tlt + AFmp*Tmp;
    end
%   movie
    if Kf(1,k) >= Ts && FFI(1,k)==0
        Spike(1,k) = 1;
    end
    


%         subplot(2,2,1);imshow(uint8(Lf));title('Output of L layer');
%         subplot(2,2,2);imshow(uint8(Pf));title('Output of P layer');
%         subplot(2,2,3);imshow(uint8(Sf));title('Output of S layer');
%         subplot(2,2,4);imshow(uint8(Gf));title('Output of G layer');
   

end


% figure(k+1);
% subplot(1,2,1),plot(Kf),title('Output of LGMD');
% subplot(1,2,2),plot(Spike,'r'),title('Spike');
% 
figure(k+1);
title('test_result');
plot(Kf);
hold on;
plot(Spike,'r');
plot(FFI,'k');
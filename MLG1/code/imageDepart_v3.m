% % v2: overlap，互相有一部分重叠
% % v3: 重新规划大小

function departedImagecombination=imageDepart_v3(imagesMatrix,partNum)
% 分割图像序列

overlapDegree = 7.5;    % 互相遮盖的部分

imagesNums=size(imagesMatrix,3);        %读取矩阵 Z尺寸
l=size(imagesMatrix,1);       %读取图像的尺寸

imagesMatrix2 = zeros((l/2)+45,(l/2)+45,imagesNums);   % 声明返回矩阵

cox=meshgrid(-(l-1)/2:(l-1)/2);       %生成l+1个点位置的横坐标,矩阵大小为(l+1)*(l+1)
coy=-cox.';                   %生成每个点位置的纵坐标
cod=rad2deg(atan2(coy,cox));  %生成d矩阵,矩阵大小为(l+1)*(l+1)；表明矩阵每一点和中心点的角度
pic=zeros(l);     %对于每帧图像都要重置0矩阵

for i=1:imagesNums
    for j=1:l
        for k=1:l
            if cod(j,k)>((partNum-1) * 22.5 - overlapDegree) && cod(j,k)<(partNum * 22.5 + overlapDegree)      %这里的22.5就是360/16=22.5，因此是将图像分割为16块
                pic(j,k)=1;
            end
        end
    end
    %     figure();imshow(mat2gray(pic));
    %     figure();imshow(mat2gray(imagesMatrix(:,:,i)));
    imagesMatrix(:,:,i)=imagesMatrix(:,:,i).*pic;
    %     figure();imshow(mat2gray(imagesMatrix(:,:,i)));
    
    
    if (1 <= partNum) && (partNum <= 4)
        imagesMatrix2(:,:,i) = imagesMatrix(1:405,316:720,i);
    elseif (5 <= partNum) && (partNum <= 8)
        imagesMatrix2(:,:,i) = imagesMatrix(1:405,1:405,i);
    elseif (-3 <= partNum) && (partNum <= 0)
        imagesMatrix2(:,:,i) = imagesMatrix(316:720,316:720,i);
    elseif (-7 <= partNum) && (partNum <= -4)
        imagesMatrix2(:,:,i) = imagesMatrix(316:720,1:405,i);        
    end
    
end





departedImagecombination=imagesMatrix2;
end
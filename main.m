clear
clc
close all;
delt_heta=15;%偏振角采集步长; Polarization angle acquisition step

Image_num=1;

cd 'Scene1_LT\';
k=1;
Iox=double(imread('clear.bmp'))/255;%读取清晰水下图像; Read clear underwater images
for i=0:delt_heta:165
    filename=[num2str(i),'_c.bmp'];
    Image1=(double(imread(filename))/255);
    Image(k,:,:)=Image1;
    k=k+1;
end 

[m,x,y]=size(Image);

for i=1:1:m/2
    Image_pair(i,1,:,:)=Image(i,:,:);
    Image_pair(i,2,:,:)=Image(i+6,:,:);
end

%相关参数的计算
Image_sum1=zeros(x,y);
Image_sum2=zeros(x,y);

for i=1:1:m
    image1(:,:)=Image(i,:,:);
    Image_sum1(:,:)=Image_sum1(:,:)+image1(:,:)*(delt_heta/180*pi);
end
for i=1:1:m/2
    image1(:,:)=Image_pair(i,1,:,:);
    image2(:,:)=Image_pair(i,2,:,:);
    Image_sum2(:,:)=Image_sum2(:,:)+(image1-image2).*(image1-image2)*(delt_heta/180*pi);
end
I=2*Image_sum1/pi;
P=sqrt(4*Image_sum2/pi)./I;% 计算偏振度; Calculate the degree of polarization

I_max1=I/2+(P.*I)/2;
I_min1=I/2-(P.*I)/2;

Pmean=mean(mean(P));

[m,n]=size(I);

Ip=I.*(P);%偏振光图像; Polarized light image
Iq=(I.*(1-P));%非偏振光图像; Polarized light image

Ix=I;
[m,n]=size(Ix);

c1= 50;%c1的取值一般介于20~120之间，50 for Scene1_LT,80 for Scene1_HT, 55 for Scene2_LT,100 for Scene2_HT

K1=1/(sqrt(2*pi)*c1);

for i=1:m
    for j=1:n
        g1(i,j)=K1*exp((-((i)^2+(j)^2))/(c1*c1));%高斯函数; Gaussian function

    end
end
g1fft=(fft2(g1));%对高斯函数进行二维傅里叶变换; Two dimensional Fourier transform of Gaussian function
Ifft=(fft2(Ix));%对原始水下图像进行二维傅里叶变换; 2D Fourier transform of original underwater image
L1fft=g1fft.*Ifft;%卷积; convolution

L1=ifft2(L1fft);%二维傅里叶逆变换; Two-dimensional inverse Fourier transform
R=log(Ix)-log(L1);%在对数域中，原图像减去低通滤波后的图像，得到高频增强的图像R;
% In the log domain, the original image subtracts the low-pass filtered image to obtain the high-frequency enhanced image R
r=exp(R);
l=Ix./r;
r_mean=r*mean(mean(l));%将入射光分量进行平均; Average the incident light components
Ip_r=r_mean.*P;
Iq_r=r_mean.*(1-P);%分别计算出均匀化照明后的偏振光和非偏振光图像;
% Calculate polarized light and unpolarized light images after uniform illumination

I_fusion=Ip_r.*(P)/Pmean+Iq_r.*Pmean./(P+0.2);
%delta 0.2 for Scene1_LT,0.25 for Scene1_HT, 0.08 for Scene2_LT,0.12 for Scene2_HT

% for i=1:m
%     for j=1:n
%         if(I_fusion(i,j)>0.8)
%             I_fusion(i,j)=0.8;
%         end
%     end
% end
% 对于场景2，部分点存在过曝的情况，因此，通过阈值消除过曝点对于归一化的影响
% For Scene 2, some points are overexposed, so the influence of overexposure points on normalization is eliminated through threshold

I_fusion=(I_fusion-min(I_fusion(:)))./(max(I_fusion(:))-min(I_fusion(:)));
%融合图像归一化
figure(1)
subplot (1,2,1)
imshow(I)
title('原始水下图像')
subplot (1,2,2)
imshow(I_fusion)
title('复原水下图像')
peme=eme(I_fusion,3);
pesnr=psnr(I_fusion,Iox);
function zxg = zixiangguan(input)

I1 = input;
% I1 = I1(:,:,1);
%figure(1) % 显示输入图像1
% colormap(gray(255));
%image(I1)
axis off
FI1 = fft2(I1);
max1 = max(FI1);
max2 = max(max1);
scale = 1.0/max2;
FI1 = FI1.*scale;
I2 = I1; % 输入图像2，输入相同的图像
% I2 = I2(:,:,1);
% 
% I2 = I2(225:257,225:257);

% figure(2) % 显示输入图像2
% colormap(gray(255));
% image(I2)
% axis off
FI2 = fft2(I2);
max1 = max(FI2);
max2 = max(max1);
scale = 1.0/max2;
FI2 = FI2.*scale;
% 进行相关性的计算
FPR = FI1.*conj(FI2); 
PR = ifft2(FPR);
PR = fftshift(PR);
max1 = max(PR);
max2 = max(max1);
scale = 1.0/max2;
PR = PR.*scale;
zxg=PR;
% figure(3) % 空间域相关显示
% colormap(gray(255));
% image(abs(256 * PR)); % 如果没有乘上256，基本都为黑
% axis off



end
function zxg = zixiangguan(input)

I1 = input;
% I1 = I1(:,:,1);
%figure(1) % ��ʾ����ͼ��1
% colormap(gray(255));
%image(I1)
axis off
FI1 = fft2(I1);
max1 = max(FI1);
max2 = max(max1);
scale = 1.0/max2;
FI1 = FI1.*scale;
I2 = I1; % ����ͼ��2��������ͬ��ͼ��
% I2 = I2(:,:,1);
% 
% I2 = I2(225:257,225:257);

% figure(2) % ��ʾ����ͼ��2
% colormap(gray(255));
% image(I2)
% axis off
FI2 = fft2(I2);
max1 = max(FI2);
max2 = max(max1);
scale = 1.0/max2;
FI2 = FI2.*scale;
% ��������Եļ���
FPR = FI1.*conj(FI2); 
PR = ifft2(FPR);
PR = fftshift(PR);
max1 = max(PR);
max2 = max(max1);
scale = 1.0/max2;
PR = PR.*scale;
zxg=PR;
% figure(3) % �ռ��������ʾ
% colormap(gray(255));
% image(abs(256 * PR)); % ���û�г���256��������Ϊ��
% axis off



end
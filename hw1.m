clearvars -except Undata
close all; clc;
if exist('Undata', 'var')==0
    load Testdata
end

L=15; % spatial domain
n=64; % Fourier modes
T = 20;
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

Un = zeros([T n n n]);

% Creation of the 4D array representing the 3D ultrasonic image
% at the 20 time slices (dimension 1)
for j=1:T
    Un(j,:,:,:)=reshape(Undata(j,:),n,n,n);
end

% Averaging process; we average along the time axis
Unt = fftshift(fftn(Un));

UntavgABS = squeeze(mean(abs(Unt),1));

[cfABS, cfI] = max(UntavgABS(:));

figure(1)
isosurface(Kx,Ky,Kz,UntavgABS/cfABS,.725)
   axis([-20 20 -20 20 -20 20]./2), grid on, %drawnow
   
xlabel('K_x')
ylabel('K_y')
zlabel('K_z')

[cfIKx, cfIKy, cfIKz] = ind2sub([n n n], cfI);

% 3D coordinates of the central frequency (presumably the marble)
cfkx = Kx(cfIKx, cfIKy, cfIKz);
cfky = Ky(cfIKx, cfIKy, cfIKz);
cfkz = Kz(cfIKx, cfIKy, cfIKz);

title(sprintf('Central Frequency Spike at (%.3f,%.3f,%.3f)', cfkx,cfky,cfkz))

% hold on
% m = scatter3(cfkx,cfky,cfkz,25,[0 0 0]);
% m.MarkerFaceColor = [1 0 0];

% 3D Gaussian filter centered at the marble's frequency in k (frq) space
s = 2^0.5; % range of 2^-(0.5, 0.5) give good results (same start, end)
Gfilter = (2*pi)^(-3/2)*s^-3 * ...
    exp(-(1/(2*s^3))*((Kx - cfkx).^2 + (Ky - cfky).^2 + (Kz - cfkz).^2));

% Replication of the filter for multiplication
G = zeros(size(Unt));
for j = 1:T
    G(j,:,:,:) = Gfilter;
end

% Filtered data
UnF = ifftn(ifftshift(G.*Unt));

figure(2)

for j = 1:T
    Ua = abs(squeeze(UnF(j,:,:,:)));
    Uam = max(Ua(:));
    isosurface(X,Y,Z,Ua/Uam,.75)
    hold on
end

axis([-20 20 -20 20 -20 20]), grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

title('Reconstructed Marble At 20 Time Slices')

% Extract the path from the data
mx = zeros([20 1]); my = mx; mz = mx;

for j = 1:T
    Unf = squeeze(UnF(j,:,:,:));
    [~, I] = max(Unf(:));
    [my(j), mx(j), mz(j)] = ind2sub(size(Unf), I);
end

fig = figure(3);
sp1 = subplot(2,3,[1 2 3]);
p = plot3(x(mx), y(my), z(mz), 'Color', [0 0 0], 'LineWidth', 1);

hold on

markersize = 75;
arrowstyle = 'plain';

h = scatter3(x(mx(1)), y(my(1)), z(mz(1)), markersize, [0 0 0], 'v');
h.MarkerFaceColor = [1 0.5 0];
g = scatter3(x(mx(end)), y(my(end)), z(mz(end)), markersize, [0 0 0], 'o');
g.MarkerFaceColor = [0 0.75 0.75];

legend('Marble''s Path', 'Start Point', 'End Point', 'Location', 'NorthWestOutside')%, 'FontSize', 16)

xlabel('x (mm)')%, 'FontSize', 16)
ylabel('y (mm)')%, 'FontSize', 16)
zlabel('z (mm)')%, 'FontSize', 16)

title('Marble''s 3D Path Over the 20 Time Slices')%, 'FontSize', 32)

axis([-20 20 -20 20 -20 20]), grid on

sp3 = subplot(2,3,4);
u = zeros(size(x)); v=u; w=u;

xvecs = diff(x(mx)); yvecs = diff(y(my)); zvecs = diff(z(mz));
xvecs = [xvecs 0]; yvecs = [yvecs 0]; zvecs = [zvecs 0];

u(mx) = xvecs; v(my) = yvecs; w(mz) = zvecs;

q = quiver(x(mx),y(my),xvecs,yvecs,0,'Color',[0 0 0],'ShowArrowHead','off');

qU = q.UData;
qV = q.VData;
qX = q.XData;
qY = q.YData;


headWidth = 2;
headLength = 3;
LineLength = 1; %0.08;

for ii = 1:length(qX)-1
    
    ah = annotation('arrow',...
        'headStyle',arrowstyle,'HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[qX(ii) qY(ii) LineLength*qU(ii) LineLength*qV(ii)]);
    
end

hold on


k = scatter(x(mx(1)), y(my(1)), markersize, [0 0 0], 'v');
k.MarkerFaceColor = [1 0.5 0];
l = scatter(x(mx(end)), y(my(end)), markersize, [0 0 0], 'o');
l.MarkerFaceColor = [0 0.75 0.75];

%legend('Marble''s Path', 'Start Point', 'End Point', 'Location', 'NorthWest')

xlabel('x (mm)')%, 'FontSize', 16)
ylabel('y (mm)')%, 'FontSize', 16)
title('Marble''s 2D Path (x,y)')%, 'FontSize', 24)
axis([-15 15 -15 15]), grid on

sp4 = subplot(2,3,5);

q2 = quiver(x(mx),z(mz),xvecs,zvecs,0,'Color',[0 0 0],'ShowArrowHead','off');

qU = q2.UData;
qW = q2.VData;
qX = q2.XData;
qZ = q2.YData;


for ii = 1:length(qX)-1
    
    ah = annotation('arrow',...
        'headStyle',arrowstyle,'HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[qX(ii) qZ(ii) LineLength*qU(ii) LineLength*qW(ii)]);
    
end

hold on


k = scatter(x(mx(1)), z(mz(1)), markersize, [0 0 0], 'v');
k.MarkerFaceColor = [1 0.5 0];
l = scatter(x(mx(end)), z(mz(end)), markersize, [0 0 0], 'o');
l.MarkerFaceColor = [0 0.75 0.75];

%legend('Marble''s Path', 'Start Point', 'End Point', 'Location', 'NorthWest')

xlabel('x (mm)')%, 'FontSize', 16)
ylabel('z (mm)')%, 'FontSize', 16)
title('Marble''s 2D Path (x,z)')%, 'FontSize', 24)
axis([-15 15 -15 15]), grid on



sp21 = subplot(2,3,6);

q3 = quiver(y(my),z(mz),yvecs,zvecs,0,'Color',[0 0 0],'ShowArrowHead','off');

qV = q3.UData;
qW = q3.VData;
qY = q3.XData;
qZ = q3.YData;


for ii = 1:length(qY)-1
    
    ah = annotation('arrow',...
        'headStyle',arrowstyle,'HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[qY(ii) qZ(ii) LineLength*qV(ii) LineLength*qW(ii)]);
    
end

hold on


k = scatter(y(my(1)), z(mz(1)), markersize, [0 0 0], 'v');
k.MarkerFaceColor = [1 0.5 0];
l = scatter(y(my(end)), z(mz(end)), markersize, [0 0 0], 'o');
l.MarkerFaceColor = [0 0.75 0.75];

%legend('Marble''s Path', 'Start Point', 'End Point', 'Location', 'NorthWest')

xlabel('y (mm)')%, 'FontSize', 16)
ylabel('z (mm)')%, 'FontSize', 16)
title('Marble''s 2D Path (y,z)')%, 'FontSize', 24)
axis([-15 15 -15 15]), grid on

xstart = x(mx(1)), ystart = y(my(1)), zstart = z(mz(1))
xend = x(mx(end)), yend = y(my(end)), zend = z(mz(end))

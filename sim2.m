% 
% XwellMultiShota.m
% J.M. Harris, May 2015
%
% This version can simulate a survey with many shots and receivers.
% Model and survey dimensions are in meters; velocities in m/sec;
%
%
% Computational size (sampels) is calculated from input model (meters) + a 40 pnt sidepad : -> 
%
%
clear all;
hold off;
tstart = tic;

% background velocity set to minmum values for entire model
cmin       = 2000;

% setting up the runtime parameters
freq   = 300;
lam    = cmin/(2*freq);
dx     = lam/4;
dz     = lam/4;

%Build Velocity Model (m)
xmin = 0;
xmax = 400;
zmin = 0;
zmax = 400;

% model dimensions (grid)
idim    = round((xmax-xmin)/dx) + 1;
jdim    = round((zmax-zmin)/dz) + 1;

% building blank models
velmod  = ones(jdim,idim).*cmin;
velmod2 = velmod;
lmod    = zeros(jdim,idim);

% Add perturbations to background velocity & attenuation to define model
%
% z z x x

ztop = 261;
zbot = 266;
xleft = 170;
xright = 175;
jz1    = round((ztop-zmin)/dz) + 1;
jz2    = round((zbot-zmin)/dz) + 1;
jx1    = round((xleft-xmin)/dx) + 1;
jx2    = round((xright-xmin)/dx) + 1;
velmod2(jz1:jz2,jx1:jx2) = 3200;


ztop = 231;
zbot = 236;
xleft = 220;
xright = 225;
jz1    = round((ztop-zmin)/dz) + 1;
jz2    = round((zbot-zmin)/dz) + 1;
jx1    = round((xleft-xmin)/dx) + 1;
jx2    = round((xright-xmin)/dx) + 1;
velmod2(jz1:jz2,jx1:jx2) = 3200;

% %add 50m thick layer with velocity of 3000 at 50m
% ztop = 50;
% zbot = 110;
% j1    = round((ztop-zmin)/dz) + 1;
% j2    = round((zbot-zmin)/dz) + 1;
% velmod2(j1:j2,:) = 2750;
% %lmod(j1:j2,:) = 0.00002;
% 
% %add dipping layer with velocity of 3500 between 90m and 110m
% ztop = 90;
% zbot = 110;
% j1    = round((ztop-zmin)/dz) + 1;
% j2    = round((zbot-zmin)/dz) + 1;
% i2    = idim;
% for j=j1:j2;
%     i1 = idim - round((j-j1)*(idim-1)/(j2-j1));
%     for i = i1:i2;
%     velmod2(j,i) = 3750;
% %    lmod(j1:j2,:) = 0.00001;
%     end
% end
% %add 20m thick layer with velocity of 3500 at 130m
% ztop = 110;
% zbot = 130;
% j1    = round((ztop-zmin)/dz) + 1;
% j2    = round((zbot-zmin)/dz) + 1;
% velmod2(j1:j2,:) = 3750;
% %lmod(j1:j2,:) = 0.00001;
% 
% %add 20m thick layer with velocity of 3500 at 130m
% ztop = 130;
% zbot = 150;
% j1    = round((ztop-zmin)/dz) + 1;
% j2    = round((zbot-zmin)/dz) + 1;
% velmod2(j1:j2,:) = 3000;
% %lmod(j1:j2,:) = 0.00002;
% 
% %add 50m thick layer with velocity of 4000 at 150m
% ztop = 150;
% zbot = 200;
% j1    = round((ztop-zmin)/dz) + 1;
% j2    = round((zbot-zmin)/dz) + 1;
% velmod2(j1:j2,:) = 4000;
% %lmod(j1:j2,:) = 0.00001;

% pixel center locations
xl = (0:idim-1).*dx;
zl = (0:jdim-1).*dz;

timCnt    = 1000;

% build receiver array
rCnt = 3;
xr = linspace(150,250,rCnt);
zr = linspace(200,200,rCnt);

% build source locations
sCnt = 5;
xs = linspace(200,200,sCnt);
zs = linspace(290,310,sCnt);



% pixel center locations
xl = (0:idim-1).*dx;
zl = (0:jdim-1).*dz;

% set display parameters
quiet = 0;
visual = 1;
snapR = 50;

% plotting the velocity model and the resulting synthetic seismogram
subplot(1,2,1);
imagesc(xl,zl,velmod2); hold on;
plot(xr,zr,'r.'); 
q1 = colorbar; q2 = get(q1,'YLabel');
set(q2,'String','V_p (m/s)');
%axis equal;
axis tight;
xlabel('X (m)'); ylabel('Depth (m)'); title('Velocity Model');


% do the modeling
[seis,tarr,s,t,snaps] = VAMultiShot(dx,dz,velmod2,lmod,xr,zr,freq,...
                                    sCnt,xs,zs,timCnt,0,quiet,visual,snapR);
tarr = 1000*tarr;

% plotting the velocity model and the resulting synthetic seismogram
subplot(1,2,1);
imagesc(xl,zl,velmod2); hold on;
plot(xr,zr,'r.'); 
q1 = colorbar; q2 = get(q1,'YLabel');
set(q2,'String','V_p (m/s)');
%axis equal;
axis tight;
xlabel('X (m)'); ylabel('Depth (m)'); title('Velocity Model');

di = 1;
for is =1:di:sCnt;
xshot = xs(is);
zshot = zs(is);
disp(cat(1,['Shot Depth = ',num2str(zshot)]));
subplot(1,2,1);
plot(xshot,zshot,'ro'); 

dseis(:,:) = seis(is,:,:);

subplot(1,2,2);
imagesc(tarr,zr,dseis); hold on;
axis([min(tarr),max(tarr),min(zl),max(zl)]);
caxis([-0.025,0.025]);
q1 = colorbar; q2 = get(q1,'YLabel');
set(q2,'String','Amplitude');
set(gca,'xgrid','on')

ylabel('Receiver Z (m)'); xlabel('Time (mS)');
title(cat(1,['Synthetic: CS = ',num2str(zshot)]));

orient landscape;

end
%save('modela','velmod2','xl','zl','seis','sCnt','xs','zs','rCnt','xr','zr','timCnt','tarr');

totalTime= toc(tstart);

disp(cat(1,['Total time =',num2str(totalTime)]));

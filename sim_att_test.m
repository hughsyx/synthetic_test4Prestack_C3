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
xmax = 1000;
zmin = 0;
zmax = 300;

% model dimensions (grid)
idim    = round((xmax-xmin)/dx) + 1;
jdim    = round((zmax-zmin)/dz) + 1;

% building blank models
velmod  = ones(jdim,idim).*cmin;
velmod2 = velmod;
lmod    = ones(jdim,idim)*0.00002;

% Add perturbations to background velocity & attenuation to define model
%
% z z x x

% ztop = 261;
% zbot = 266;
% xleft = 200;
% xright = 205;
% jz1    = round((ztop-zmin)/dz) + 1;
% jz2    = round((zbot-zmin)/dz) + 1;
% jx1    = round((xleft-xmin)/dx) + 1;
% jx2    = round((xright-xmin)/dx) + 1;
% velmod2(jz1:jz2,jx1:jx2) = 3200;

% for i=1:10
%     ztop = round(rand(1)*40+230);
%     zbot = ztop+5;
%     xleft = round((rand(1)*2-1)*50+500-2);
%     xright = xleft+5;
%     jz1    = round((ztop-zmin)/dz) + 1;
%     jz2    = round((zbot-zmin)/dz) + 1;
%     jx1    = round((xleft-xmin)/dx) + 1;
%     jx2    = round((xright-xmin)/dx) + 1;
%     velmod2(jz1:jz2,jx1:jx2) = 2400;
% end


% pixel center locations
xl = (0:idim-1).*dx;
zl = (0:jdim-1).*dz;

timCnt    = 1000;

% build receiver array
rCnt = 91;
xr = linspace(500,500,rCnt);
zr = linspace(180,260,rCnt);

% build source locations
% sCnt = 100;
% xs = round((rand(1,sCnt)*2-1)*40+200);
% zs = round(rand(1,sCnt)*40);

sCnt = 1;
xs = linspace(500,500,sCnt);
zs = linspace(100,160,sCnt);



% pixel center locations
xl = (0:idim-1).*dx;
zl = (0:jdim-1).*dz;

% set display parameters
quiet = 0;
visual = 1;
snapR = 100;

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
%[seis2,tarr2,s2,t2,snaps2] = VAMultiShot(dx,dz,velmod2,lmod,xs,zs,freq,...
%                                    rCnt,xr,zr,timCnt,0,quiet,visual,snapR);
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

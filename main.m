% Core procedure for correction of PACT images in non-uniform medium:
% Isolating and revealing aberrations after medium is measured with 
% the reflection matrix Rxx
% Copyright (C) 2026 Tianhua Zhou
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% basic params
Nelem = 128; % receive element count
TxNum = 65;
TxAngles = -16:0.5:16;
elemwidth = 0.3;
receivers_xpos = 0:elemwidth*1e-3:(Nelem-1)*elemwidth*1e-3;
receivers_xpos = receivers_xpos - mean(receivers_xpos);
receivers_zpos = zeros(1,Nelem);
receivers_pos = [receivers_xpos.' receivers_zpos.'];
elementNormal = repmat([0,1],Nelem,1);
fs = 40;
fc = 7.5;

%% reconstruction grid setup
dx = 1e-4; % [m]
dz = 1e-4; % [m]
xrange = max(receivers_xpos) - min(receivers_xpos); % [m]
zrange = 60e-3; % [m]
Nx = round(xrange / dx) + 1;
Nz = round(zrange / dz) + 1;
xvec = -xrange/2:dx:xrange/2;
zvec = 0:dz:zrange;
[x,z] = meshgrid(xvec,zvec);
sos = 1480; % [m/s]
sos_PW = 1540; % [m/s] sos used to convert steering angles to transmit delays

rMParams.xin = 1:length(xvec);
rMParams.xout = rMParams.xin;
rMParams.xinm = xvec(rMParams.xin);
rMParams.xoutm = xvec(rMParams.xout);
rMParams.Nin = length(rMParams.xin);
rMParams.Nout = length(rMParams.xout);

zstart = 50;
zend = 500;
zrem = zend - zstart + 1;

%% calculation of projection matrix T0
kc = 2*pi*(fc*1e6)/sos;
kvector = -11e3:50:11e3;
thetavector = asind(kvector ./ kc);
T0 = exp(1i*kvector.'.*rMParams.xoutm); % T0: [kout, xout]
figure,imagesc(abs(T0*T0')); % check completeness of projection
ax = gca;
ax.DataAspectRatio = [1,1,1];

%% calculate valid range of k for each x based on receiving angle and z
leftmostxm = rMParams.xoutm(1);
rightmostxm = rMParams.xoutm(end);
leftHAngle = atand((rightmostxm - x) ./ z);
rightHAngle = - atand((x - leftmostxm) ./ z);
leftHAngle(leftHAngle > 44) = 44;
rightHAngle(rightHAngle < -44) = -44;
leftHAngle(1:zstart-1,:) = [];
rightHAngle(1:zstart-1,:) = [];
leftThetaInd = zeros(size(leftHAngle));
rightThetaInd = zeros(size(rightHAngle));
for i = 1:zrem
    for j = 1:length(rMParams.xin)
        indLeft = find(thetavector > leftHAngle(i,j),1);
        if isempty(indLeft)
            indLeft = length(thetavector);
        end
        indRight = find(thetavector < rightHAngle(i,j),1,'last');
        if isempty(indRight)
            indRight = 1;
        end
        leftThetaInd(i,j) = indLeft;
        rightThetaInd(i,j) = indRight;
    end
end
for i = 120:zrem
    leftThetaInd(i,:) = leftThetaInd(120,:);
    rightThetaInd(i,:) = rightThetaInd(120,:);
end

%% projection of reflection matrix into the far field
for i = 1:zrem
    Rkk(:,:,i) = T0 * Rxx(:,:,i) * T0.'; % Rxx is reflection matrix
end
figure,imagesc(abs(Rkk(:,:,250)));
ax = gca;
ax.DataAspectRatio = [1,1,1];
title('Raw Rkk @ a depth');
Rkkf = Rkk;
Rxxf = Rxx;

%% projection into dual Fourier-spatial basis
for i = 1:zrem
    Rkx_original(:,:,i) = T0 * Rxxf(:,:,i);
    Rxk_original(:,:,i) = Rxxf(:,:,i) * T0.';
end
Rkx_original = Rkx_original ./ abs(Rkx_original);
Rxk_original = Rxk_original ./ abs(Rxk_original);

%% apply bounds to Rkx
Rkx = Rkx_original;
Rxk = Rxk_original;
for i = 1:zrem
    for j = 1:size(Rxx,1)
        kMask = ones(size(T0,1),1);
        kMask(1:(rightThetaInd(i,j) - 1)) = 0;
        kMask((leftThetaInd(i,j) + 1):end) = 0;
        Rkx(~kMask,j,i) = 1;
%         Rkx(:,j,i) = kMask .* Rkx(:,j,i);
    end
end

%%
figure,imagesc([min(rMParams.xoutm),max(rMParams.xoutm)],[min(kvector),max(kvector)],angle(Rkx(:,:,250)));
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
title('Rkx');
xlabel('xin');
ylabel('kout');

%% distortion matrix
for i = 1:zrem
    Dkx(:,:,i) = Rkx_original(:,:,i) .* conj(T0);
    Dxk(:,:,i) = Rxk_original(:,:,i) .* T0';
end
for i = 1:zrem
    for j = 1:size(Rxx,1)
        kMask = ones(size(T0,1),1);
        kMask(1:(rightThetaInd(i,j) - 1)) = 0;
        kMask((leftThetaInd(i,j) + 1):end) = 0;
        Dkx(~kMask,j,i) = 0;
%         Rkx(:,j,i) = kMask .* Rkx(:,j,i);
        Dxk(j,~kMask,i) = 0;
    end
end

%% display phase of D matrix at single depth
figure,imagesc([min(rMParams.xoutm),max(rMParams.xoutm)],[min(kvector),max(kvector)],angle(Dkx(:,:,250)));
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
title('D');

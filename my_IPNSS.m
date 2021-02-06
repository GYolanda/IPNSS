derad = pi/180;                         % 角度->弧度
radeg = 180/pi;                         % 弧度->角度
twpi = 2*pi;
Melm = 16;                              % 天线数
kelm = 12;
L = 1024;                               % 快拍数
theta = [0 3.75];                       % 信号入射角度
dd = 0.15;                              % 振元间距
d = 0:dd:(Melm - 1)*dd;
f = [250e6, 250e6]';                    % 信号频率
snr1 = 10;                              % 信号1信噪比
snr2 = 20;                              % 信号2信噪比
fs = 2.6e6;                             % 采样频率
t = (0:1:L - 1)/fs;

A = exp(-1j*twpi*d.'*sin(theta*derad)); % 方向向量
% 构造接收信号
S = exp(1j*twpi*f*t);
S(1, :) = sqrt(10^(snr1/10)) * S(1, :);
S(2, :) = sqrt(10^(snr2/10)) * S(2, :);
x = A * S;                              % 不含噪声接收信号
noise = (randn(Melm, L) + 1j*randn(Melm, L))/sqrt(2);     % 噪声
X = x + noise;                          % 含噪声接收信号
Rxxm = X*X'/L;
issp = 2;                               % 设置平衡算法模式

% 空间平滑
if issp == 1
    Rxx = ssp(Rxxm, kelm);
elseif issp == 2
    Rxx = mssp(Rxxm, kelm);
else
    Rxx = Rxxm;
    kelm = Melm;
end

[EV, D] = eig(Rxx);                     % 特征值分解
EVA = diag(D)';
[EVA, I] = sort(EVA);                   % 特征值从小到大排序
EVA = fliplr(EVA);                      % 左右翻转，从大到小排序
EV = fliplr(EV(:, I));                  % 对特征向量排序

eta = trace(Rxx)/kelm;
% 构造IPNSS谱函数
angle = zeros(1441, 1);
G = zeros(1441, 1);
for iang = 1:1441
    angle(iang) = (iang - 721)/8;
    phim = derad*angle(iang);
    d2 = 0:dd:(kelm - 1)*dd;
    a = exp(-1j*twpi*d2*sin(phim)).';
    
    Rxx2 = Rxx + eta*(a*a');
    [EV2, D2] = eig(Rxx2);              % 特征值分解
    EVA2 = diag(D2)';
    [EVA2, I2] = sort(EVA2);            % 特征值从小到大排序
    EVA2 = fliplr(EVA2);                % 左右翻转，从大到小排序
    EV2 = fliplr(EV2(:, I2));           % 对特征向量排序
    
    tmp = 0;
    for i = (size(theta, 2)+1):kelm
        tmp = tmp + abs(EVA2(i) - EVA(i));
    end
    G(iang) = 1/tmp;
end

% 作图
Gmax = max(G);
G = 10*log10(G/Gmax);
g = plot(angle, G);
set(g, 'Linewidth', 1)
xlabel('angle (degree)');
ylabel('magtitude (dB)');
axis([-90 90 -50 0])
set(gca, 'XTick', (-90:30:90))
grid on
% -----------------------------------------------------------------%
function crs = ssp(cr, K)           % 前向平滑
M = size(cr, 1);
N = M - K + 1;
crs = zeros(K, K);
for in = 1:N
    crs = crs + cr(in:in+K-1, in:in+K-1);
end
crs = crs/N;
end

function crs = mssp(cr, K)          % 前后向平滑
M = size(cr, 1);
N = M - K + 1;
J = fliplr(eye(M));
crfb = (cr + J*cr.'*J)/2;
crs = zeros(K, K);
for in = 1:N
    crs = crs + crfb(in:in+K-1, in:in+K-1);
end
crs = crs/N;
end
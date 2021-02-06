derad = pi/180;                         % �Ƕ�->����
radeg = 180/pi;                         % ����->�Ƕ�
twpi = 2*pi;
Melm = 16;                              % ������
kelm = 12;
L = 1024;                               % ������
theta = [0 3.75];                       % �ź�����Ƕ�
dd = 0.15;                              % ��Ԫ���
d = 0:dd:(Melm - 1)*dd;
f = [250e6, 250e6]';                    % �ź�Ƶ��
snr1 = 10;                              % �ź�1�����
snr2 = 20;                              % �ź�2�����
fs = 2.6e6;                             % ����Ƶ��
t = (0:1:L - 1)/fs;

A = exp(-1j*twpi*d.'*sin(theta*derad)); % ��������
% ��������ź�
S = exp(1j*twpi*f*t);
S(1, :) = sqrt(10^(snr1/10)) * S(1, :);
S(2, :) = sqrt(10^(snr2/10)) * S(2, :);
x = A * S;                              % �������������ź�
noise = (randn(Melm, L) + 1j*randn(Melm, L))/sqrt(2);     % ����
X = x + noise;                          % �����������ź�
Rxxm = X*X'/L;
issp = 2;                               % ����ƽ���㷨ģʽ

% �ռ�ƽ��
if issp == 1
    Rxx = ssp(Rxxm, kelm);
elseif issp == 2
    Rxx = mssp(Rxxm, kelm);
else
    Rxx = Rxxm;
    kelm = Melm;
end

[EV, D] = eig(Rxx);                     % ����ֵ�ֽ�
EVA = diag(D)';
[EVA, I] = sort(EVA);                   % ����ֵ��С��������
EVA = fliplr(EVA);                      % ���ҷ�ת���Ӵ�С����
EV = fliplr(EV(:, I));                  % ��������������

eta = trace(Rxx)/kelm;
% ����IPNSS�׺���
angle = zeros(1441, 1);
G = zeros(1441, 1);
for iang = 1:1441
    angle(iang) = (iang - 721)/8;
    phim = derad*angle(iang);
    d2 = 0:dd:(kelm - 1)*dd;
    a = exp(-1j*twpi*d2*sin(phim)).';
    
    Rxx2 = Rxx + eta*(a*a');
    [EV2, D2] = eig(Rxx2);              % ����ֵ�ֽ�
    EVA2 = diag(D2)';
    [EVA2, I2] = sort(EVA2);            % ����ֵ��С��������
    EVA2 = fliplr(EVA2);                % ���ҷ�ת���Ӵ�С����
    EV2 = fliplr(EV2(:, I2));           % ��������������
    
    tmp = 0;
    for i = (size(theta, 2)+1):kelm
        tmp = tmp + abs(EVA2(i) - EVA(i));
    end
    G(iang) = 1/tmp;
end

% ��ͼ
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
function crs = ssp(cr, K)           % ǰ��ƽ��
M = size(cr, 1);
N = M - K + 1;
crs = zeros(K, K);
for in = 1:N
    crs = crs + cr(in:in+K-1, in:in+K-1);
end
crs = crs/N;
end

function crs = mssp(cr, K)          % ǰ����ƽ��
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
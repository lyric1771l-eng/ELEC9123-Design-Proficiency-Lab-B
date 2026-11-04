clear; clc; close all;
%% ======== 1. Define =======
% 1.1 Eb/N0 setting, conversion between linear and dB
EbN0_dB = 0:1:15; % Eb/N0 (dB) range
EbN0_lin = 10.^(EbN0_dB/10); % linear Eb/N0
numBits = 1e6; % Monte Carlo number of bits
fprintf('Eb/N0 range = [%d:%d] dB, number of bits = %g\n\n', EbN0_dB(1), EbN0_dB(end), numBits);

% 1.2 Pre-allocate for storage
Pb_awgn_theory = zeros(size(EbN0_lin));
Pb_lap_theory = zeros(size(EbN0_lin));
Pb_ray_theory = zeros(size(EbN0_lin));

Pb_awgn_sim = zeros(size(EbN0_lin));
Pb_lap_sim = zeros(size(EbN0_lin));
Pb_ray_sim = zeros(size(EbN0_lin));
%% ====== 2. Theoretical BER calculation ======
% - AWGN: Pb = Q(sqrt(2*Eb/N0))
% - Laplacian: Pb = 0.5 * exp(-sqrt(2*Eb/N0))

% 2.1 Define Function Q
Q = @(x) 0.5 * erfc(x ./ sqrt(2));

% 2.2 Traverse different Eb/N0 values and calculate the theoretical BER
for idx = 1:numel(EbN0_lin)
gamma_b = EbN0_lin(idx);

% 2.2.1 AWGN theoretical BER
Pb_awgn_theory(idx) = Q( sqrt(2*gamma_b) );

% 2.2.2 Laplacian theoretical BER
Pb_lap_theory(idx) = 0.5 * exp( - sqrt(2 * gamma_b) );
% 2.2.3Rayleigh (slow flat fading, no CSI, coherent detection)
Pb_ray_theory(idx) = 0.5 * (1 - sqrt(gamma_b./(1 + gamma_b)));
end

%% ====== 3. Monte Carlo Simulation ========
% 3.0 Generate random bits {0,1}, map to BPSK: 0 -> -1, 1 -> +1
% 3.1 AWGN Channel: Add Gaussian noise, make decisions and count the BER
% 3.2 Laplacian Channel: Generate Laplacian noise using the inverse transform method, make a decision and count the BER

% 3.0 Generate random bits and map to BPSK symbols
bits = randi([0,1], 1, numBits); % 0/1 bits
symbols = 2*bits - 1; % BPSK: 0->-1, 1->+1

% 3.1 AWGN and Laplacian noise simulation: Traverse Eb/N0
for idx = 1:numel(EbN0_lin)
gamma_b = EbN0_lin(idx);

% ========== 3.1 AWGN Simulation ==========
% sigma = sqrt( N0/2 ). For BPSK, sigma = 1/sqrt(2*Eb/N0)
sigma_awgn = sqrt(1/(2 * gamma_b));

% 3.1.1 Add Gaussian noise
noise_awgn = sigma_awgn * randn(1, numBits);

% 3.1.2 Add noise and make decisions
rx_awgn = symbols + noise_awgn;
bits_hat_awgn = (rx_awgn > 0); % Decision: >0 → 1, otherwise 0
err_awgn = bits_hat_awgn ~= bits; % Error bit vector
Pb_awgn_sim(idx) = mean(err_awgn); % BER

% ========== 3.2 Laplacian Simulation ==========
% Idea: Use inverse transform sampling method to generate Laplacian distribution:
% If U ~ Uniform(-0.5, 0.5), then Laplace(0, b) = -b * sign(U) * log(1 - 2|U|)
% Where b = N0/sqrt(2) = 1/(sqrt(2*Eb/N0))

b = sqrt(1/(2 * gamma_b)); % Laplacian parameter b
U = rand(1, numBits) - 0.5; % U ~ Uniform(-0.5,0.5)
noise_lap = - b * sign(U) .* log(1 - 2*abs(U));

% 3.2.1 Add noise and make decisions
rx_lap = symbols + noise_lap;
bits_hat_lap = (rx_lap > 0);
err_lap = bits_hat_lap ~= bits;
Pb_lap_sim(idx) = mean(err_lap);

% ==== 3.3 Rayleigh Fading Channel ====
% Rayleigh fading: |h|~Rayleigh(σ=1/sqrt(2)), E[|h|^2]=1
h = (randn(1, numBits) + 1j*randn(1, numBits))/sqrt(2); % Complex Rayleigh fading
rx_ray = h .* symbols; % Apply fading
sigma_ray = sigma_awgn; % Noise std same as AWGN (unit energy)
noise_ray = sigma_ray * (randn(1, numBits)+ 1j *randn(1, numBits));

% Received signal with fading and noise
rx_ray = rx_ray + noise_ray;
% Coherent receiver: divide by channel gain amplitude (assume perfect CSI)
rx_ray_eq = real(rx_ray ./ h); % Equalize fading (real part)
bits_hat_ray = (rx_ray_eq > 0);
err_ray = bits_hat_ray ~= bits;
Pb_ray_sim(idx) = mean(err_ray);
end

%% ==== 4. Laplacian Distribution Verification =======
% - 4.1 First select the b corresponding to any Eb/N0 (such as 3 dB), and generate a larger number of Laplacian noise samples
% - 4.2 Draw a sample histogram, superimpose the theoretical PDF; draw the empirical CDF, superimpose the theoretical CDF

% 4.1 Select Laplacian noise parameter b for Eb/N0 = 3 dB
chosen_idx = find(EbN0_dB == 3);
gamma_chosen = EbN0_lin(chosen_idx);
b = sqrt(1/(2 * gamma_chosen));
numSamples = 1e6; % Larger sample size for distribution verification

% 4.2 Generate Laplacian noise samples
U_val = rand(1, numSamples) - 0.5;
lap_samples = - b * sign(U_val) .* log(1 - 2*abs(U_val));

% 4.3 Define theoretical PDF and CDF functions
pdf_lap_theory = @(x) (1./(2*b)) .* exp(-abs(x)/b);
cdf_lap_theory = @(x) (x < 0) .* (0.5 * exp(x/b)) + (x >= 0) .* (1 - 0.5 * exp(-x/b));

% 4.4 Plot histogram + theoretical PDF
% same lap_samples
figure;
histogram(lap_samples, 'Normalization', 'pdf', ...
'BinWidth', 0.2, ...
'FaceColor', [0.2 0.6 0.9], ...
'EdgeColor', 'k'); 
hold on;
x_vals = linspace(min(lap_samples), max(lap_samples), 1000);
plot(x_vals, pdf_lap_theory(x_vals), 'r-', 'LineWidth', 1.5);
xlabel('x');
ylabel('Probability Density');
title('Laplacian Noise: Sample Histogram and Theoretical PDF (Eb/N0 = 3 dB)');
legend('Sample Histogram','Theoretical PDF');
grid on;

% 4.5 Plot empirical CDF + theoretical CDF
figure;% use date same as PDF 
[f_empirical, x_empirical] = ecdf(lap_samples); plot(x_empirical, f_empirical, 'b', 'LineWidth', 1.5); hold on;

x_vals = linspace(min(lap_samples), max(lap_samples), 1000);
plot(x_vals, cdf_lap_theory(x_vals), 'r--', 'LineWidth', 2);

xlabel('x'); ylabel('Cumulative Probability');
title('Laplacian Noise: Empirical CDF and Theoretical CDF (Eb/N0 = 3 dB)');
legend('Empirical CDF', 'Theoretical CDF');
grid on;

%% ======== 5. Plot BER Comparison Curves ========
figure;
semilogy(EbN0_dB, Pb_awgn_theory, 'b-', 'LineWidth',3); hold on;
semilogy(EbN0_dB, Pb_awgn_sim, 'bo--','LineWidth',1.2);
semilogy(EbN0_dB, Pb_lap_theory, 'r-', 'LineWidth',3);
semilogy(EbN0_dB, Pb_lap_sim, 'ro--','LineWidth',1.2);
semilogy(EbN0_dB, Pb_ray_theory, 'y-', 'LineWidth',3);
semilogy(EbN0_dB, Pb_ray_sim, 'yo--','LineWidth',1.2);

xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BPSK BER: AWGN vs Laplacian vs Rayleigh');
legend('AWGN Theory','AWGN Simulation','Laplacian Theory','Laplacian Simulation','Rayleigh Theory','Rayleigh Simulation','Location','SouthWest');
grid on;
ylim([1e-5, 1]);

%% ===== 6. Print Simulation vs Theory Table =========
%Output the BER values of each type and the simulation-theory deviation percentage at each Eb/N0

fprintf('\n%-8s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
'Eb/N0(dB)', 'AWGN Th', 'AWGN Sim', 'AW Dev(%)','Lap Th','Lap Sim','Lap Dev(%)');
fprintf('%s\n', repmat('-',1,80));

for idx = 1:numel(EbN0_dB)
% Deviation calculation
dev_aw = 100 * abs(Pb_awgn_sim(idx) - Pb_awgn_theory(idx)) / Pb_awgn_theory(idx);
dev_lp = 100 * abs(Pb_lap_sim(idx) - Pb_lap_theory(idx)) / Pb_lap_theory(idx);
fprintf('%-8d | %-12.4e | %-12.4e | %-12.2f | %-12.4e | %-12.4e | %-12.2f\n', ...
EbN0_dB(idx), ...
Pb_awgn_theory(idx), Pb_awgn_sim(idx), dev_aw, ...
Pb_lap_theory(idx), Pb_lap_sim(idx), dev_lp);
end
%% ===7. Rayleigh Outage Probability Analysis =====
C_th = 1.2; % Required capacity (bps/Hz)
gamma_th = 2^C_th - 1; % SNR threshold for outage
Pout_ray_theory = zeros(size(EbN0_lin));
Pout_ray_sim = zeros(size(EbN0_lin));

for idx = 1:numel(EbN0_lin)
gamma_bar = EbN0_lin(idx);

% Theoretical outage probability for Rayleigh fading:
Pout_ray_theory(idx) = 1 - exp(-gamma_th/gamma_bar);

% Simulation
h = (randn(1, numBits) + 1j*randn(1, numBits))/sqrt(2);
gamma_inst = abs(h).^2 * gamma_bar;
outage_event = log2(1 + gamma_inst) < C_th;
Pout_ray_sim(idx) = mean(outage_event);
end

figure;
semilogy(EbN0_dB, Pout_ray_theory, 'k-', 'LineWidth',2); hold on;
semilogy(EbN0_dB, Pout_ray_sim, 'ko--','LineWidth',1.5);
xlabel('E_b/N_0 (dB)');
ylabel('Outage Probability');
title(sprintf('Rayleigh Fading Channel Outage Probability (Capacity = %.1f bps/Hz, SNR_{th}=%.2f)', C_th, gamma_th));
legend('Theory','Simulation','Location','SouthWest');
grid on;
%% ======= 8.Rician fading and AWGN channels ==========
%% === parameters ===
U=1;
R = 3;
hd = 1;
alpha = 2.2;
N0 = 1;
K = 5;

Pb = zeros(size(EbN0_dB)); %  BER
Pout_theory = zeros(size(EbN0_lin));
Pout_sim = zeros(size(EbN0_lin));

for i = 1:length(EbN0_lin)
Pb(i) = ber_rician_theory(EbN0_lin(i), K);
gammaBar = EbN0_lin(i);
Pout_theory(i) = outage_rician_theory(EbN0_lin(i), K, C_th);
end

% === Rician fading PDF ===
function f = pdf_rician(x, gammaBar, K)
A = (K + 1) / gammaBar;
B = exp( - K - A * x );
arg = 2 * sqrt( K * A * x );
f = A .* B .* besseli(0, arg);

% safe
f(~isfinite(f)) = 0;
f(x < 0) = 0;
end

% === theory BER  ===
function Pb = ber_rician_theory(gammaBar, K)
ric_pdf = @(x) pdf_rician(x, gammaBar, K);
Pb = integral(@(x) 0.5 .* erfc(sqrt(x)) .* ric_pdf(x), 1e-10, Inf, ...
'RelTol',1e-6, 'AbsTol',1e-12);
end
%=== simulation BER===
BER_sim = zeros(size(EbN0_lin));

for idx = 1:length(EbN0_lin)
EbN0 = EbN0_lin(idx);
% === BPSK  ===
bits = randi([0, 1], 1, numBits);
symbols = 2*bits - 1; % BPSK: 0 → -1, 1 → +1

% ===  Rician Fading Channel ===
mu = sqrt(K / (K + 1));
sigma = sqrt(1 / (2 * (K + 1)));
h = mu + sigma * (randn(1, numBits) + 1j * randn(1, numBits));

% === Noise ===
noise_var = 1 / (2 * EbN0);
noise = sqrt(noise_var) * (randn(1, numBits) + 1j * randn(1, numBits));

% === Receiver signal ===
rx = h .* symbols + noise;

% === Coherent Detection ===
rx_eq = real(rx ./ h);
bits_hat = rx_eq > 0;

% === BER Caculation ===
BER_sim(idx) = mean(bits_hat ~= bits);
end

%Theory outrange
function Pout = outage_rician_theory(gammaBar, K, C_th)
x = sqrt(2 * K);
y = sqrt(2 * (K + 1) * (2^C_th - 1) / gammaBar);
Pout = 1 - marcumq(x, y);
end
%Simulation outrange
Pout_sim = zeros(size(EbN0_lin));
for i = 1:length(EbN0_lin)
EbN0 = EbN0_lin(i);

% Rician Channel
mu = sqrt(K / (K + 1));
sigma = sqrt(1 / (2 * (K + 1)));
h = mu + sigma * (randn(1, numBits) + 1j * randn(1, numBits));

% Real-time Rate
SNR_inst = EbN0 * abs(h).^2;
rate_inst = log2(1 + SNR_inst);

% Outage 
Pout_sim(i) = mean(rate_inst < C_th);
end

figure('Color','w', 'Position',[100,100,900,600]); % 
subplot(1,2,1); hold on;
semilogy(EbN0_dB, Pb, 'b-', 'LineWidth', 2,'DisplayName', 'Theory (Rician)');
semilogy(EbN0_dB, BER_sim, 'ro--', 'LineWidth', 2, 'DisplayName', 'Simulation');
set(gca, 'YScale', 'log');
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('Location', 'southwest'); grid on; grid minor;
title(sprintf('Rician Fading BER theory(K = %d)', K));
grid on;

subplot(1,2,2); hold on;
semilogy(EbN0_dB, Pout_theory, 'b-', 'LineWidth',2, 'DisplayName', 'Theory');
semilogy(EbN0_dB, Pout_sim, 'ro--', 'LineWidth',2, 'DisplayName', 'Simulation');
set(gca, 'YScale', 'log');
xlabel('E_b/N_0 (dB)');
ylabel('Outage Probability');
legend('Location','southwest'); grid on; grid minor;
title(sprintf('Outage under Rician fading (K = %d, C_{th} = %.2f)', K, C_th));


%% ======= 9.randomly deploy users over Rayleigh fading and AWGN channels =========================


% Initialize the result matrix
Pout_ray_theory_user = zeros(U, length(EbN0_dB));
Pout_ric_theory_user = zeros(U, length(EbN0_dB));
Pout_awgn_theory = zeros(U, length(EbN0_dB));
Pout_ray = zeros(U, length(EbN0_lin));
Pout_awgn = zeros(U, length(EbN0_lin));

Pb_ray_theory_user = zeros(U, length(EbN0_dB));
Pb_ric_theory_user = zeros(U, length(EbN0_dB));
BER_ray_sim = zeros(U, length(EbN0_dB));
Pout_ric = zeros(U, length(EbN0_lin)); BER_ric_sim = zeros(U, length(EbN0_lin));

% Random User
theta = 2 * pi * rand(1,U);
r = R * sqrt(rand(1,U));
x = r .* cos(theta); y = r .* sin(theta);
d = sqrt(x.^2 + y.^2 + hd^2); % Distance
path_loss = d.^(-alpha); % Pathloss
K=5;
%% === Main Loop ===
for u = 1:U
for idx = 1:length(EbN0_lin)
snr_user = EbN0_lin(idx) * path_loss(u);
gammaBar_u = EbN0_lin(idx) * path_loss(u);
pdfA = @(x) (1/gammaBar_u) .* exp(-x/gammaBar_u) .* (x >= 0);
% Theory outage
Pout_ray_theory_user(u, idx) = integral(pdfA, 0, gamma_th, 'RelTol',1e-6);
% --- Rician PDF with Random User---
pdfG = @(x) real(min(1e100, ...
(K+1)/gammaBar_u .* ...
exp(-K - (K+1)*x./gammaBar_u) .* ...
besseli(0, 2*sqrt(K*(K+1)*x./gammaBar_u)) )) .* (x>=0);


% --- Simulation outage ---
Pout_ric_theory_user(u, idx) = integral(pdfG, 0, gamma_th, 'RelTol',1e-6);

% Theory  BER with Random User
Pb_ray_theory_user(u, idx) = integral(@(x) 0.5 .* erfc(sqrt(x)) .* pdfA(x), 0, Inf);
Pb_ric_theory_user(u, idx) = integral(@(x) 0.5 .* erfc(sqrt(x)) .* pdfG(x), 0, Inf);


%Simulation Rayleigh outage
h = (randn(1,numBits) + 1j*randn(1,numBits)) / sqrt(2);
rate_inst = log2(1 + snr_user * abs(h).^2);
Pout_ray(u, idx) = mean(rate_inst < C_th);



% ---------------BER Simulation --------------------
% --- Rayleigh ---
sigma_ray = sqrt(1/(2*snr_user));
noise_ray = sigma_ray * (randn(1,numBits) + 1j*randn(1,numBits));
rx_ray = h .* symbols + noise_ray;
rx_eq = real(rx_ray ./ h);
bits_hat_ray = rx_eq > 0;
BER_ray_sim(u, idx) = mean(bits_hat_ray ~= bits);

% ---------------Rician Simulation ------ 
% === Rician Channel ===
mu_ric = sqrt(K / (K + 1));
sigma_ric = sqrt(1 / (2 * (K + 1)));
h_ric = mu_ric + sigma_ric * (randn(1,numBits) + 1j*randn(1,numBits));

% === Simulation outage (Rician) ===
rate_ric = log2(1 + snr_user * abs(h_ric).^2);
Pout_ric(u, idx) = mean(rate_ric < C_th);

% === Simulation BER (Rician) ===
noise_ric = sqrt(1/(2 * snr_user)) * (randn(1, numBits) + 1j * randn(1, numBits));
rx_ric = h_ric .* symbols + noise_ric;
rx_eq_ric = real(rx_ric ./ h_ric);
bits_hat_ric = rx_eq_ric > 0;
BER_ric_sim(u, idx) = mean(bits_hat_ric ~= bits);

end
end

%% === Outage ===
figure; hold on;

for u = 1:U
semilogy(EbN0_dB, Pout_ray_theory_user(u,:), 'go-', ...
'DisplayName', sprintf('User %d Raytheory', u));
semilogy(EbN0_dB, Pout_ray(u, :), 'g*-', ...
'DisplayName', sprintf('User %d Raysim', u));
semilogy(EbN0_dB, Pout_ric_theory_user(u,:), 'b>-', ...
'DisplayName', sprintf('User %d Rictheory', u));
semilogy(EbN0_dB, Pout_ric(u,:), 'bx--', ...
'DisplayName', sprintf('User %d RicSim', u));
set(gca, 'YScale', 'log');

end
grid on; ylim([1e-2 1]);
xlabel('E_b/N_0 (dB)'); ylabel('Outage Probability');
title('Outage Probability vs E_b/N_0'); legend('Location','southwest');

%% === BER (Random User)===
figure('Color','w'); hold on;
colors = lines(U);
markers = {'o','s','^','d','v'};

% User
for u = 1:U
semilogy(EbN0_dB, Pb_ray_theory_user(u,:), 'g*-', ...
'DisplayName', sprintf('User %d - RayTheory', u));
semilogy(EbN0_dB, BER_ray_sim(u,:),'go-',...
'DisplayName', sprintf('User %d - RaySim', u));
semilogy(EbN0_dB, Pb_ric_theory_user(u,:), 'b>-', ...
'DisplayName', sprintf('User %d - RicTheory', u));
semilogy(EbN0_dB, BER_ric_sim(u,:), 'bx--', ...
'DisplayName', sprintf('User %d - RicSim', u));
set(gca, 'YScale', 'log');

end

xlabel('E_b/N_0 (dB)', 'FontSize',12);
ylabel('Bit Error Rate (BER)', 'FontSize',12);
title('BER vs E_b/N_0 for one random User', 'FontSize',14, 'FontWeight','bold');
legend('Location','northeast', 'FontSize',10);
grid on; grid minor;
ylim([1e-2 0.5]); xlim([min(EbN0_dB) max(EbN0_dB)]);
set(gca, 'FontSize',11);


%% === Random User(Color= Average Rayleigh outage)===
figure;
scatter(x, y, 120, mean(Pout_ray,2), 'filled');
colormap turbo; colorbar;
xlabel('X (m)'); ylabel('Y (m)');
title('User Deployment (color = Avg Rayleigh Outage)');
axis equal;
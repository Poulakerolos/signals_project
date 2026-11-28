% Phase 0: Q1 to Q5
clear; 
close all; 
clc;

%% Q1: Read both files.

[A, frA] = audioread('ResearcherA.wav');
[B, frB] = audioread('ResearcherB.wav');
infoA = audioinfo('ResearcherA.wav');
infoB = audioinfo('ResearcherB.wav');

fprintf('A: fr = %d Hz, %d samples, %d channel(s), %.2f s\n', frA, size(A,1), size(A,2), size(A,1)/frA);
fprintf('B: fr = %d Hz, %d samples, %d channel(s), %.2f s\n', frB, size(B,1), size(B,2), size(B,1)/frB);

% if a file is stereo: are the two channels the same (within one step)?
if size(A,2) == 2
    stepA = 2 / 2^infoA.BitsPerSample;                 % quantisation step
    fprintf('A: max difference between channels = %g (step = %g)\n', max(abs(A(:,1) - A(:,2))), stepA);
end
if size(B,2) == 2
    stepB = 2 / 2^infoB.BitsPerSample;
    fprintf('B: max difference between channels = %g (step = %g)\n', max(abs(B(:,1) - B(:,2))), stepB);
end

% keep the first channel only
A = A(:,1);
B = B(:,1);
fr = frA;                                              % both files have the same rate

%% ---------------- Q2: time axis, plot, speaking intervals ----------------
tA = (0:length(A)-1)' / fr;                            % time in seconds
tB = (0:length(B)-1)' / fr;

% Find who speaks: cut the signal into 0.1 s pieces and measure the power of each piece
pieceLen = 0.1;                                        % length of a piece in seconds
frame = round(pieceLen * fr);                          % samples per piece
x = {A, B};
active = cell(1,2);                                    % active{i}(m) = 1 if person i speaks in piece m

for i = 1:2
    nPieces = floor(length(x{i}) / frame);
    power = zeros(nPieces, 1);
    for m = 1:nPieces
        piece = x{i}((m-1)*frame+1 : m*frame);
        power(m) = mean(piece.^2);
    end
    active{i} = power > 0.01 * max(power);             % louder than 20 dB below the maximum

    % remove isolated single pieces (tiny blips / noise)
    for m = 2:nPieces-1
        if active{i}(m) && ~active{i}(m-1) && ~active{i}(m+1)
            active{i}(m) = false;
        end
    end
end

% Start and end time of every speaking interval
starts = cell(1,2);
ends   = cell(1,2);
for i = 1:2
    d = diff([0; active{i}; 0]);                       % +1 = speech starts, -1 = speech stops
    starts{i} = (find(d == 1) - 1) * pieceLen;         % seconds
    ends{i}   = (find(d == -1) - 1) * pieceLen;
end

% Do the intervals overlap?  (compare piece by piece)
n = min(length(active{1}), length(active{2}));
both = active{1}(1:n) & active{2}(1:n);                % 1 where both "speak"
overlapTime = sum(both) * pieceLen;                    % total overlap in seconds
fprintf('Total overlap = %.1f s out of %.2f s\n', overlapTime, length(A)/fr);
if overlapTime == 0
    disp('The speakers never overlap.');
else
    disp('Overlap found at these times (s):');
    disp(((find(both) - 1) * pieceLen)');
end

% Plot both signals and mark the speaking intervals with thick lines
figure; hold on; grid on;
plot(tA, A, 'b');
plot(tB, B, 'r');
for k = 1:length(starts{1})
    plot([starts{1}(k) ends{1}(k)], [1.05 1.05], 'b', 'LineWidth', 4);
end
for k = 1:length(starts{2})
    plot([starts{2}(k) ends{2}(k)], [1.12 1.12], 'r', 'LineWidth', 4);
end
xlabel('Time (s)'); ylabel('Amplitude');
title('Researcher A (blue) and B (red); thick lines = speaking');
legend('Researcher A', 'Researcher B');
ylim([-1 1.2]);

%% ---------------- Q3: spectrum with the FFT ----------------
figure;
names  = {'Researcher A', 'Researcher B'};
colors = {'b', 'r'};

for i = 1:2
    N = length(x{i});
    X = fft(x{i});

    half = X(1 : floor(N/2) + 1);                      % one-sided: 0 ... fr/2
    mag  = 2 * abs(half) / N;                          % amplitude
    magdB = 20 * log10(mag + eps);                     % in dB (eps avoids log(0))
    f = (0 : length(half)-1)' * fr / N;                % frequency axis in Hz

    fprintf('%s: frequency resolution = fr/N = %.4f Hz\n', names{i}, fr/N);

    subplot(2,1,1); hold on;
    plot(f, magdB, colors{i});

    subplot(2,1,2); hold on;
    plot(f, magdB, colors{i});
end

subplot(2,1,1); grid on; xlim([0 fr/2]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Spectrum, 0 to fr/2'); legend(names);

subplot(2,1,2); grid on; xlim([0 1000]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Zoom, 0 to 1 kHz'); legend(names);

%% ---------------- Q4: pitch of a 50 ms voiced segment ----------------
q4_showMarks = false;      % false = clean plots; true = mark the harmonic peaks used (to check)

segStart = [66.600, 173.075];                    % start time of each segment (s): A, B
Lseg = round(0.05 * fr);                         % 50 ms in samples
Nfft = 2^16;                                     % zero-padding for a smoother spectrum
w    = 0.5 - 0.5*cos(2*pi*(0:Lseg-1)'/(Lseg-1)); % Hann window
tms  = (0:Lseg-1)' / fr * 1000;                  % time axis in ms
fq   = (0:Nfft/2)' * fr / Nfft;                  % frequency axis in Hz
F0a  = zeros(1,2);                               % F0 from the period of the waveform
F0b  = zeros(1,2);                               % F0 from the spacing of the harmonics

figure;
for i = 1:2
    n0  = round(segStart(i) * fr) + 1;
    seg = x{i}(n0 : n0 + Lseg - 1);

    % ---- (a) period of the waveform ----
    % The period is the time shift at which the waveform best matches itself
    % (autocorrelation). We search for it between 70 Hz and 300 Hz.
    s  = seg - mean(seg);
    r  = real(ifft(abs(fft(s, 2*Lseg)).^2));     % r(lag+1) = similarity at that lag
    lagMin = round(fr/300);
    lagMax = round(fr/70);
    [~, k] = max(r(lagMin+1 : lagMax+1));
    lag = lagMin + k - 1;
    y1 = r(lag);  y2 = r(lag+1);  y3 = r(lag+2);
    lagRef = lag + 0.5*(y1 - y3)/(y1 - 2*y2 + y3);   % sharper peak position
    T0     = lagRef / fr;                            % period in s
    F0a(i) = 1 / T0;

    % ---- (b) spacing of the harmonics ----
    S   = fft(seg .* w, Nfft);
    SdB = 20*log10(2*abs(S(1:Nfft/2+1)) / sum(w) + eps);

    K  = floor(1000 / F0a(i));                   % harmonics below 1 kHz
    hk = (1:K)';                                 % harmonic numbers 1, 2, 3, ...
    hf = zeros(K,1);                             % their frequencies
    idx = zeros(K,1);
    for m = 1:K
        win = find(fq > (m-0.4)*F0a(i) & fq < (m+0.4)*F0a(i));   % look near m*F0
        [~, j] = max(SdB(win));                  % highest point = the harmonic peak
        idx(m) = win(j);
        hf(m)  = fq(idx(m));
    end
    c      = polyfit(hk, hf, 1);                 % straight line: frequency = F0 * number
    F0b(i) = c(1);                               % its slope is the spacing = F0

    fprintf('%s: (a) T0 = %.3f ms -> F0 = %.1f Hz | (b) harmonic spacing F0 = %.1f Hz\n', ...
            names{i}, T0*1000, F0a(i), F0b(i));

    % ---- plots (no data cursor labels) ----
    subplot(2,2,2*i-1);
    plot(tms, seg, colors{i}); grid on;
    xlabel('Time (ms)'); ylabel('Amplitude');
    title(sprintf('%s: 50 ms segment, T_0 = %.2f ms', names{i}, T0*1000));

    subplot(2,2,2*i);
    plot(fq, SdB, colors{i}); grid on; xlim([0 1500]); hold on;
    if q4_showMarks, plot(fq(idx), SdB(idx), 'ko'); end  %#ok<UNRCH>
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    title(sprintf('%s: spectrum, F_0 = %.1f Hz', names{i}, F0b(i)));
end
saveas(gcf, 'Q4_pitch.png');                     % saves the figure in the script's folder
%% ---------------- Q5: cumulative energy E(f) ----------------
figure; hold on; grid on;
h = zeros(1,2);
for i = 1:2
    N    = length(x{i});
    P    = abs(fft(x{i})).^2;                  % energy of every frequency bin
    half = floor(N/2) + 1;                     % bins from 0 Hz to fr/2
    wt   = 2 * ones(half, 1);                  % each positive bin also stands for its negative twin
    wt(1) = 1;                                 % 0 Hz has no twin
    if mod(N,2) == 0, wt(end) = 1; end         % neither does fr/2
    E    = 100 * cumsum(wt .* P(1:half)) / sum(P);     % cumulative energy in %
    f    = (0 : half-1)' * fr / N;             % frequency of each bin in Hz

    % the numbers asked for
    f99   = f(find(E >= 99, 1));               % frequency below which 99 % of the energy lies
    k34   = round(3400 * N / fr) + 1;          % bin number of 3.4 kHz
    k4    = round(4000 * N / fr) + 1;          % bin number of 4 kHz
    fprintf('%s: 99%% of the energy is below %.0f Hz\n', names{i}, f99);
    fprintf('%s: energy above 3.4 kHz = %.2f %%, above 4 kHz = %.2f %%\n', names{i}, 100-E(k34), 100-E(k4));
    fprintf('%s: check, E at fr/2 = %.4f %%\n', names{i}, E(end));

    % plot (only 2000 log-spaced points, otherwise the plot is very slow)
    kStart = ceil(50 * N / fr) + 1;
    pick   = unique(round(logspace(log10(kStart), log10(half), 2000)));
    h(i)   = plot(f(pick), E(pick), colors{i}, 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log');
plot([3400 3400], [0 100], 'k--');  text(3400, 8, ' 3.4 kHz');
plot([4000 4000], [0 100], 'k:');   text(4000, 18, ' 4 kHz');
xlim([50 fr/2]); ylim([0 100]);
xlabel('Frequency f (Hz)'); ylabel('E(f) (%)');
title('Cumulative energy'); legend(h, names, 'Location', 'southeast');
saveas(gcf, 'Q5_energy.png');
save('phase0_data.mat', 'A', 'B', 'fr');

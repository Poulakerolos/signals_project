
% Phase 0: Q1, Q2, Q3 
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

save('phase0_data.mat', 'A', 'B', 'fr');
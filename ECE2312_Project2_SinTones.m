%  ECE 2312 - Project 2
%  Matthew Madamba

%********************************************
%         Sine Tone Generation 
% *******************************************

% Produce a sin tone

% Define frequencies
Fs = 48000;            % Sampling rate from project 1 (Hz)
sine_frequency = 5000; % Given frequency of the sine tone (Hz)

% Sine tone duration to match recorded speech file - 5s duration
duration = 5; 

% Create time vector
t = linspace(0, duration, Fs * duration);

% Generate sine tone - sin(2pi*f*t)
sine_tone = sin(2 * pi * sine_frequency * t);

% Play the sine tone
sound(sine_tone, Fs);
pause(duration);  % Wait for the duration of the sine tone

% Save the sine tone as a WAV file
filename = 'team[MatthewMadamba]-sinetone.wav';
audiowrite(filename, sine_tone, Fs);


% Define parameters for the spectrogram
windowLength = 1024;                 % Length of the window for computing the spectrogram 
overlap = round(windowLength * 0.9); % Overlap 
nfft = 2048;                         % Number of FFT points

% Generate spectrogram of the sine tone
figure;
spectrogram(sine_tone, hamming(windowLength), overlap, nfft, Fs, 'yaxis'); 
title('Spectrogram of Sine Tone');
xlabel('Time (seconds)');
ylabel('Frequency (kHz)');
colorbar;
ylim([0 10]);  % Limit y-axis to 10 kHz


%*****************************************************************
%                   Chirp Signal Generation
%*****************************************************************

% Produce a chirp signal

% Varied frequency range
f0 = 0;               % Starting frequency (Hz)
f1 = 8000;            % Ending frequency (Hz)

% Generate Chirp signal using sin function
instantaneous_frequency = ((f1-f0)/(2*duration)).*t;
chirp_signal = sin(2*pi*instantaneous_frequency.*t);

% Play the chirp signal
sound(chirp_signal, Fs);
pause(duration);  % Wait for the duration of the chirp tone

% Save the chirp signal as a WAV file
filename = 'team[MatthewMadamba]-chirp.wav';
audiowrite(filename, chirp_signal, Fs);

% Plot the spectrogram of the chirp signal
windowLength = 1024;                 % Length of the window for computing the spectrogram 
overlap = round(windowLength * 0.9); % Overlap 
nfft = 2048;                         % Number of FFT points

figure;
spectrogram(chirp_signal, hamming(windowLength), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of Chirp Signal');
xlabel('Time (seconds)');
ylabel('Frequency (kHz)');
colorbar;
ylim([0 10]);  % Limit y-axis to 10 kHz


%*****************************************************************
%                     Fun With Sine Tone
%*****************************************************************

% Define frequencies for each tone in sequence
f_cetk_1 = 329.63 * 10;  % E
f_cetk_2 = 369.99 * 10;  % F#
f_cetk_3 = 293.66 * 10;  % D
f_cetk_4 = 146.83 * 10;  % D (one octave lower)
f_cetk_5 = 220.00 * 10;  % A

% Define the durations for each tone
duration_tone_1 = 0.75;  % Duration of tone 1 in seconds
duration_tone_2 = 0.9;  % Duration of tone 2 in seconds
duration_tone_3 = 1;  % Duration of tone 3 in seconds
duration_tone_4 = 0.9;  % Duration of tone 4 in seconds
duration_tone_5 = 1.75;  % Duration of tone 5 in seconds

% Define time vectors for each tone
t_tone_1 = linspace(0, duration_tone_1, Fs * duration_tone_1);
t_tone_2 = linspace(0, duration_tone_2, Fs * duration_tone_2);
t_tone_3 = linspace(0, duration_tone_3, Fs * duration_tone_3);
t_tone_4 = linspace(0, duration_tone_4, Fs * duration_tone_4);
t_tone_5 = linspace(0, duration_tone_5, Fs * duration_tone_5);

% Generate sine tone for each frequency
cetk_1 = sin(2 * pi * f_cetk_1 * t_tone_1);
cetk_2 = sin(2 * pi * f_cetk_2 * t_tone_2);
cetk_3 = sin(2 * pi * f_cetk_3 * t_tone_3);
cetk_4 = sin(2 * pi * f_cetk_4 * t_tone_4);
cetk_5 = sin(2 * pi * f_cetk_5 * t_tone_5);

% Concatenate all tones into one sequence
cetk_sequence = [cetk_1, cetk_2, cetk_3, cetk_4, cetk_5];

% Define the total duration of the concatenated sequence
total_duration = duration_tone_1 + duration_tone_2 + ...
                 duration_tone_3 + duration_tone_4 + ...
                 duration_tone_5;

% Play the concatenated sequence
sound(cetk_sequence, Fs);

% Save the concatenated sequence as a WAV file
filename = 'team[MatthewMadamba]-cetk_sequence.wav';
audiowrite(filename, cetk_sequence, Fs);

% Plot the spectrogram of the concatenated sequence
windowLength = 1024;                 % Length of the window for computing the spectrogram 
overlap = round(windowLength * 0.9); % Overlap 
nfft = 2048;                         % Number of FFT points

figure;
spectrogram(cetk_sequence, hamming(windowLength), overlap, nfft, Fs, 'yaxis');
title('Spectrogram of CETK Sequence');
xlabel('Time (seconds)');
ylabel('Frequency (kHz)');
colorbar;
ylim([0 10]);  % Limit y-axis to 10 kHz


%*****************************************************************
%                     Combining Sound Files
%*****************************************************************

% Load the speech signal
[speech_signal, Fs_speech] = audioread('the_quick_brown_fox.wav');

% Generate the sine tone with a frequency of 5000 Hz
duration_speech = length(speech_signal) / Fs_speech;
t_speech = linspace(0, duration_speech, length(speech_signal));
sine_frequency = 5000;
sine_tone = sin(2 * pi * sine_frequency * t_speech);

% Ensure both signals have the same length
sine_tone = sine_tone(1:length(speech_signal));

% Add the sine tone to the speech signal
speech_with_sine = speech_signal + sine_tone';

% Scale down the entire combined signal to prevent clipping
max_amplitude = max(abs(speech_with_sine));
scaling_factor = 0.9 / max_amplitude; % Adjust this factor as needed
speech_with_sine_scaled = scaling_factor * speech_with_sine;

% Play the resulting signal
sound(speech_with_sine, Fs_speech);

% Save the resulting signal as a WAV file
filename = 'team[MatthewMadamba]-speechchirp.wav';
audiowrite(filename, speech_with_sine, Fs_speech);

% Plot the spectrogram of the resulting signal
windowLength = 1024;                 % Length of the window for computing the spectrogram 
overlap = round(windowLength * 0.9); % Overlap 
nfft = 2048;                         % Number of FFT points

figure;
spectrogram(speech_with_sine, hamming(windowLength), overlap, nfft, Fs_speech, 'yaxis');
title('Spectrogram of Speech with Added Sine Tone');
xlabel('Time (seconds)');
ylabel('Frequency (kHz)');
colorbar;
ylim([0 10]);  % Limit y-axis to 10 kHz

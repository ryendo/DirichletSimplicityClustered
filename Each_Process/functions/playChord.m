function playChord(frequencies, duration)
    % playChord: Plays a chord consisting of specified frequencies and saves it to a file.
    %
    % Input:
    % frequencies - An array of frequencies that make up the chord (Hz)
    % duration - Duration of the playback (seconds)
    %
    % Example:
    % playChord([261.63, 329.63, 392.00], 2);

    % Set the sample rate
    Fs = 44100; % Sample rate (Hz)

    % Create the time axis
    t = 0:1/Fs:duration;

    % Generate and combine sample data for each frequency
    y = zeros(1, length(t));
    for i = 1:length(frequencies)
        y = y + sin(2*pi*frequencies(i)*t);
    end

    % Normalize the volume
    y = y / max(abs(y));

    % Play the chord
    sound(y, Fs);

    % Generate filename based on frequencies
    freq_str = sprintf('%.2f_', frequencies);
    freq_str = freq_str(1:end-1); % Remove trailing underscore
    filename = ['chord_' freq_str '.wav'];

    % Save the audio to a file
    % audiowrite(filename, y, Fs);
end

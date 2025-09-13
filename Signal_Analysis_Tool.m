function SignalAnalyzerOnline
    % --- Create Main Figure ---
    hFig = figure('Name','Signal Analyzer (Online GUI)', ...
                  'NumberTitle','off', ...
                  'Position',[100 100 900 600]);

    % --- Axes for plots ---
    ax1 = subplot(3,2,1,'Parent',hFig); title(ax1,'Original Signal');
    ax2 = subplot(3,2,2,'Parent',hFig); title(ax2,'Filtered Signal');
    ax3 = subplot(3,2,3,'Parent',hFig); title(ax3,'Manual DFT');
    ax4 = subplot(3,2,4,'Parent',hFig); title(ax4,'Recursive FFT');
    ax5 = subplot(3,2,5,'Parent',hFig); title(ax5,'Z-Plane of FIR Filter');

    % --- Buttons ---
    uicontrol('Style','pushbutton','String','🎙️ Record Signal', ...
              'Position',[50 20 120 40], ...
              'Callback',@recordSignal);

    uicontrol('Style','pushbutton','String','📊 Analyze Signal', ...
              'Position',[200 20 120 40], ...
              'Callback',@analyzeSignal);

    uicontrol('Style','pushbutton','String','🔊 Play Audio', ...
              'Position',[350 20 120 40], ...
              'Callback',@playAudio);

    % --- Label (Classification Result) ---
    hLabel = uicontrol('Style','text','String','Awaiting Signal...', ...
                       'Position',[500 20 300 40], ...
                       'FontSize',12,'FontWeight','bold');

    % --- Global Data ---
    Fs = 8000;      % Sampling rate
    Duration = 3;   % Record duration
    x = [];
    filteredSignal = [];
    filterCoeffs = [];

    % ================= Callback Functions =================
    function recordSignal(~,~)
        rec = audiorecorder(Fs, 16, 1);
        recordblocking(rec, Duration);
        x = getaudiodata(rec)'; % row vector
        t = (0:length(x)-1)/Fs;
        plot(ax1, t, x);
        xlabel(ax1,'Time (s)'); ylabel(ax1,'Amplitude');
        grid(ax1,'on');
        title(ax1,'Original Audio Signal');
    end

    function analyzeSignal(~,~)
        if isempty(x)
            warndlg('Please record a signal first.','Warning');
            return;
        end

        % Step 1: FIR Filter
        N = 50; fc = 1000;
        normFc = fc/(Fs/2);
        h = fir1(N, normFc, hamming(N+1));
        filterCoeffs = h;

        % Step 2: Convolution
        filteredSignal = manualConvolution(x,h);
        t = (0:length(filteredSignal)-1)/Fs;
        plot(ax2,t(1:length(x)),filteredSignal(1:length(x)));
        xlabel(ax2,'Time (s)'); ylabel(ax2,'Amplitude');
        title(ax2,'Filtered Signal'); grid(ax2,'on');

        % Step 3: Manual DFT
        dftInput = filteredSignal(1:min(512,end));
        [magDFT,~,fDFT] = manualDFT(dftInput,Fs);
        plot(ax3,fDFT,magDFT);
        xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'|X[k]|');
        title(ax3,'Manual DFT'); grid(ax3,'on');

        % Step 4: Recursive FFT
        fftInput = dftInput;
        if length(fftInput)<512
            fftInput = [fftInput,zeros(1,512-length(fftInput))];
        end
        Xfft = recursiveFFT(fftInput);
        f_fft = (0:511)*(Fs/512);
        plot(ax4,f_fft,abs(Xfft));
        xlabel(ax4,'Frequency (Hz)'); ylabel(ax4,'|X[k]|');
        title(ax4,'Recursive FFT'); grid(ax4,'on');

        % Step 5: Z-Plane
        z = roots(h);
        cla(ax5); hold(ax5,'on');
        theta = linspace(0,2*pi,300);
        plot(ax5,cos(theta),sin(theta),'--k'); % unit circle
        plot(ax5,real(z),imag(z),'bo','MarkerFaceColor','b');
        title(ax5,'Z-Plane'); xlabel(ax5,'Real'); ylabel(ax5,'Imag');
        axis(ax5,'equal'); grid(ax5,'on'); hold(ax5,'off');

        % Step 6: Frequency Classification
        fftMag = abs(Xfft(1:256));
        f_half = f_fft(1:256);
        [pks,locs] = findpeaks(fftMag,'MinPeakHeight',max(fftMag)/10);
        if isempty(pks)
            set(hLabel,'String','⚠️ No dominant frequency found.');
        else
            domFreq = f_half(locs(1));
            if domFreq<300
                set(hLabel,'String','🟢 Low Frequency (Male Voice)');
            elseif domFreq<800
                set(hLabel,'String','🟡 Mid Frequency (Female Voice)');
            elseif domFreq<2000
                set(hLabel,'String','🔴 High Frequency (Instrument)');
            else
                set(hLabel,'String','⚠️ Very High Frequency (Noise)');
            end
        end
    end

    function playAudio(~,~)
        if isempty(filteredSignal)
            warndlg('No filtered audio found. Analyze first.','Error');
            return;
        end
        sound(filteredSignal,Fs);
    end

    % ================= Utility Functions =================
    function y = manualConvolution(x,h)
        Nx = length(x); Nh = length(h);
        y = zeros(1,Nx+Nh-1);
        for n=1:length(y)
            for k=1:Nh
                if (n-k+1>0)&&(n-k+1<=Nx)
                    y(n) = y(n)+h(k)*x(n-k+1);
                end
            end
        end
    end

    function [mag,phase,f] = manualDFT(x,Fs)
        N = length(x);
        X = zeros(1,N);
        for k=1:N
            for n=1:N
                X(k) = X(k)+x(n)*exp(-1j*2*pi*(k-1)*(n-1)/N);
            end
        end
        mag = abs(X(1:N/2));
        phase = angle(X(1:N/2));
        f = (0:N/2-1)*(Fs/N);
    end

    function X = recursiveFFT(x)
        N = length(x);
        if N==1
            X=x;
        else
            even = recursiveFFT(x(1:2:end));
            odd = recursiveFFT(x(2:2:end));
            X = zeros(1,N);
            for k=1:N/2
                W = exp(-1j*2*pi*(k-1)/N);
                X(k) = even(k)+W*odd(k);
                X(k+N/2) = even(k)-W*odd(k);
            end
        end
    end
end

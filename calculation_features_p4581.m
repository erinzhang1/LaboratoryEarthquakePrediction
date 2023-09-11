runacpath='D:\UCL\MSc_proj\original_data\p4581acfiles\run1\WF_'; % path where the .ac files live

% acoustic parameters
load('ACTimeVector_4581.mat'); % load in the acoustic time vector data;
ts=ts_adjusted; % offset for time vector

acSettings = load('p4581_AE.mat');              % load acoustic settings
numSFpfile = acSettings.numFrames/2;            % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;               % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;         % number of WF per file and per channel
numCH = length(acSettings.channels2save);       % number of channels
WFlength = acSettings.Nsamples;                 % segment length
ts = ts/1e6;                                    % from microsec to sec
fs = 1/ts;                                      % acoustic sampling rate
%clear acSettings

% time vector for each waveform
timeWF = (0:WFlength-1)'*ts;
AETime = NaN(WFlength*numWFpfilepCH,1);

base = 3800;
result = [["time", "var1", "var2", "skew1", "skew2", "kurt1", "kurt2", "entropy1", "entropy2", "entropy change1", "entropy change2"]];
num = 0;
H1_prev = 0.000000000001;
H2_prev = 0.000000000001;

idxAcTime = find(acTime > base,1,'first'); %find the next index over from the time at which you want to analyze   
    
filenumber = ceil(idxAcTime/numWFpfilepCH); % file number for WF   
     
% build time vector for each file   
for jj = 1:numWFpfilepCH
    AETime((jj-1)*WFlength+1:jj*WFlength) = acTime((filenumber-1)*numWFpfilepCH+jj) + timeWF;
end    
ACfilename = [runacpath num2str(filenumber) '.ac'];
fid = fopen(ACfilename,'r');
ACdata = fread(fid,'int16');   
fclose(fid);
% reshape to get one column per channel
ACdata = reshape((ACdata),[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel

AETime_all = AETime;
ACdata_all = ACdata;
pointer = 1;
length = floor(1/((AETime(end) - AETime(1))/size(ACdata,1)));
shift = floor(length/10);

while num < 3100
    disp(num)
    size_AC_left = size(ACdata_all(pointer:end, :));
    if size_AC_left(1) >= length
        windows = ACdata_all(pointer: pointer + length, :);
        time = AETime_all(pointer + floor(length/2));
        pointer = pointer + shift + 1;
    else
        temp = ACdata_all(pointer:end, :);
        temp_time = AETime_all(pointer:end, :);
        filenumber = filenumber + 1;
        for jj = 1:numWFpfilepCH
            AETime((jj-1)*WFlength+1:jj*WFlength) = acTime((filenumber-1)*numWFpfilepCH+jj) + timeWF;
        end    
            
        ACfilename = [runacpath num2str(filenumber) '.ac'];
        fid = fopen(ACfilename,'r');
        ACdata = fread(fid,'int16');   
        fclose(fid);
        % reshape to get one column per channel
        ACdata = reshape((ACdata),[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
        ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
        ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel
        
        ACdata_all = [temp; ACdata];
        AETime_all = [temp_time; AETime];
        time = AETime_all(floor(length/2));
        windows = ACdata_all(1: length, :);
        pointer = 1 + shift;
    end
    % calculate entropy
    symbols1 = unique(windows(:, 1));
    probabilities1 = histcounts(windows(:, 1), [symbols1; max(symbols1)+1]) / numel(windows(:, 1));
    H1 = -sum(probabilities1 .* log2(probabilities1));

    symbols2 = unique(windows(:, 2));
    probabilities2 = histcounts(windows(:, 2), [symbols2; max(symbols2)+1]) / numel(windows(:, 2));
    H2 = -sum(probabilities2 .* log2(probabilities2));
    
    H1_delta = abs((H1 - H1_prev)/H1_prev) * 100;
    H2_delta = abs((H2 - H2_prev)/H2_prev) * 100;
    result = [result; [time, var(windows), skewness(windows), kurtosis(windows), H1, H2, H1_delta, H2_delta]];
    num = num + 1;
    H1_prev = H1;
    H2_prev = H2;
end
writematrix(result, "variance_p4581_old.csv")
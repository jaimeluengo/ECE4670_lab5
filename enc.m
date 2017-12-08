clear



%% Parametres to tweak
numbits=2e5; %length of bitsteam to send
P = 0.00125;%Average power constraint
bits=randi([0 1],numbits,1); %bitstream
M=128; % modulation index of higher gain frequency band
N=16;% lower modulation index of lower gain frequency bands
n_data_symble=1850; %number of samples of data per symbol
n_lowqam=25;%number of samples with N modulation index at low frequency
n_highqam=125;%number of samples with N modulation index at high frequency
n_lowf=172; %number of samples zero-padded in the low frequency band
n_highf=1058;%number of sa`mples zero-padded in the high frequency band
n_prefix=120;%number of samples in the cyclic prefix
gamma = P*30;%scaled power for data
gammat= P*34;%scaled power for training

%% Subsequent parametres
symbol_size = (1+2*(n_highf+n_lowf+n_data_symble));% number of samples per OFDM symbol(without cyclic prefix)
n_prime=n_data_symble-n_highqam-n_lowqam; %number of samples with M modulation index at high gain frequency band
n_symbols = numbits/(n_prime*log2(M)+(n_highqam+n_lowqam)*log2(N));%number of 'data' OFDM symbles
n_tsymbols =double(uint8(n_symbols/3));%number of training OFDM symbles
block_size= numbits/n_symbols;%number of bits
%the relationship between number of bits per sample is determined by the
%modulation index of qam



%% Pre-allocate memory
%to do


%% Generation of training symbles
% t=[0:1/4000:1-1/4000];
 t=rand(n_data_symble,1);
 rand_realizations =ones(n_data_symble,1).*exp(j*t*2*pi); %same module but
 %randomized phases

% rand_realizations= ones(n_data_symble,1);

% rand_realizations=[1:n_data_symble]';

x2_train=[zeros(n_lowf,1);rand_realizations;zeros(n_highf,1)];%zero-padding
x3_train=[0;x2_train;conj(fliplr(x2_train')')];%prepend DC 0 and postpend flipped conjugate
x4_train=ifft(x3_train)*sqrt(length(x3_train)); %make training symble real by taking ifft
x5_train=[x4_train(end-n_prefix+1:end);x4_train]*gamma;%satisfy power constraint
x=[];%initialize output vector

x0_p=[];%delete later

%% encoding OFDM symbols
for i=1:n_symbols %fills output vector symbol by symble
    x7_prime=[];
    if(mod(i,3)==0)%introduce training symbol at positions multiple of 3
        x=[x; x5_train];
    end
    x0=bits((i-1)*block_size+1:i*block_size); %extract bits to send
    
    %encode bits using qam with different modulation index depending on the
    %gain at the frequency that they are encoded to.
    
    for k=1:n_lowqam %takes bits in groups of log2(N) and converts them to decimals
        x0_low(k,1)=bi2de(fliplr(x0((k-1)*log2(N)+1:k*log2(N))'));
    end
    for k=1:n_prime %same but groups of log2(M)
        x0_prime(k,1)=bi2de(fliplr(x0(n_lowqam*log2(N)+(k-1)*log2(M)+1:n_lowqam*log2(N)+k*log2(M))'));
    end
    for k=1:n_highqam
         x0_high(k,1)=bi2de(fliplr(x0(n_lowqam*log2(N)+n_prime*log2(M)+(k-1)*log2(N)+1:n_lowqam*log2(N)+n_prime*log2(M)+k*log2(N))'));
    end

    x7_prime=[x7_prime;x0_prime];% delete later

    
    % encodes data symbol in the form of decimal numbers using QAM
    x1=[qammod(x0_low,N)]/log2(N);
    x1=[x1;qammod(x0_prime,M)/log2(M)];
    x1=[x1;qammod(x0_high,N)/log2(N)];
    x1=x1*gamma;% satisfy power constraint
    x1_rand=x1.*rand_realizations; %avoid clipping

    % generate OFDM symbol as in training
    x2=[zeros(n_lowf,1);x1_rand;zeros(n_highf,1)];
    x3=[0;x2;conj(fliplr(x2')')];
    x4=ifft(x3)*sqrt(length(x3));
    x0_p=[x0_p;x0_prime];
    x5=[x4(end-n_prefix+1:end);x4];
    x=[x;x5];
end

% check if power constraint is violated
avg_power = (x'*x)/length(x);
avg_pwr_flag = (avg_power <= P);

Fs=44100;
audiowrite('tx.wav',x,Fs);

%% Fake channel
% load('IR0.mat', 'impulse');
% h = impulse';
% y = conv(x, impulse);
% % no noise
% y = y(1:length(x));
% % plus random silence
% % create random L0 and L1 pauses from uniform distribution 
% L0 = round(1000*rand(1,1));
% L1 = round(1000*rand(1,1));
% y = [zeros(L0,1) ; y ; zeros(L1,1)];
% % % add noise
% % y = y + 0.0001*randn(length(y),1);



% actual channel:
create random L0 and L1 pauses from uniform distribution 
a = 0.25;
b = 3;
L0 = (b-a).*rand(1,1) + a;
L1 = (b-a).*rand(1,1) + a;
% change these lines before submission!!!!!
cmd1 = sprintf('~/Desktop/ECE4670/ccplay/ccplay --prepause ');
cmd2 = sprintf('%f --postpause %f --channel audio0 tx.wav rx.wav', L0, L1);
cmd = [cmd1 cmd2];
system(cmd);

%% Decoder
y = audioread('rx0.wav');

%% trigger to avoid prepause of receiver
square = y.^2;
W = 10;
s = movmean(square, W);
len_y = length(y);
thresh = 1e-6;
tau = -1;
for k = 1:len_y
    if (abs(s(k))>thresh)
        start = k-tau;  
        break;
    end
end
L1=start;


stem(y)
disp(start)

hold on
plot(start,0,'-o')
len_sent=length(x);%hardcode later!!
y = y(start:(start+len_sent-1));%removed silence

% remember to normalize by gamma


%% training
indexes=[];%indexes used to delete training indexes later from y vector
for (u=1:n_tsymbols)
    y0_train(:,u) = y((2+4*(u-1))*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*(3+4*(u-1)));
    %since we are seniding training symbols at positions multiple of three
    %we need to extract them in the form of [2+4*(u-1):3+4*(u-1)
    indexes=[indexes;[(2+4*(u-1))*(symbol_size+n_prefix)+1:(n_prefix+symbol_size)*(3+4*(u-1))]'];
    % Convert training symbol to its original vector doing the fft and
    % extraction
    y1_train (:,u) = fft(y0_train(:,u))/sqrt(symbol_size);
    y2_train (:,u) = y1_train(n_lowf+2:n_data_symble+1+n_lowf,u)/gamma;
    
    %Gain of the channel is output/input in absolute value
    gain (:,u) = abs(y2_train(:,u)./rand_realizations);% vector channel gain
    phase_y=angle(y2_train(:,u));
    phase_r=angle(rand_realizations);
    %Phase shift of the channel is angle(output)-angle(input)
    phase (:,u)= wrapToPi(phase_y-fliplr(phase_r)); %mantain between [-pi pi]
end

y(indexes)=[]; %delete training indexes from y
outbits = [];

%generate a vector of positions of data symbols until which a corresponding
%data symbol with index f is used
for k=1:n_tsymbols
    if(k==n_tsymbols)
        u(k)=4+3*k;%The last training symbol is used for all the remaing data symbols
    else
        u(k)=1+3*k;%index of last data symbol for which traing symbol k is used
    end
end
f=1;%initialization of training symble index

%% Decoding
for i=1:n_symbols
    if(i>=u(f))
        f=f+1; %f is the training symbol used for decoding at each time
    end
    y3=[];
    %more indexing... used for extraction of data symble without cyclic
    %prefix
    a=n_prefix+(n_prefix+symbol_size)*(i-1)+1;
    b=(n_prefix+symbol_size)*i;
    y0 = y(a:b);
    %reverse operations as of encoding
    y1 = fft(y0)/sqrt(symbol_size);
    y2 = y1(n_lowf+2:n_lowf+n_data_symble+1);
    y3=y2./gain(:,f).*exp(-j*phase(:,f))./rand_realizations/gamma;
    %separate samples into its corresponding qam frequency bands
    y3_lowqam=y3(1:n_lowqam)*log2(N);
    y3_prime=y3(n_lowqam+1:n_data_symble-n_highqam)*log2(M);
    y3_highqam=y3(n_data_symble-n_highqam+1:end)*log2(N);
    
    %for each subset of samples decode using qamdemod and convert to binary
    %dont forget binary base (log2(M) or log2(N)
    z_lowqam=qamdemod(y3_lowqam,N);
    for k=1:n_lowqam
        y4_lowqam(log2(N)*(k-1)+1:log2(N)*k,1)=(fliplr(de2bi(z_lowqam(k),log2(N))))';
    end
    outbits=[outbits;y4_lowqam];

     z_prime=qamdemod(y3_prime,M);
    for k=1:n_prime
        y4_prime(log2(M)*(k-1)+1:(log2(M))*k,1)=(fliplr(de2bi(z_prime(k),log2(M))))';
    end
    outbits=[outbits;y4_prime];

    z_highqam=qamdemod(y3_highqam,N);
    for k=1:n_highqam
        y4_highqam(log2(N)*(k-1)+1:log2(N)*k,1)=(fliplr(de2bi(z_highqam(k),log2(N))))';
    end
    outbits=[outbits;y4_highqam];
end


% Count correct ones and plot the errors
correct=sum(outbits==bits);
d=find(outbits~=bits);
t=length(d);
wrong=zeros(t,1);
wrong(d)=1;
figure(1)
plot(wrong)

%plot channel gain vs qam allocation and zero-padding
figure(2)
scale=(n_lowf+n_data_symble+n_highf)/length(gain(:,1));
alloc=[zeros(uint64(n_lowf/scale),1);ones(uint64(n_lowqam/scale),1);2*ones(uint64(n_prime/scale),1);ones(uint64(n_highqam/scale),1);zeros(uint64(n_highf/scale),1)];
plot(alloc);
hold on
plot(gain(:,1));

% FIgure of merit
Ne=numbits-correct;
Pr=max(1,800*avg_power);
data_rate=2e5*44.1e3/length(x);
perf=(data_rate*(1-Ne/1e5)^10)/Pr;


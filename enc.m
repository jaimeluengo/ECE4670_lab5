function x = enc(bits)

%% Parameters to tweak
numbits = 2e5;  %length of bitsteam to send
M = 64;  % modulation index of higher gain frequency band
N = 16;  % lower modulation index of lower gain frequency bands
n_data_symble = 2125;  %number of samples of data per symbol
n_lowqam = 50;  % number of samples with N modulation index at low frequency
n_highqam = 75;  % number of samples with N modulation index at high frequency
n_lowf = 112;  % number of samples zero-padded in the low frequency band
n_highf = 908;  % number of sa`mples zero-padded in the high frequency band
n_prefix = 120;  % number of samples in the cyclic prefix

%% Subsequent parameters
symbol_size = (1+2*(n_highf+n_lowf+n_data_symble));  % number of samples per OFDM symbol(without cyclic prefix)
n_prime = n_data_symble-n_highqam-n_lowqam;  % number of samples with M modulation index at high gain frequency band
n_symbols = numbits/(n_prime*log2(M)+(n_highqam+n_lowqam)*log2(N));  % number of 'data' OFDM symbles
n_tsymbols = double(uint8(n_symbols/3));  % number of training OFDM symbles
block_size = numbits/n_symbols;  %number of bits
% the relationship between number of bits per sample is determined by the
% modulation index of QAM
P = 0.00125;  % Average power constraint
gamma = P*30;  % scaled power
gammat = P*40;

%% Generation of training symbols
rng(17);
t = rand(n_data_symble,1);
rand_realizations = ones(n_data_symble,1).*exp(1j*t*2*pi);  % same module but

% zero-padding
x2_train = [zeros(n_lowf,1) ; rand_realizations ; zeros(n_highf,1)];
% prepend DC 0 and postpend flipped conjugate
x3_train = [0;x2_train ; conj(fliplr(x2_train')')];
% make training symbol real by taking ifft
x4_train = ifft(x3_train)*sqrt(length(x3_train)); 
% satisfy power constraint
x5_train = [x4_train(end-n_prefix+1:end);x4_train]*gammat;
%initialize output vector
x=[];
%delete later
x0_p=[];

%% encoding OFDM symbols
for i = 1:n_symbols  % fills output vector symbol by symbol
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


Fs=44100;
audiowrite('tx.wav',x,Fs);

return 

function outbits = dec()

clear all
%% Parametres to tweak
numbits=2e5; %length of bitsteam to send
M=64; % modulation index of higher gain frequency band
N=16;% lower modulation index of lower gain frequency bands
n_data_symble=2125; %number of samples of data per symbol
n_lowqam=50;%number of samples with N modulation index at low frequency
n_highqam=75;%number of samples with N modulation index at high frequency
n_lowf=112; %number of samples zero-padded in the low frequency band
n_highf=908;%number of sa`mples zero-padded in the high frequency band
n_prefix=120;%number of samples in the cyclic prefix

%% Subsequent parametres
symbol_size = (1+2*(n_highf+n_lowf+n_data_symble));% number of samples per OFDM symbol(without cyclic prefix)
n_prime=n_data_symble-n_highqam-n_lowqam; %number of samples with M modulation index at high gain frequency band
n_symbols = numbits/(n_prime*log2(M)+(n_highqam+n_lowqam)*log2(N));%number of 'data' OFDM symbles
n_tsymbols =double(uint8(n_symbols/3));%number of training OFDM symbles
block_size= numbits/n_symbols;%number of bits
%the relationship between number of bits per sample is determined by the
%modulation index of qam
P = 0.00125;%Average power constraint
gamma = P*30;%scaled power
gammat = P*40;

rng(17);
t=rand(n_data_symble,1);
rand_realizations =ones(n_data_symble,1).*exp(j*t*2*pi); %same module but

%%

y = audioread('rx.wav');


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


% stem(y)
% disp(start)
% 
% hold on
% plot(start,0,'-o')
len_sent=length(x);%hardcode later!!
y = y(start:(start+len_sent-1));%removed silence

% remember to normalize by gamma


%% training
indexes=[];%indexes used to delete training indexes later from y vector
for (u=1:n_tsymbols)
    y0_train(:,u) = y((2+5*(u-1))*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*(3+5*(u-1)));
    %since we are seniding training symbols at positions multiple of three
    %we need to extract them in the form of [2+5*(u-1):3+5*(u-1)
    indexes=[indexes;[(2+5*(u-1))*(symbol_size+n_prefix)+1:(n_prefix+symbol_size)*(3+5*(u-1))]'];
    % Convert training symbol to its original vector doing the fft and
    % extraction
    y1_train (:,u) = fft(y0_train(:,u))/sqrt(symbol_size);
    y2_train (:,u) = y1_train(n_lowf+2:n_data_symble+1+n_lowf,u)/gammat;
    
    %Gain of the channel is output/input in absolute value
    gain (:,u) = abs(y2_train(:,u)./rand_realizations);% vector channel gain
    
    phase_y=angle(y2_train(:,u));
    phase_r=angle(rand_realizations);
    x_ro=[1:length(n_data_symble)]';
    
    %phase has clock skew
    for(t=1:2)
        if(max(gain(:,1))<0.6)
            Pol=[7.9859e-04;0.075]*11*t/12;
        else
            Pol=[3.3459e-04;0.0025]*11*t/12;
        end
%         Pol=[6.65058823529413e-05;0.007416666666667]*t;
        phase_dif=Pol(1)*x_ro+Pol(2);
        phase(:,2+4*(u-1)+1-t)=wrapToPi(phase_y-fliplr(phase_r)+phase_dif);
        phase(:,2+4*(u-1)+t)=wrapToPi(phase_y-fliplr(phase_r)-phase_dif);
    end

end

y(indexes)=[]; %delete training indexes from y
outbits = [];

%generate a vector of positions of data symbols until which a corresponding
%data symbol with index f is used
f=1;%initialization of training symble index
g=1; %data symbol position
%% Decoding
for i=1:n_symbols
    if(i~=1&&mod(i-1,4)==0)
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
    y3=y2./gain(:,f).*exp(-j*phase(:,i))./rand_realizations/gamma;

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
    g=g+1;
end


return

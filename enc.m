clear
numbits=2e5;
bits=randi([0 1],numbits,1);
M=32; % modulation index
N=16;% lower modulation index
n_data_symble=2625; %number of samples of data per symbol
n_lowqam=25;
n_highqam=600;
n_lowf=112; %number of samples in the low frequency band
n_highf=708;%number of samples in the low frequency band
n_prefix=120;%number of samples in the cyclic prefix
symbol_size = (1+2*(n_highf+n_lowf+n_data_symble));
n_symbols = numbits/((n_data_symble-n_lowqam-n_highqam)*log2(M)+(n_highqam+n_lowqam)*log2(N));
n_tsymbols =double(uint8(n_symbols/3));
block_size= numbits/n_symbols;
P = 0.00125;
gamma = P*34;


% t=[0:1/4000:1-1/4000];
 t=rand(n_data_symble,1);
 rand_realizations =ones(n_data_symble,1).*exp(j*t*2*pi);

% rand_realizations= ones(n_data_symble,1);

% rand_realizations=[1:n_data_symble]';

x2_train=[zeros(n_lowf,1);rand_realizations;zeros(n_highf,1)];
x3_train=[0;x2_train;conj(fliplr(x2_train')')];
x4_train=ifft(x3_train)*sqrt(length(x3_train));
x5_train=[x4_train(end-n_prefix+1:end);x4_train]*gamma;
x=[];
x0_p=[];
for i=1:n_symbols
    x7_prime=[];
    if(mod(i,3)==0)
        x=[x; x5_train];
    end
    x0=bits((i-1)*block_size+1:i*block_size);   %*sqrt(2);
    for k=1:n_lowqam
        x0_low(k,1)=bi2de(fliplr(x0((k-1)*(log2(M)-1)+1:k*(log2(M)-1))'));
    end
    for k=1:n_data_symble-n_highqam-n_lowqam
        x0_prime(k,1)=bi2de(fliplr(x0(n_lowqam*log2(N)+(k-1)*log2(M)+1:n_lowqam*log2(N)+k*log2(M))'));
    end
    for k=1:n_highqam
         x0_high(k,1)=bi2de(fliplr(x0(n_lowqam*log2(N)+(n_data_symble-n_highqam-n_lowqam)*log2(M)+(k-1)*log2(N)+1:n_lowqam*log2(N)+(n_data_symble-n_highqam-n_lowqam)*log2(M)+k*log2(N))'));
    end

    x7_prime=[x7_prime;x0_prime];
    x1=[qammod(x0_low,2^(log2(M)-1))]/(log2(M)-1);
    x1=[x1;qammod(x0_prime,M)/log2(M)];
    x1=[x1;qammod(x0_high,2^(log2(M)-1))/(log2(M)-1)];
    x1=x1*gamma;
    x1_rand=x1.*rand_realizations;
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

%% channel
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
% % add noise
% % y = y + 0.0001*randn(length(y),1);



% actual channel:
% create random L0 and L1 pauses from uniform distribution 
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
y = audioread('rx.wav');


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
%disp(start)

% hold on
% plot(start,0,'-o')
len_sent=length(x);%hardcode later!!
y = y(start:(start+len_sent-1));%removed silence

% remember to normalize by gamma


%training
indexes=[];
for (u=1:n_tsymbols)
    y0_train(:,u) = y((2+4*(u-1))*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*(3+4*(u-1)));
    indexes=[indexes;[(2+4*(u-1))*(symbol_size+n_prefix)+1:(n_prefix+symbol_size)*(3+4*(u-1))]'];
    y1_train (:,u) = fft(y0_train(:,u))/sqrt(symbol_size);
    y2_train (:,u) = y1_train(n_lowf+2:n_data_symble+1+n_lowf,u)/gamma; % 113:41112

    gain (:,u) = abs(y2_train(:,u)./rand_realizations);% vector channel gain
    phase_y=angle(y2_train(:,u));
    phase_r=angle(rand_realizations);
    phase (:,u)= wrapToPi(phase_y-fliplr(phase_r)); %vector of channel phase shift
end

y(indexes)=[];
% ya=y(1:2*(symbol_size+n_prefix));
% for i=1:n_tsymbols
%     if(i==n_tsymbols)
%         ya=[ya;y((n_prefix+symbol_size)*(3+4*(i-1))+1:(n_prefix+symbol_size)*(4+4*(i-1)))];
%     else
%          ya=[ya;y((n_prefix+symbol_size)*(3+4*(i-1))+1:(6+4*(i-1))*(symbol_size+n_prefix))];
%     end
% end
outbits = [];

% load('phase_chan.mat','phase_chan');
% phase=phase_chan;
% stem(y)
% hold on 
% plot(k,0,'-o')

%generate a vector of positions of data symbols until which a corresponding
%data symbol with index f is used
for k=1:n_tsymbols
    if(k==n_tsymbols)
        u(k)=4+3*k;%The last training symbol
    else
        u(k)=1+3*k;
    end
end
f=1;

%decode the data symbols
for i=1:n_symbols
    if(i>=u(f))
        f=f+1; %f is the training symbol used for decoding at each time
    end
    y3=[];
    a=n_prefix+(n_prefix+symbol_size)*(i-1)+1;
    b=(n_prefix+symbol_size)*i;
    y0 = y(a:b);
    y1 = fft(y0)/sqrt(symbol_size);
    y2 = y1(n_lowf+2:n_lowf+n_data_symble+1);
    y3=y2./gain(:,f).*exp(-j*phase(:,f))./rand_realizations/gamma;
    y3_lowqam=y3(1:n_lowqam)*(log2(M)-1);
    y3_prime=y3(n_lowqam+1:n_data_symble-n_highqam)*log2(M);
    y3_highqam=y3(n_data_symble-n_highqam+1:end)*(log2(M)-1);

    z_lowqam=qamdemod(y3_lowqam,2^(log2(M)-1));
    for k=1:n_lowqam
        y4_lowqam((log2(M)-1)*(k-1)+1:(log2(M)-1)*k,1)=(fliplr(de2bi(z_lowqam(k),log2(M)-1)))';
    end
    outbits=[outbits;y4_lowqam];

     z_prime=qamdemod(y3_prime,M);
    for k=1:n_data_symble-n_highqam-n_lowqam
        y4_prime(log2(M)*(k-1)+1:(log2(M))*k,1)=(fliplr(de2bi(z_prime(k),log2(M))))';
    end
    outbits=[outbits;y4_prime];

    z_highqam=qamdemod(y3_highqam,2^(log2(M)-1));
    for k=1:n_highqam
        y4_highqam((log2(M)-1)*(k-1)+1:(log2(M)-1)*k,1)=(fliplr(de2bi(z_highqam(k),log2(M)-1)))';
    end
    outbits=[outbits;y4_highqam];
end


% figure(2)
correct=sum(outbits==bits);
d=find(outbits~=bits);
t=length(d);
wrong=zeros(t,1);
wrong(d)=1;
plot(wrong);

clear
numbits=2e5;
bits=randi([0 1],numbits,1);
M=32; % modulation index
n_data_symble=4000; %number of samples of data per symbol
n_lowf=112; %number of samples in the low frequency band
n_highf=908;%number of samples in the low frequency band
n_prefix=120;%number of samples in the cyclic prefix
symbol_size = (1+2*(n_highf+n_lowf+n_data_symble));
n_symbols = numbits/log2(M)/n_data_symble;
block_size= numbits/n_symbols;
package_size=block_size/log2(M);
P = 0.00125;
gamma = P*34;


% t=[0:1/4000:1-1/4000];
 t=rand(n_data_symble,1);
 rand_realizations =ones(n_data_symble,1).*exp(j*t*2*pi);

%rand_realizations= ones(4e3,1);

% rand_realizations=[1:4e3]';

x2_train=[zeros(n_lowf,1);rand_realizations;zeros(n_highf,1)];
x3_train=[0;x2_train;conj(fliplr(x2_train')')];
x4_train=ifft(x3_train)*sqrt(length(x3_train));
x5_train=[x4_train(end-n_prefix+1:end);x4_train]*gamma;
x=[];
x0_p=[];
for i=1:n_symbols
    if(i==3 ||i==6 || i==9)
        x=[x; x5_train];
    end
    x0=bits((i-1)*block_size+1:(i-1)*block_size+block_size);   %*sqrt(2);
    for k=1:package_size
        x0_prime(k,1)=bi2de(fliplr(x0((k-1)*log2(M)+1:k*log2(M))'));
    end
    x1=qammod(x0_prime,M)*gamma;
    x1_rand=x1.*rand_realizations/5;
%     x1=x1;
%     x1=x0.*rand_realizations;
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
load('IR0.mat', 'impulse');
h = impulse';
y = conv(x, impulse);
% no noise
y = y(1:length(x));
% plus random silence
% create random L0 and L1 pauses from uniform distribution 
L0 = round(1000*rand(1,1));
L1 = round(1000*rand(1,1));
y = [zeros(L0,1) ; y ; zeros(L1,1)];
% add noise
%y = y + 0.0001*randn(length(y),1);

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
thresh = 2e-8;
tau = -1;
for k = 1:len_y
    if (abs(s(k))>thresh)
        start = k-tau;  
        break;
    end
end
L1=start;

stem(y)
%disp(start)

hold on
plot(start,0,'-o')
len_sent=length(x);%hardcode later!!
y = y(start:(start+len_sent-1));%removed silence

% remember to normalize by gamma

%training 1
y0_train(:,1) = y(2*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*3);
y1_train (:,1) = fft(y0_train(:,1))/sqrt(symbol_size);
y2_train (:,1) = y1_train(n_lowf+2:package_size+1+n_lowf,1)/gamma; % 113:41112

gain (:,1) = abs(y2_train(:,1)./rand_realizations);% vector channel gain
phase_y=angle(y2_train(:,1));
phase_r=angle(rand_realizations);
phase (:,1)= wrapToPi(phase_y-fliplr(phase_r)); %vector of channel phase shift


%training 2
y0_train(:,2) = y(6*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*7);
y1_train(:,2) = fft(y0_train(:,2))/sqrt(symbol_size);
y2_train(:,2) = y1_train(n_lowf+2:package_size+1+n_lowf,2)/gamma; % 113:41112

gain(:,2) = abs(y2_train(:,2)./rand_realizations);% vector channel gain
phase_y=angle(y2_train(:,2));
phase_r=angle(rand_realizations);
phase(:,2)= wrapToPi(phase_y-fliplr(phase_r)); %vector of channel phase shift


%training 3
y0_train(:,3) = y(10*(symbol_size+n_prefix)+1+n_prefix:(n_prefix+symbol_size)*11);
y1_train(:,3) = fft(y0_train(:,3))/sqrt(symbol_size);
y2_train(:,3) = y1_train(n_lowf+2:package_size+1+n_lowf,3)/gamma; % 113:41112

gain(:,3) = abs(y2_train(:,3)./rand_realizations);% vector channel gain
phase_y=angle(y2_train(:,3));
phase_r=angle(rand_realizations);
phase(:,3)= wrapToPi(phase_y-fliplr(phase_r)); %vector of channel phase shift

y=[y(1:2*(symbol_size+n_prefix));
    y((n_prefix+symbol_size)*3+1:6*(symbol_size+n_prefix));
    y((n_prefix+symbol_size)*7+1:10*(symbol_size+n_prefix));
    y((n_prefix+symbol_size)*11+1:end)];%removed training
outbits = [];

% load('phase_chan.mat','phase_chan');
% phase=phase_chan;
% stem(y)
% hold on 
% plot(k,0,'-o')



for i=1:n_symbols
    if(i<4)
        y3=[];
        a=n_prefix+(n_prefix+symbol_size)*(i-1)+1;
        b=(n_prefix+symbol_size)*i;
        y0 = y(a:b);
        y1 = fft(y0)/sqrt(symbol_size);
        y2 = y1(n_lowf+2:n_lowf+n_data_symble+1);
        y3=y2./gain(:,1).*exp(-j*phase(:,1))./rand_realizations*5/gamma;
        z=qamdemod(y3,M);
        for k=1:package_size
            y4(5*(k-1)+1:5*k,1)=(fliplr(de2bi(z(k),log2(M))))';
        end
        outbits=[outbits;y4];
    elseif(i<7)
        y3=[];
        a=n_prefix+(n_prefix+symbol_size)*(i-1)+1;
        b=(n_prefix+symbol_size)*i;
        y0 = y(a:b);
        y1 = fft(y0)/sqrt(symbol_size);
        y2 = y1(n_lowf+2:n_lowf+n_data_symble+1);
        y3=y2./gain(:,2).*exp(-j*phase(:,2))./rand_realizations*5/gamma;
        z=qamdemod(y3,M);
        for k=1:package_size
            y4(5*(k-1)+1:5*k,1)=(fliplr(de2bi(z(k),log2(M))))';
        end
        outbits=[outbits;y4]; 
    else
        y3=[];
        a=n_prefix+(n_prefix+symbol_size)*(i-1)+1;
        b=(n_prefix+symbol_size)*i;
        y0 = y(a:b);
        y1 = fft(y0)/sqrt(symbol_size);
        y2 = y1(n_lowf+2:n_lowf+n_data_symble+1);
        y3=y2./gain(:,3).*exp(-j*phase(:,3))./rand_realizations*5/gamma;
        z=qamdemod(y3,M);
        for k=1:package_size
            y4(5*(k-1)+1:5*k,1)=(fliplr(de2bi(z(k),log2(M))))';
        end
        outbits=[outbits;y4];
        end

end
figure(2)
correct=sum(outbits==bits);
d=find(outbits~=bits);
t=length(d);
wrong=zeros(t,1);
wrong(d)=1;
plot(wrong);

% ind_wrong = [];
% for i = 1:length(bits)
%     if outbits(i) ~= bits(i) 
%         ind_wrong = [ind_wrong ; 1];
%     else
%         ind_wrong = [ind_wrong; 0];
%     end 
% end
% figure
% plot(ind_wrong(1:4000))

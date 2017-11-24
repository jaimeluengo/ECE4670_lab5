function outbits = dec(y,L1)
P = 0.00125*4820/4000;
n_plus = 120;
n_low=112;
n_high=708;
L = 9641;
num_bits = 2e5;

m = 3;
t = 2;

% estimate c threshold 
c = y(1:11)' * y(1:11) + y(end - 10: end)' * y(end - 10: end) * 30 / 20;

% for i=1:length(y)
%     if(abs(y(i))>0.005);
%         break
%     end
% end
% L1=i;
% % estimate L1
s = 0;
k = 0;
while s <= c
    s = 0;
    k = k + 1;
    for i = 1:m
        if k+1-i < 1
            break;
        end
        s = s+(y(k+1-i))^2;
    end
    s = s/m;
end
L1 = k-t;


y0 = y(L1+n_plus+1: L1+n_plus+L);
y0 = fft(y0)/sqrt(L);
y1 = y0(113:4113);

gain = abs(y1/sqrt(2*P));
bits = [];

counter = 0;
i=1;
while (counter < num_bits)
    y0 = y(L1+n_plus+(n_plus+L)*i+1: L1+n_plus+L+(n_plus+L)*i);
    y0 = fft(y0)/sqrt(L);
    y1 = y0(n_low+2:n_low+4001);
    z = conj(y1).*y1;
    for k= 1:4000
        if abs(y1(k)) >= sqrt(2*P)/2*gain(k);
        %if z(k) >= sqrt(2*P)/2*gain(k);
            bits = [bits ; 1];
        else
            bits = [bits; 0];
        end
        counter = counter+1;
        if counter >= num_bits
            break;
        end
    end
    i= i+1;
end
outbits = bits;
return

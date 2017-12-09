function outbits = sdec()
    y = audioread('rx.wav');
    P = 0.00125*4820/4000;
    n_plus = 120;
    n_low = 112;
    L = 9641;
    num_bits = 2e5;

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
    L1 = start;


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
            if abs(y1(k)) >= sqrt(2*P)/2*gain(k)
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

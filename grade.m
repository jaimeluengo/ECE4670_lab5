clear all
perf_ms = [];
perf_m = [];
rounds = 1;
for i = 1:rounds
    standard = 0;

    numbits=2e5; %length of bitsteam to send
    bits=randi([0 1],numbits,1); %bitstream

    if (standard == 1)
        x = senc(bits);
    else
        x = enc(bits);
    end

    % create random L0 and L1 pauses from uniform distribution 
    a = 0.25;
    b = 3;
    L0 = (b-a).*rand(1,1) + a;
    L1 = (b-a).*rand(1,1) + a;
    % change these lines before submission!!!!!
    cmd1 = sprintf('~/Desktop/ECE4670/ccplay/ccplay --prepause ');
    cmd2 = sprintf('%f --postpause %f --channel audio? tx.wav rx.wav', L0, L1);
    cmd = [cmd1 cmd2];
    system(cmd);

    if (standard == 1)
        outbits = sdec();
    else
        outbits = dec();
    end


    % Count correct ones and plot the errors
    correct = sum(outbits==bits);
    % d=find(outbits~=bits);
    % t=length(d);
    % wrong=zeros(t,1);
    % wrong(d)=1;
    % figure(1)
    % plot(wrong)

    avg_power = (x'*x)/length(x);
    % FIgure of merit
    Ne=numbits-correct;
    Pr=max(1,800*avg_power);
    data_rate=2e5*44.1e3/length(x);
    if (standard == 1)
        perf_s = (1-Ne/1e5)^10;
        perf_ms = [perf_ms ; perf_s];
    else
        perf = (data_rate*(1-Ne/1e5)^10)/Pr;
        perf_m = [perf_m ; perf];
    end 
end

if (standard == 1)
    mean_perf = mean(perf_ms);
else 
    mean_perf = mean(perf_m);
end

    
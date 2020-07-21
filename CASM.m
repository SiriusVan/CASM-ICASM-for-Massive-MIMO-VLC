function Establishment_of_Constellation()
clear all;   clc;
Ot=[0.116 0.116 0.5];%centre of the LED plane
Or=[0.116 0.116 1.8];%centre of the PD plane
v=0.058;
A=0.01^2;
semiangle=10;
k1=-log(2)/log(cosd(semiangle/2));
Xt=[Ot(1)-1.5*v Ot(1)-0.5*v Ot(1)+0.5*v Ot(1)+1.5*v; Ot(1)-1.5*v Ot(1)-0.5*v Ot(1)+0.5*v Ot(1)+1.5*v; Ot(1)-1.5*v Ot(1)-0.5*v Ot(1)+0.5*v Ot(1)+1.5*v; Ot(1)-1.5*v Ot(1)-0.5*v Ot(1)+0.5*v Ot(1)+1.5*v];
Yt=Xt.';
Xr=[Or(1)-1.5*v Or(1)-0.5*v Or(1)+0.5*v Or(1)+1.5*v; Or(1)-1.5*v Or(1)-0.5*v Or(1)+0.5*v Or(1)+1.5*v; Or(1)-1.5*v Or(1)-0.5*v Or(1)+0.5*v Or(1)+1.5*v; Or(1)-1.5*v Or(1)-0.5*v Or(1)+0.5*v Or(1)+1.5*v];
Yr=[Or(2)-1.5*v Or(2)-0.5*v Or(2)+0.5*v Or(2)+1.5*v; Or(2)-1.5*v Or(2)-0.5*v Or(2)+0.5*v Or(2)+1.5*v; Or(2)-1.5*v Or(2)-0.5*v Or(2)+0.5*v Or(2)+1.5*v; Or(2)-1.5*v Or(2)-0.5*v Or(2)+0.5*v Or(2)+1.5*v].';
H=zeros(16,16);
for i=1:16%receiving PD index
    for j=1:16%transmitting LED index
        horizontal=sqrt((Xr(i)-Xt(j))^2+(Yr(i)-Yt(j))^2);
        if horizontal/(Or(3)-Ot(3))<tand(semiangle)
            angle=atand(horizontal/(Or(3)-Ot(3)));
            H(i,j)=cosd(angle)^(k1+1)*((k1+1)*A)/(2*pi*(Or(3)-Ot(3))^2);
        end
    end
end
%%%%%%%Now we start to choose the constellation
k1=2;%Number of LEDs in a constellation
L1=nchoosek(16,k1);
c=1:1:16;
C=nchoosek(c,k1);
LEDICIMAX=zeros(1,L1);
L1
for m1=1:L1
    LEDCOMBINATION=nchoosek(C(m1,:),2);
    LEDC=nchoosek(k1,2);
    for i=1:LEDC
        for j=1:16
        if inf>H(j,LEDCOMBINATION(i,1))/H(j,LEDCOMBINATION(i,2))>LEDICIMAX(m1)
            LEDICIMAX(m1)=H(j,LEDCOMBINATION(i,1))/H(j,LEDCOMBINATION(i,2));
        end
        if inf>H(j,LEDCOMBINATION(i,2))/H(j,LEDCOMBINATION(i,1))>LEDICIMAX(m1)
            LEDICIMAX(m1)=H(j,LEDCOMBINATION(i,2))/H(j,LEDCOMBINATION(i,1));
        end
        end
    end
end
LEDICIMAX=LEDICIMAX.';
Constelations=[LEDICIMAX,C];
FINALCHOICE=[];
m1=0;
threshold1=2;
for i=1:L1
if Constelations(i,1)<threshold1
    m1=m1+1;
    FINALCHOICE(m1,:)=Constelations(i,:);
end
end
x=[];
for j = 1: m1-1
  for i= 1:m1-j
    if FINALCHOICE(i,1)>FINALCHOICE(i+1,1)
    x=FINALCHOICE(i,:);
    FINALCHOICE(i,:) = FINALCHOICE(i+1,:);
    FINALCHOICE(i+1,:) = x;
    end
  end
end
m1=2^(floor(log2(m1)));
m1
FINALCHOICE=FINALCHOICE(1:m1,:);
FINALLED=FINALCHOICE(:,2:k1+1);
%%%%Now we begin to simulate the Modulation part
Mqam1 = 4;   %QAM Modulation Order
bit_SMsym = log2(Mqam1^k1*m1);      % number of bit per spatial modulation sysmbol
Symbols1=3000;
Nbits = bit_SMsym*Symbols1;  % Number of bits to be simulated.
hmodem =modem.qammod('M',Mqam1,  'SymbolOrder', 'Gray','InputType', 'bit');%QAM Modulation
hdemodem =modem.qamdemod('M', Mqam1,'SymbolOrder','Gray','OutputType','bit');%QAM Demodulation
bit_T = zeros(Nbits,1);%bit_T = randi([0 1],Nbits,1);%Transmitted bits
Nt    = log2(m1); %Bits transmitted by the combination of LEDS
Nobit = k1*log2(Mqam1); %The bits transmitted by QAM
x = reshape(bit_T,Nobit+Nt,[]);%
ant_no =  bi2de([x(1:Nt,:)].' , 'left-msb') + 1;
for i=1:k1
digMod(i,:) = modulate(hmodem,x(Nt+1+(i-1)*log2(Mqam1):Nt+i*log2(Mqam1),:));
end
signal = zeros(16 , Symbols1);
for j=1:k1 
for i = 1 : Symbols1
signal(FINALLED(ant_no(i),j),i)=digMod(j,i);
end
end
%%%%%%%%%Now we begin to simulate the Transmission Part
HS=H*signal;
SNR = 40 : 2 :80;
Eac = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No= (Eac)*10.^(-SNR/10);    % noise variance SNR=10lg£¨E/N£©No= (Eac)*10.^(-SNR/10);  
L_SNR=length(SNR);
ber1= zeros (L_SNR,1);
bit_R=zeros(Nbits ,1);%Received bits
for ii=1:L_SNR
    for j = 1 : size(HS ,2)
    noise = sqrt(.5)*(randn(16, 1) + 1i*randn(16 , 1))* sqrt(No(ii));
    Reception=HS(:,j)+noise;%Reception=HS(:,j)+noise;
    MEAN=zeros(2,m1);
      for jj=1:m1
         Aberration=zeros(1,2^(Nobit));
        for jjj=0:(2^(Nobit)-1)
           nobit=de2bi(jjj,Nobit,'left-msb');
           nobit=reshape(nobit,log2(Mqam1),[]);
           Mod = modulate(hmodem,nobit(1:log2(Mqam1),:));
           TESTsignal=zeros(1,16);
           for i = 1 : k1
           TESTsignal(FINALLED(jj,i))=Mod(i);
           end
           TESTreception=H*TESTsignal.';
           Aberration(jjj+1)=mean((Reception-TESTreception).^2);
        end 
          [MEAN(2,jj),MEAN(1,jj)] = min(Aberration);
      end
    [~,jj]=min(MEAN(2,:));
    jjj=MEAN(1,jj);
    bit_R(1+(j-1)*bit_SMsym:j*bit_SMsym)=[de2bi(jj-1,Nt,'left-msb') de2bi(jjj-1,Nobit,'left-msb')];
    end
    [~,ber1(ii,1)] = biterr(bit_T,bit_R);
end
%%%%%%%Now we start to choose the constellation
k2=2;%Number of LEDs in a constellation
L2=nchoosek(16,k2);
c=1:1:16;
C=nchoosek(c,k2);
LEDICIMAX=zeros(1,L2);
L2
for m2=1:L2
    LEDCOMBINATION=nchoosek(C(m2,:),2);
    LEDC=nchoosek(k2,2);
    for i=1:LEDC
        for j=1:16
        if inf>H(j,LEDCOMBINATION(i,1))/H(j,LEDCOMBINATION(i,2))>LEDICIMAX(m2)
            LEDICIMAX(m2)=H(j,LEDCOMBINATION(i,1))/H(j,LEDCOMBINATION(i,2));
        end
        if inf>H(j,LEDCOMBINATION(i,2))/H(j,LEDCOMBINATION(i,1))>LEDICIMAX(m2)
            LEDICIMAX(m2)=H(j,LEDCOMBINATION(i,2))/H(j,LEDCOMBINATION(i,1));
        end
        end
    end
end
LEDICIMAX=LEDICIMAX.';
Constelations=[LEDICIMAX,C];
FINALCHOICE=[];
m2=0;
threshold2=2;
for i=1:L2
if Constelations(i,1)<threshold2
    m2=m2+1;
    FINALCHOICE(m2,:)=Constelations(i,:);
end
end
x=[];
for j = 1: m2-1
  for i= 1:m2-j
    if FINALCHOICE(i,1)<FINALCHOICE(i+1,1)
    x=FINALCHOICE(i,:);
    FINALCHOICE(i,:) = FINALCHOICE(i+1,:);
    FINALCHOICE(i+1,:) = x;
    end
  end
end
m2=2^(floor(log2(m2)));
m2
FINALCHOICE=FINALCHOICE(1:m2,:);
FINALLED=FINALCHOICE(:,2:k2+1);
%%%%Now we begin to simulate the Modulation part
Mqam2 = 4;   %QAM Modulation Order
bit_SMsym = log2(Mqam2^k2*m2);      % number of bit per spatial modulation sysmbol
Symbols2=3000;
Nbits = bit_SMsym*Symbols2;  % Number of bits to be simulated.
hmodem =modem.qammod('M',Mqam2,  'SymbolOrder', 'Gray','InputType', 'bit');%QAM Modulation
hdemodem =modem.qamdemod('M', Mqam2,'SymbolOrder','Gray','OutputType','bit');%QAM Demodulation
bit_T = zeros(Nbits,1);%bit_T = randi([0 1],Nbits,1);%Transmitted bits
Nt    = log2(m2); %Bits transmitted by the combination of LEDS
Nobit = k2*log2(Mqam2); %The bits transmitted by QAM
x = reshape(bit_T,Nobit+Nt,[]);%
ant_no =  bi2de([x(1:Nt,:)].' , 'left-msb') + 1;
Calculation=[1:1:m2; zeros(1,m2)]; 
for i=1:Symbols2
   Calculation(2,ant_no(i))=Calculation(2,ant_no(i))+1;
end
for i=1:m2-1
    for j=1:m2-i
        if Calculation(2,j)<Calculation(2,j+1)
            xyz=Calculation(:,j+1);
            Calculation(:,j+1)=Calculation(:,j);
            Calculation(:,j)=xyz;
        end
    end
end
Transfer=zeros(m2,k2+1);
for i=1:m2
     Transfer(Calculation(1,i),:)=FINALCHOICE(i,:);
end
FINALLED=Transfer(:,2:k2+1);
for i=1:k2
digMod(i,:) = modulate(hmodem,x(Nt+1+(i-1)*log2(Mqam2):Nt+i*log2(Mqam2),:));
end
signal = zeros(16 , Symbols2);
for j=1:k2 
for i = 1 : Symbols2
signal(FINALLED(ant_no(i),j),i)=digMod(j,i);
end
end
%%%%%%%%%Now we begin to simulate the Transmission Part
HS=H*signal;
SNR = 40 : 2 :80;
Eac = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No= (Eac)*10.^(-SNR/10);    % noise variance SNR=10lg£¨E/N£©No= (Eac)*10.^(-SNR/10);  
L_SNR=length(SNR);
ber2= zeros (L_SNR,1);
bit_R=zeros(Nbits ,1);%Received bits
for ii=1:L_SNR
    for j = 1 : size(HS ,2)
    noise = sqrt(.5)*(randn(16, 1) + 1i*randn(16 , 1))* sqrt(No(ii));
    Reception=HS(:,j)+noise;%Reception=HS(:,j)+noise;
    MEAN=zeros(2,m2);
      for jj=1:m2
         Aberration=zeros(1,2^(Nobit));
        for jjj=0:(2^(Nobit)-1)
           nobit=de2bi(jjj,Nobit,'left-msb');
           nobit=reshape(nobit,log2(Mqam2),[]);
           Mod = modulate(hmodem,nobit(1:log2(Mqam2),:));
           TESTsignal=zeros(1,16);
           for i = 1 : k2
           TESTsignal(FINALLED(jj,i))=Mod(i);
           end
           TESTreception=H*TESTsignal.';
           Aberration(jjj+1)=mean((Reception-TESTreception).^2);
        end 
          [MEAN(2,jj),MEAN(1,jj)] = min(Aberration);
      end
    [~,jj]=min(MEAN(2,:));
    jjj=MEAN(1,jj);
    bit_R(1+(j-1)*bit_SMsym:j*bit_SMsym)=[de2bi(jj-1,Nt,'left-msb') de2bi(jjj-1,Nobit,'left-msb')];
    end
    [~,ber2(ii,1)] = biterr(bit_T,bit_R);
end
%%%%Now we begin to draw the figures
semilogy(SNR,ber1(:,1),'-+',SNR,ber2(:,1),'r:*');
grid on;
xlabel('$$OSNR(dB)$$','Interpreter','latex')
ylabel('Average BER','Interpreter','latex')
title('BER Comparison Between VLC CASM and Improved CASM With Special Signals')%title('Bit Error Rate of VLC CASM Using ML Detector With Different Parameters')
legend(['CASM ',num2str(Mqam1),'-QAM k=',num2str(k1),' Threshold=',num2str(threshold1),' m=',num2str(m1)] , ['Improved CASM ',num2str(Mqam2),'-QAM k=',num2str(k2),' Threshold=',num2str(threshold2),' m=',num2str(m2)], 'Location','NorthEast')
end

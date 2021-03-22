clc;
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
L=nchoosek(16,k1);
c=1:1:16;
C=nchoosek(c,k1);
LEDICIMAX=zeros(1,L);
L
for m1=1:L
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
threshold1=1.2;
for i=1:L
if Constelations(i,1)<threshold1
    m1=m1+1;
    FINALCHOICE(m1,:)=Constelations(i,:);
end
end
x=[];
for j = 1: m1-1
  for i= 1:m1-1
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
Mpam1 = 16;   %PAM Modulation Order
bit_SMsym = log2(Mpam1^k1*m1);      % number of bit per spatial modulation sysmbol
Symbols1=1000;
Nbits = bit_SMsym*Symbols1;  % Number of bits to be simulated.
%hmodem =qammod('M',Mqam,  'SymbolOrder', 'Gray','InputType', 'bit');%QAM Modulation
%hdemodem =qamdemod('M', Mqam,'SymbolOrder','Gray','OutputType','bit');%QAM Demodulation

bit_T = randi([0 1],Nbits,1);%Transmitted bits
%bit_T = ones(Nbits,1);
Nt    = log2(m1); %Bits transmitted by the combination of LEDS
Nobit = k1*log2(Mpam1); %The bits transmitted by PAM
x = reshape(bit_T,Nobit+Nt,[]);%
ant_no =  bi2de(x(1:Nt,:).' , 'left-msb') + 1;
digMod=zeros(k1,size(x,2));
for i=1:k1
    for j=1+(i-1)*log2(Mpam1):i*log2(Mpam1)
        digMod(i,:)=digMod(i,:)+1+(2^(log2(Mpam1)-1-(j-(1+(i-1)*log2(Mpam1)))))*x(Nt+j,:);
    end
%digMod(i,:) = modulate(hmodem,x(Nt+1+(i-1)*log2(Mqam):Nt+i*log2(Mqam),:));
end
signal = zeros(16 , Symbols1);
for j=1:k1 
for i = 1 : Symbols1
signal(FINALLED(ant_no(i),j),i)=digMod(j,i);
end
end
%%%%%%%%%Now we begin to simulate the Transmission Part
HS=H*signal;
SNR = 50 : 1 :90;
%Eac = (mean(hmodem.Constellation .* hmodem.Constellation));
%Eac=10;
Eac = mean(mean(digMod.*digMod))/200;
No= (Eac)*10.^(-SNR/10);    % noise variance SNR=10lg?¡§E/N??No= (Eac)*10.^(-SNR/10);  
L_SNR=length(SNR);
ber1= zeros (L_SNR,1);
bit_R=zeros(Nbits ,1);%Received bits
%%%%%Optimization begins£¡
%MEAN=zeros(2,m);
COMPARISON=zeros(16,2^(Nobit),m1);
SAMPLE=zeros(k1,2^(Nobit)); %Nobit = k*log2(Mpam);
 for jjj=0:(2^(Nobit)-1)
           nobit=de2bi(jjj,Nobit,'left-msb');
           nobit=reshape(nobit,log2(Mpam1),[]);
           %Mod = modulate(hmodem,nobit(1:log2(Mpam),:));
             Mod=zeros(1,k1);
            for jjjj=1:log2(Mpam1)
             Mod=Mod+1+(2^(log2(Mpam1)-jjjj))*nobit(jjjj,:);
            end
            SAMPLE(:,jjj+1)=Mod.';
 end
 
 for jj=1:m1
     for jjj=0:(2^(Nobit)-1)
     TESTsignal=zeros(1,16);
           for i = 1 : k1
           TESTsignal(FINALLED(jj,i))=SAMPLE(i,jjj+1);
           end
     COMPARISON(:,jjj+1,jj)=TESTsignal.';
     end
     COMPARISON(:,:,jj)=H*COMPARISON(:,:,jj);
 end
 Aberration=zeros(1,2^(Nobit),m1);
 PAM=zeros(1,1,m1);
 
%%%%%%%%%%%Final Section


 for ii=1:L_SNR
    %noise = sqrt(.5)*(randn(16, 1) + 1i*randn(16 , 1))* sqrt(No(ii)); 
    noise = sqrt(.5)*(randn(16, Symbols1) )* sqrt(No(ii));
    Reception=HS+noise;
    %Reception=reshape(Reception,16,[]);
    for j = 1 : Symbols1
    RECEPTION=repmat(Reception(:,j),1,2^(Nobit),m1);
    Aberration=mean(((RECEPTION-COMPARISON).^2).^0.1);
    [~,PAM]=min(Aberration);
    [~,jj]=min(min(Aberration));
    jjj=PAM(1,1,jj);
    bit_R(1+(j-1)*bit_SMsym:j*bit_SMsym)=[de2bi(jj-1,Nt,'left-msb') de2bi(jjj-1,Nobit,'left-msb')];
    end   
    [~,ber1(ii,1)] = biterr(bit_T,bit_R);
 end
%%%%%Optimization ends£¡
semilogy(SNR,ber1(:,1),'color',[0,0.75,0.75],'linestyle','-','LineWidth',1);
grid on;
xlabel('$$OSNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('Bit Error Rate of VLC CASM using ML detector')
legend([num2str(Mpam1),'-PAM k=',num2str(k1),' Threshold=',num2str(threshold1),'Symbols=',num2str(Symbols1)] ,  'Location','NorthEast')

  %%%%%
   Reception=HS;
    %Reception=reshape(Reception,16,[]);
    for j = 1 : Symbols1
    RECEPTION=repmat(Reception(:,j),1,2^(Nobit),m1);
    Aberration=mean((RECEPTION-COMPARISON).^2);
    [ABCD,PAM]=min(Aberration);
    [ABCDE,jj]=min(min(Aberration));
    jjj=PAM(1,1,jj);
    bit_R(1+(j-1)*bit_SMsym:j*bit_SMsym)=[de2bi(jj-1,Nt,'left-msb') de2bi(jjj-1,Nobit,'left-msb')];
    end   
    [~,ber1] = biterr(bit_T,bit_R);
    ber1
%%%%%%%


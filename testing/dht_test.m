clear all

R = 10
jmodes = 128
ii = 215
N=256
n = ii
h=[]


[H,kk,rr,I,KK,RR]=dht([],R,jmodes,ii)

%save(['output/ker_' num2str(floor(ii))],'I','kk','rr','KK','RR')
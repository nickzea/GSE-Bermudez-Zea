function[Pcal, Qcal]=calcPQ(MV,G,B,teta,Nb)
Pcal=zeros(Nb,1);
Qcal=zeros(Nb,1);
for i=1:Nb
        for k=1:Nb
            Pcal(i)=(MV(i))*(MV(k))*(G(i,k)*cos(teta(i)-teta(k))+B(i,k)*sin(teta(i)-teta(k)))+Pcal(i);
            Qcal(i)=(MV(i))*(MV(k))*(G(i,k)*sin(teta(i)-teta(k))-B(i,k)*cos(teta(i)-teta(k)))+Qcal(i);
        end

end
end 
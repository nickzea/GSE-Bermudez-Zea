close all
clear
%Base de datos del sistema
BMva=100;
[bus, gen, branch]=datossistema();

% Datos de barras
Nb  = size(bus,1);           %Número de barras.
tipo  = bus(:,2);            %Tipo de barra.
MV    = ones(Nb,1);            %Magnitud de Tensión de barra.
V     = bus(:,10);
teta  = zeros(Nb,1);            %Angulo da tensión de barra.
          
% Datos de linea 
nlin = size(branch,1);          %Numero de Líneas
DE   = branch(:,1);             %Barra DE (k).
PARA = branch(:,2);             %Barra PARA (m).
r    = branch(:,3);
x    = branch(:,4);             %Reactancias en pu.
b    = branch(:,5);             %Susceptancia shunt de la línea.
a    = branch(:,9);             %Valor del tap de los transformadores.

slack = find(tipo == 3);   %% Barra slack
pv  = find(tipo == 2 );    %% índices de barras PV
npv = length(pv);
pq  = find(tipo == 1);     %% índices da barras PQ
npq = length(pq);
pvq = find(tipo==1 | tipo==2);
npvq = length(pvq);

%% Construccion de la Matriz Ybus


y=zeros(24,24);
B_y=zeros(24,24);
A= zeros(24,24);

for e=1:1:nlin
    A(DE(e),PARA(e))=a(e);
    if a(e)==0
        y(DE(e),PARA(e))=1/(r(e)+1j*x(e))+y(DE(e),PARA(e));   
        y(PARA(e),DE(e))= y(DE(e),PARA(e));
    else
        y(DE(e),PARA(e))=1/((r(e)+1j*x(e))*a(e));
        y(PARA(e),DE(e))= y(DE(e),PARA(e)); 
        if V(DE(e))>V(PARA(e))
        y(PARA(e),PARA(e,1))=((1-a(e))/a(e)^2)*(1/(r(e)+1j*x(e)))+ y(PARA(e),PARA(e));
        y(DE(e),DE(e))=((a(e)-1)/a(e))*(1/(r(e)+1j*x(e)))+y(DE(e),DE(e));
        end
        if V(DE(e))<V(PARA(e))
        y(DE(e),DE(e))=((1-a(e))/a(e)^2)*(1/(r(e)+1j*x(e)))+y(DE(e),DE(e));
        y(PARA(e),PARA(e))=((a(e)-1)/a(e,1))*(1/(r(e)+1j*x(e)))+ y(PARA(e),PARA(e));
       end
    end
     B_y(DE(e),DE(e))=b(e)/2+B_y(DE(e),DE(e));
     B_y(PARA(e),PARA(e))=b(e)/2+B_y(PARA(e),PARA(e));
end

Y=-y+diag(diag(y))+diag(sum(y,2))+1j*B_y;

B=imag(Y);

%% Potencias y Voltaje
PV=zeros(2,7);
QV=zeros(2,7);
for z=1:1:7

Pg=zeros(Nb,1);
Pd=zeros(Nb,1);
Qg=zeros(Nb,1);
Qd=zeros(Nb,1);

for e=1:1:33
    j=gen(e,1);
    Pg(j)=Pg(j)+gen(e,2);
    Qg(j)=Qg(j)+gen(e,3);
    

    MV(j)=gen(e,6);
    
    if e<=length(bus)
    k=bus(e,1);
    Pd(k)=Pd(k)+bus(e,3);
    Qd(k)=Qd(k)+bus(e,4);
    end
end 

Pd=Pd*(1+(z-1)/10);
Qd=Qd*(1+(z-1)/10);
Pesp = (Pg - Pd)/100;
Qesp = (Qg - Qd)/100;



%% Newton Raphson

% Inicializar contador de iteraciones y tolerancia
tol=1;
it = 0; 

%Fast decoupled (XX)
%G=zeros(Nb,Nb);

%Fast decoupled (BB)
%G=real(Y);

tic
while tol> 1e-5

slack = find(tipo == 3);   %% Barra slack
pv  = find(tipo == 2 );    %% índices de barras PV
npv = length(pv);
pq  = find(tipo == 1);     %% índices da barras PQ
npq = length(pq);
pvq = find(tipo==1 | tipo==2);
npvq = length(pvq);
cont=0;

%Fast decoupled (XB)
G=zeros(Nb,Nb);
%G=real(Y);

[~, Qcal1]=calcPQ(MV,G,B,teta,Nb);

    % Crear Jacobiano
    % Matriz H : Derivada de potencia activa respecto al ángulo
    HJ = zeros(npvq,npvq);
    for k=1:npvq
        for m=1:npvq
            if pvq(k)==pvq(m)
                HJ(k,m)=-B(pvq(k),pvq(k))*MV(pvq(k))^2-Qcal1(pvq(k));
                 
            else
                HJ(k,m)=MV(pvq(k))*MV(pvq(m))*(G(pvq(k),pvq(m))*sin(teta(pvq(k))-teta(pvq(m)))-B(pvq(k),pvq(m))*cos(teta(pvq(k))-teta(pvq(m)))); 
            end
        end
    end
 

%Fast decoupled (BX) 
%G=zeros(Nb,Nb);
G= real(Y);

[~, Qcal2]=calcPQ(MV,G,B,teta,Nb);

        % Matriz L : Derivada de potencia reactiva respecto a la tensión
    LJ = zeros(npq,npq);
    for k=1:npq
        for m=1:npq
            if pq(k)==pq(m)
                LJ(k,m)=(1/MV(pq(k)))*Qcal2(pq(k))-B(pq(k),pq(k))*MV(pq(k));
            else
                LJ(k,m)=MV(pq(k))*(G(pq(k),pq(m))*sin(teta(pq(k))-teta(pq(m)))-B(pq(k),pq(m))*cos(teta(pq(k))-teta(pq(m))));              
            end
        end
    end 

G= real(Y);
[Pcal, Qcal]=calcPQ(MV,G,B,teta,Nb);
   
    %% flujo AC
    fxang=teta(pvq);
    fxvol=MV(pq);
    delp=Pesp(pvq)-Pcal(pvq);
    delq=Qesp(pq)-Qcal(pq);
    
    delang=inv(HJ)*delp;
    delvol=inv(LJ)*delq;
    fxvol=fxvol+delvol;
    fxang=fxang+delang;
    teta(pvq)=fxang;
    MV(pq)=fxvol; 
    
 %% Error
tol=max(abs([delp;delq]));
it=it+1;
end
toc

%% Resultados
%Teta en grados
%teta=teta*180/pi

%Calcula la potencia reactiva y activa para cada nodo con los V y theta
%obtenidos
[Pcal, Qcal]=calcPQ(MV,G,B,teta,Nb);

%Fasor V
for i=1:Nb
V(i)=MV(i)*exp(1j*teta(i));
end


%Perdidas de línea y Potencias injectadas
S=zeros(Nb,Nb);
Lij=zeros(nlin,1);
for m=1:1:nlin
    i=DE(m);
    j=PARA(m);
    
    S(i,j)=(V(i)*(conj(V(j))-conj(V(i)))*conj(Y(i,j)))*BMva;
    S(j,i)=(V(j)*conj((V(i)-V(j))*Y(j,i)))*BMva;
    
    Lij(m,1)=(S(i,j)+S(j,i));

end

%P Q generadas
Qg=(Qcal)*BMva+Qd;
Pg=(Pcal)*BMva+Pd;



%Tabla de resultados

fprintf('Converge en '); fprintf('%2.0f', it); fprintf(' iteraciones con una tolerancia de '); fprintf('%s.', tol);
disp('                                                                                       ');
disp('-----------------------------------------------------------------------------------------');
disp('                              Newton Raphson Loadflow Analysis ');
disp('-----------------------------------------------------------------------------------------');
disp('| Bus |    V   |  Angle  |        Load        |     Generation     |');
disp('| No  |   pu   |  Degree |    MW   |   MVar   |    MW   |  Mvar    |');
for m = 1:Nb
    fprintf('%3g', bus(m,1)); fprintf('  %8.4f', MV(m)); fprintf('   %8.4f', teta(m));
    fprintf('  %8.3f', Pd(m)); fprintf('   %8.3f', Qd(m)); 
    fprintf('  %8.3f', Pg(m)); fprintf('   %8.3f', Qg(m)); 
    disp('                                                                                       ');
end

disp('---------------------------------------------------------------------');
fprintf('Total:'); fprintf('                    %8.3f',sum(Pd)); fprintf('   %8.3f', sum(Qd)); fprintf('  %8.3f', sum(Pg)); fprintf('   %8.3f', sum(Qg));

disp('                                                                                       ');
fprintf('Pérdidas P(MW):'); fprintf('%8.3f',(sum(Pg)-sum(Pd)));
disp('                                                                                       ');
disp('-------------------------------------------------------------------------------------------------------------------');
disp('                                                   Branch Data ');
disp('-------------------------------------------------------------------------------------------------------------------');
disp('| Brnch|   From |  To     |       From Bus Injection     |     To Bus Injection    |      Pérdidas de línea      |  ');
disp('| No   |   Bus  |  Bus    |    MW         |   MVar       |    MW      |  Mvar      |      MW      |    Mvar      | ');

P_inj=zeros(nlin,1);
Q_inj=zeros(nlin,1);
P_jni=zeros(nlin,1);
Q_jni=zeros(nlin,1);
for m = 1:nlin
    i=DE(m);
    j=PARA(m);
    
    P_inj(m)=real(S(i,j));
    Q_inj(m)=imag(S(i,j));
    P_jni(m)=real(S(j,i));
    Q_jni(m)=imag(S(j,i));
    
    fprintf('%3g', m); fprintf('  %8.0f', i); fprintf('   %8.0f', j);
    fprintf('    %8.3f', P_inj(m)); fprintf('      %8.3f',  Q_inj(m)); 
    fprintf('       %8.3f', P_jni(m)); fprintf('       %8.3f', Q_jni(m)); 
    fprintf('       %8.3f', real(Lij(m))); fprintf('       %8.3f',imag(Lij(m))); 
    disp('                                                                                       ');

end
disp('---------------------------------------------------------------------------------------------------------------');
fprintf('Total:'); fprintf('                                                                                 %8.3f',sum(real(Lij))); fprintf('       %8.3f',sum(imag(Lij)));
disp('                                                                                       ');

PV(1,z)=sum(Pd);
PV(2,z)=MV(10);

QV(1,z)=sum(Qg);
QV(2,z)=MV(10);

end

subplot(2,1,1)
plot(PV(1,:),PV(2,:),'o');
grid on
title('\fontsize{16}Curva de Estabilidad de Voltaje Bus 10')
ylabel('|V| p.u.');
xlabel('Potencia Activa Carga del Sistema MW'); 

subplot(2,1,2)
plot(QV(1,:),QV(2,:),'o');
grid on
title('\fontsize{16}Curva de Estabilidad de Voltaje Bus 10')
ylabel('|V| p.u.');
xlabel('Potencia Reactiva Carga del Sistema MVar'); 

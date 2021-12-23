clf;
clear all;
clc;
%%Extended Kalman Filter for linear Pendulum
%Initial Condition
g=10
l=1
H0=[1 0]
R=[0.05]
Pminus=l
I=[1 0
    0 1]
Xprediction1step=[0
    0]
Q=[0 0
    0 0.01]
Z(1,1)=0.02
iterasi=400
%%Generating Zero Mean White Noise for X2 state
for i=1:+1:iterasi
    W(1,i)=sqrt(0.05).*randn(1)+0
end
%%Generating Zero Mean White Noise for Measurement
for i=1:+1:iterasi
    V(1,i)=sqrt(0.01).*randn(1)+0
end
%%Iteration
for k=1:+1:iterasi
    H1=H0.' %Fungsi transpose matriks H
    K1=inv(H0*Pminus*H1+R) %Fungsi invers dari matriks
    K=Pminus*H1*K1 %Kalman Gain
    P=(I-K*H0)*Pminus %Fungsi Pk
    Xestimate=Xprediction1step+K*(Z(1,k)-H0*Xprediction1step) %Fungsi X estimasi
    Z(1,k+1)=[1 0]*Xestimate+V(1,k) %Fungsi hasil pengukuran
    X1prediction1step=Xestimate(1,1)+Xestimate(2,1) %Fungsi X1 prediksi untuk 1 step kedepan
    X2prediction1step=(g/l)*sin(Xestimate(1,1))+(1/l)*W(1,k) %Fungsi X2 prediksi untuk 1 step kedepan
    Xprediction1step=[X1prediction1step %Fungsi X
        X2prediction1step]
    %%Menyimpan X1 dan X2 untuk dilakukan ploting nantinya
    X1(1,k)=Xestimate(1,1)
    X2(1,k)=Xestimate(2,1)
    X3(1,k)=X1prediction1step
    X4(1,k)=X2prediction1step
    %%Fungsi untuk mencari matriks jacobian X
    jacobianFx=[1 1
        (g/l)*cos(Xestimate(1,1)) 0]
    jacobianFxTranspose=jacobianFx.'
    Pprediction=jacobianFx*P*jacobianFxTranspose+Q %Update parameter P prediksi
    Pminus=Pprediction %Update parameter P berdasarkan k sebelumnya untuk iterasi selanjutnya
end
figure(1)
plot(X1)
hold on
plot(X3);plot(Z);
hold off
title('X1');legend('X1(k|k)','X1(k|k-1)','Hasil Pengukuran Z')
figure(2)
plot(X2)
hold on
plot(X4)
hold off
title('X2');legend('X2(k|k)','X2(k|k-1)')
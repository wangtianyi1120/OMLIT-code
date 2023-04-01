% Reflection and transmission coefficient calculator of lossy multilayer 

%This program produces the complex reflection and transmission coefficients
%of a multilayer stack given the angle of incidence, polarization,
%wavelength, complex-refractive index of each layer, and thickness of each
%layer. The program assumes lossless dielectric incident and exit media. 
%The program also assumes that all layers are non-magnetic. Magnetic media
%can be handled by more general theory (see technical report cited below). 

%This program was inspired by the technical report written by K. Pascoe, 
%"Reflectivity and Transmissivity through Layered, Lossy Media: A
%User-Friendly Approach," 2001. See the following link for the technical
%paper:
%http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix=html&identifier=ADA389099

function [r,t,R,T,A] = multilayer_model(l,d,n,t0,polarization)
%l = free space wavelength, nm
%d = layer thickness vector, nm
%n = layer complex refractive index vector �㸴������ʸ��
%t0= angle of incidence
%polarization should be 0 for TE (s-polarized), otherwise TM (p-polarized)

Z0=376.730313; %impedance of free space, Ohms ���ɿռ��迹

%the line below had mistakenly been a Z0/n instead of a n/Z0 in version 1!
Y=n./Z0; %admittance in terms of impedance of free space and refractive index, assuming non-magnetic media
g=1i*2*pi*n/l; %propagation constant in terms of free space wavelength and refractive index

%all of the calculations rely on cosine of the complex angle, but we can
%only find the sine of the complex angle from snells law. So we use the
%fact that cos(asin(x))=sqrt(1-x^2)
%t=asin(n(1)./n*sin(t0)), complex theta for each layer
ct=sqrt(1-(n(1)./n*sin(t0)).^2); %cosine theta

if polarization==0
    eta=Y.*ct; %tilted admittance, TE case
else
    eta=Y./ct; %tilted admittance, TM case
end

delta=1i*g.*d.*ct;

M=zeros(2,2,length(d));

for j=1:length(d)
    M(1,1,j)=cos(delta(j));
    M(1,2,j)=1i./eta(j).*sin(delta(j));
    M(2,1,j)=1i*eta(j).*sin(delta(j));
    M(2,2,j)=cos(delta(j));
end    

M_t=[1,0;0,1]; %M total
for j=2:(length(d)-1)
    M_t=M_t*M(:,:,j);
end

r=(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))-(M_t(2,1)+M_t(2,2)*eta(end)))/(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)));
t=2*eta(1)/(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)));

R=abs(r)^2;
T=real(eta(end)/eta(1))*abs(t)^2;
A=(4*eta(1)*real((M_t(1,1)+M_t(1,2)*eta(end))*conj(M_t(2,1)+M_t(2,2)*eta(end))-eta(end)))/abs(eta(1)*(M_t(1,1)+M_t(1,2)*eta(end))+(M_t(2,1)+M_t(2,2)*eta(end)))^2;





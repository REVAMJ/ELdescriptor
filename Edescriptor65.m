function GDdescriptor = Edescriptor65(patches)
numPatches=size(patches,2);
%weighting functions parameters
sig1=3.;
sig2=5.;
sig3=8.75;
%sig1=3;
%sig2=5.5;
%sig3=9.75;
r1=14.5;
r2=31.5;
Areas=17;
binNum=8;
GDdescriptor=zeros(Areas*binNum,numPatches);
sigma=2.4;
%computation of basic filters
GDfilters = getFilters(sigma);
%weighting functions
w0 = createGaussianFilter (sig1); 
warea1 = createGaussianFilter (sig2);
warea2 = createGaussianFilter (sig3);

[meje1,centri1]=computePoolingCenters (r1,8,65,warea1,pi/8);

for i=1:8
    dx=centri1(i,1)-round(centri1(i,1));
    dy=centri1(i,2)-round(centri1(i,2));
    w1{i} = createGaussianFilter1 (sig2,dx,dy);
end
[meje2,centri2]=computePoolingCenters (r2,8,65,warea2,0);

for i=1:8
    dx=centri1(i,1)-round(centri1(i,1));
    dy=centri1(i,2)-round(centri1(i,2));
    w2{i} = createGaussianFilter1 (sig3,dx,dy);
    norM=imfilter(ones(65,65),w2{i});
    w2{i}=w2{i}/norM(round(centri2(i,2)),round(centri2(i,1)));
end

for k=1:numPatches
    %image filtering 
    [GDmaps,Angle] = GDerivMaps(patches{k},GDfilters);
    %compute k-th descriptor
    GDdescriptor(:,k)= byParts(GDmaps,Angle,w0,w1,w2,meje1,meje2);
    
    %cut off high peaks
    for cikel=1:10      
        povprecna=mean(GDdescriptor(:,k));
        tr=2.6*povprecna;
        above=GDdescriptor(:,k)>tr;
        GDdescriptor(:,k)=above.*tr+(1-above).*GDdescriptor(:,k);
    end
    %L2 descriptor normalization
    nor=sqrt(sum(GDdescriptor(:,k).^2));
    if nor==0
        GDdescriptor(:,k)=0;
    else        
        GDdescriptor(:,k)=GDdescriptor(:,k)/nor;
    end
    %L1 descriptor normalization
    nor=sum(GDdescriptor(:,k));
    if nor==0
        GDdescriptor(:,k)=0;
    else
    GDdescriptor(:,k)=sqrt(GDdescriptor(:,k)/nor); 
    end
end
end
function [meje,centri]=computePoolingCenters (r,n,dim,filter,ang)
pCenters=zeros(n,2);
centri=zeros(n,2);
meje=zeros(n,8);
d=(size(filter,1)-1)/2;
c=dim/2+0.5;
%centers
for i=1:n
    centri(i,1) = c+(r*cos(2*pi/n*(i-1)+ang)); %x
    centri(i,2) = c-(r*sin(2*pi/n*(i-1)+ang)); %y
    pCenters(i,1)=round(centri(i,1)); %x
    pCenters(i,2)=round(centri(i,2)); %y
end

%Filter borders
for i=1:n
    ys=pCenters(i,2)-d;
    ye=pCenters(i,2)+d;
    xs=pCenters(i,1)-d;
    xe=pCenters(i,1)+d;  
    dx=1;  %filter start at
    dy=1;
    ddx=2*d+1;  %filter end at
    ddy=2*d+1;
    
    if xs<1
        dx=-xs+2; %filter
        xs=1;   %patch   
    end
    if ys<1
        dy=-ys+2; %filter
        ys=1;     %patch
    end

    if xe>dim
        ddx=2*d+1-(xe-dim); %filter
        xe=dim;      %patch
    end
    if ye>dim
        ddy=2*d+1-(ye-dim); %filter
        ye=dim; %patch
    end

    meje(i,1)=xs;
    meje(i,2)=xe;
    meje(i,3)=ys;
    meje(i,4)=ye;
    meje(i,5)=dx;
    meje(i,6)=ddx;
    meje(i,7)=dy;
    meje(i,8)=ddy;

end
end
    
function Descriptor = byParts (A,kot,w0,w1,w2,m1,m2)
binNum=8;
tip=zeros(17,binNum);
kot1=smer8(kot(:,:,1));
for i=1:8
    kot1(:,:,i)=kot1(:,:,i).*A(:,:,1);
end
d1=size(w0,1);
ds=(65-d1)/2+1;
de=ds+d1-1;
%first circle
for i=1:8 
    ww=w1{i};    
    w=ww(m1(i,7):m1(i,8),m1(i,5):m1(i,6));   
      %first derivative 
    for k=1:8 
        tip(i,k)=sum(sum((w.*kot1(m1(i,3):m1(i,4),m1(i,1):m1(i,2),k))));
    end
end

%second circle
for ii=9:16 
    i=ii-8;
    ww=w2{i};
    w=ww(m2(i,7):m2(i,8),m2(i,5):m2(i,6));
       %first derivative 
    for k=1:8
        tip(ii,k)=sum(sum(w.*kot1(m2(i,3):m2(i,4),m2(i,1):m2(i,2),k)));
    end
end
%center
% first derivative
for k=1:8
    tip(17,k)=sum(sum((w0.*kot1(ds:de,ds:de,k))));
end
Descriptor=reshape(tip,17*binNum,1);
end
function smer = smer8(kot)
[dim1,dim2]=size(kot);
kot=(kot+pi)/(pi/4);
predal=floor(kot)+1;

kot=(predal==9).*(kot-8)+(1-(predal==9)).*kot;
predal=(predal==9)+(1-(predal==9)).*predal;
smer=zeros(dim1,dim2,8);
for i=1:7
    smer(:,:,i)=smer(:,:,i)+(predal==i).*(i-kot);
    smer(:,:,i+1)=(predal==i).*(kot-(i-1));
end
    smer(:,:,8)=smer(:,:,8)+(predal==8).*(8-kot); 
    smer(:,:,1)=smer(:,:,1)+(predal==8).*(kot-7);    
end


function filt = createGaussianFilter (sigma)
   % Gaussian filter with specified sigma 
   polmer=round(2.7*sigma);
   [ x, y ] = meshgrid(-polmer:polmer,-polmer:polmer);
   filt = exp(-((x).^2/sigma^2 + (y).^2/sigma^2)/2);
   filt = filt / sum(sum(filt));
end
function [A,kot] = GDerivMaps (im,filtri)
im = double(im);
[dim1,dim2]=size(im);

A=zeros(dim1,dim2,4);
kot=zeros(dim1,dim2,4);

%first derivative 
G1x = imfilter(im, filtri{1}, 'replicate');
G1y = imfilter(im, filtri{2}, 'replicate');

kot(:,:,1)=atan2(G1y,G1x);
A(:,:,1)=cos(kot(:,:,1)).*G1x+sin(kot(:,:,1)).*G1y;
end
function filt = createGaussianFilter1 (sigma,dx,dy)
   % Gaussian filter with specified sigma 
   polmer=round(2.7*sigma);
   [ x, y ] = meshgrid(-polmer:polmer,-polmer:polmer);
   filt = exp(-((x-dx).^2/sigma^2 + (y-dy).^2/sigma^2)/2);
   filt = filt / sum(sum(filt));
end
function [filtri] = getFilters(sigma)
sze = round(3*sigma);
[x,y] = meshgrid(-sze:sze);
filtri=[];

%Gaussian filter
g0 = exp(-(x.^2/sigma^2 + y.^2/sigma^2)/2);
g0=g0/sum(sum(g0));

%first derivatives
g1x=-g0.*x/sigma^2;
g1y=-g0.*y/sigma^2;
%ww=sqrt(sum(sum(g1x.*g1x)));
filtri{1}=g1x;
filtri{2}=g1y;
end


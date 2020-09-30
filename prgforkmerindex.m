
%This prg ignore gaps and start cgr from midpoint fo the begining of a new
%protien . This program gives the k-mer index values only
filename='sequencecov2refseq.txt';
filename1='BatcoronavirusRATG13.txt';
fid = fopen(filename);
j=0;
while true
  X = fgetl(fid);
  if X~=-1
    for i=1:length(X)
      if ~strncmp(X,'>',1)
          j=j+1;   
      end 
    end
  end 
  if ~ischar(X); break; end 
end
fclose(fid);
fid = fopen(filename);
Z=zeros(j,1);
j=1;
while true
  X = fgetl(fid);
  if X~=-1
    for i=1:length(X)
      if ~strncmp(X,'>',1)
          Z(j,1)=X(i)-64;
          j=j+1;   
      end 
    end
  end 
  if ~ischar(X); break; end 
end
fclose(fid);

fid1 = fopen(filename1);
j=0;
while true
  X1 = fgetl(fid1);
  if X1~=-1
    for i=1:length(X1)
      if ~strncmp(X1,'>',1)
          j=j+1;   
      end 
    end
  end 
  if ~ischar(X1); break; end 
end
fclose(fid1);
fid1 = fopen(filename1);
Z1=zeros(j,1);
j=1;
while true
  X1 = fgetl(fid1);
  if X1~=-1
    for i=1:length(X1)
      if ~strncmp(X1,'>',1)
          Z1(j,1)=X1(i)-64;
          j=j+1;   
      end 
    end
  end 
  if ~ischar(X1); break; end 
end
fclose(fid1);

nvrtx = 4;
dnacodontotal=[730;713;711;717;330;313;390;373;371;377;111;117;171;177;2220;2203;2201;2207;2077;520;503;501;507;350;333;331;337;320;303;301;307;920;903;901;907;750;733;731;737;790;773;771;777;2050;2033;2031;2037;2030;2013;2090;2073;311;317;150;133;131;137;130;113;190;173;2011;2017;2071];
letters='NBAP';
coarsedate = Z;%reads data 1 corresponds to a 20 to t 7 to g and 3 to c
coarsedate1=Z1;

%%%%give the order
b03=3.0;
b1=4.0;b2=5.0;b3=6.0;b4=7.0;b5=8.0;b6=9.0;
bit03=2^b03;
bit1=2^b1;bit2=2^b2;
bit3=2^b3;bit4=2^b4;
bit5=2^b5;bit6=2^b6;

% Define the vertices
vrtxthree = [0 0;bit03 0;0 bit03;bit03 bit03];
vrtx =    [0 0;bit1 0;0 bit1;bit1 bit1];
vrtxfive = [0 0;bit2 0;0 bit2;bit2 bit2];
vrtxsix= [0 0;bit3 0;0 bit3;bit3 bit3];
vrtxsev =    [0 0;bit4 0;0 bit4;bit4 bit4];
vrtxeit = [0 0;bit5 0;0 bit5;bit5 bit5];
vrtxnine= [0 0;bit6 0;0 bit6;bit6 bit6];
sz=size(coarsedate);
niter = sz(1,1);  % Specify the number of iterations
sz1=size(coarsedate1);
niter1 = sz1(1,1); 


i=1;k=0;
while i<= niter
TF=int32(isnan(coarsedate(i)));
if TF==0
    
k=k+1;i=i+1;
else
    i=i+1;
end
end
niter=k;
coarse=zeros(niter,1);

i=1;k=0;
while i<= niter1
TF=int32(isnan(coarsedate1(i)));
if TF ==0
    
k=k+1;i=i+1;
else
    i=i+1;
end
end

niter1=k;
coarse1=zeros(niter1,1);

i=1;j=1;
while i<=sz(1,1)
    TF=int32(isnan(coarsedate(i)));
    if TF == 0
coarse(j)=coarsedate(i);
i=i+1;j=j+1;
    else
  i=i+1;      
    end
end

i=1;j=1;
while i<=sz1(1,1)
   TF=int32(isnan(coarsedate1(i)));
    if TF == 0
coarse1(j)=coarsedate1(i);
i=i+1;j=j+1;
    else
  i=i+1;      
    end
end


A=int32(niter);
B=int32(3);
codonniter=idivide(A,B);
codoncoarsevalue=zeros(codonniter,1);

t=1;u=2;v=3;i=1;
while v <= niter
codoncoarsevalue(i)=(coarse(t)*100) +(coarse(u)*10)+ coarse(v);
 i=i+1;t=t+3;u=u+3;v=v+3;
end



codonsize=size(codoncoarsevalue);
codonsz=codonsize(1,1);
codonabnpcount=0;
for i=1:codonsz
    for j=1:64
        if codoncoarsevalue(i)==dnacodontotal(j)
            codonabnpcount=codonabnpcount+1;
        end
    end
end
    


codonniter=codonabnpcount;
codoncoarse=zeros(codonsz,1);nostopcodon=0;
codonbin=[1;2;3;4;5];
for i=1:codonsz
    for j= 1:4
        if codoncoarsevalue(i)==dnacodontotal(j)
            codoncoarse(i)=codonbin(3);
        end
    end
    for j=5:14
         if codoncoarsevalue(i)==dnacodontotal(j)
            codoncoarse(i)=codonbin(2);
         end
    end
  for j=15:43
         if codoncoarsevalue(i)==dnacodontotal(j)
            codoncoarse(i)=codonbin(1);
         end
  end
  for j=44:61
         if codoncoarsevalue(i)==dnacodontotal(j)
            codoncoarse(i)=codonbin(4);
         end
         
  end
  for j=62:64
         if codoncoarsevalue(i)==dnacodontotal(j)
            codoncoarse(i)=codonbin(5);
            nostopcodon=nostopcodon+1;
         end
  end
  
         
end   

i=1;j=1;
codoncoarsenew=zeros(codonniter,1);
while i<=codonsz
    if  codoncoarse(i)~=0
        codoncoarsenew(j)=codoncoarse(i);
        i=i+1;
        j=j+1;
    else
        i=i+1;
    end
end
coarse=codoncoarsenew;
    niter=codonniter;
    

no1=0;no2=0;no3=0;no4=0;no5=0;
for i=1:niter
    if coarse(i)==1
        no1=no1+1;
    elseif coarse(i)==2
        no2=no2+1;
    elseif coarse(i)==3
        no3=no3+1;
    elseif coarse(i)==4
        no4=no4+1;
    else
        no5=no5+1;
    end
end
percentno1=no1*100/(niter-no5);
percentno2=no2*100/(niter-no5);
percentno3=no3*100/(niter-no5);
percentno4=no4*100/(niter-no5);


%three bit

pointthree = zeros(niter,2) ;
pointthree(1,1)=bit03/2;pointthree(1,2)=bit03/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    if vIdx~=5
    pointthree(i,:) = vrtxthree(vIdx,:) - (vrtxthree(vIdx,:) - pointthree(i-1,:))/2;
    else
     pointthree(i,1)=bit03/2;pointthree(i,2)=bit03/2;   
    end
end


threebit= zeros(bit03,bit03);
for i=1:niter
    for k=1:bit03
        if (pointthree(i,1)<=k && pointthree(i,1)>(k-1))
            for j=1:bit03
                if( pointthree(i,2)<=j && pointthree(i,2)>(j-1))
                    threebit(k,j)=threebit(k,j)+1;
                end
            end
            
            
        end
    end
end
threebit(bit03/2,bit03/2)=threebit(bit03/2,bit03/2)-nostopcodon;
threebitp=threebit;% will be filled or notfilled ie 0 or 1.
threebitp(logical(threebitp))=1;


% four bit cgr

point = zeros(niter,2);
point(1,1)=bit1/2;point(1,2)=bit1/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    if vIdx~=5
    point(i,:) = vrtx(vIdx,:) - (vrtx(vIdx,:) - point(i-1,:))/2;
    else
     point(i,1)=bit1/2;point(i,2)=bit1/2;
    end
    
end

frbit= zeros(bit1,bit1);
for i=1:niter
    for k=1:bit1
        if (point(i,1)<=k && point(i,1)>(k-1))
            for j=1:bit1
                if( point(i,2)<=j && point(i,2)>(j-1))
                    frbit(k,j)=frbit(k,j)+1;
                end
            end
            
            
        end
    end
end
frbit(bit1/2,bit1/2)=frbit(bit1/2,bit1/2)-nostopcodon;
frbitp=frbit;%frbitp will be filled or notfilled ie 0 or 1.
frbitp(logical(frbitp))=1;

%five bit

pointfive = zeros(niter,2) ;
pointfive(1,1)=bit2/2;pointfive(1,2)=bit2/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    if vIdx~=5
    pointfive(i,:) = vrtxfive(vIdx,:) - (vrtxfive(vIdx,:) - pointfive(i-1,:))/2;
    else
      pointfive(i,1)=bit2/2;pointfive(i,2)=bit2/2;
    end
end

fivebit= zeros(bit2,bit2);
for i=1:niter
    for k=1:bit2
        if (pointfive(i,1)<=k && pointfive(i,1)>(k-1))
            for j=1:bit2
                if( pointfive(i,2)<=j && pointfive(i,2)>(j-1))
                    fivebit(k,j)=fivebit(k,j)+1;
                end
            end
            
            
        end
    end
end

fivebit(bit2/2,bit2/2)=fivebit(bit2/2,bit2/2)-nostopcodon;
fivebitp=fivebit;% will be filled or notfilled ie 0 or 1.
fivebitp(logical(fivebitp))=1;

%sixbit
pointsix = zeros(niter,2) ;
pointsix(1,1)=bit3/2;pointsix(1,2)=bit3/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
    if vIdx~=5
    pointsix(i,:) = vrtxsix(vIdx,:) - (vrtxsix(vIdx,:) - pointsix(i-1,:))/2;
    else
      pointsix(i,1)=bit3/2;pointsix(i,2)=bit3/2;
    end  
end

sixbit= zeros(bit3,bit3);
for i=1:niter
    for k=1:bit3
        if (pointsix(i,1)<=k && pointsix(i,1)>(k-1))
            for j=1:bit3
                if( pointsix(i,2)<=j && pointsix(i,2)>(j-1))
                    sixbit(k,j)=sixbit(k,j)+1;
                end
            end
            
            
        end
    end
end

sixbit(bit3/2,bit3/2)=sixbit(bit3/2,bit3/2)-nostopcodon;
sixbitp=sixbit;% will be filled or notfilled ie 0 or 1.
sixbitp(logical(sixbitp))=1;

%sevenbit
pointseven = zeros(niter,2) ;
pointseven(1,1)=bit4/2;pointseven(1,2)=bit4/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
     if vIdx~=5
    pointseven(i,:) = vrtxsev(vIdx,:) - (vrtxsev(vIdx,:) - pointseven(i-1,:))/2;
     else
     pointseven(i,1)=bit4/2;pointseven(i,2)=bit4/2;
     end
end

sevenbit= zeros(bit4,bit4);
for i=1:niter
    for k=1:bit4
        if (pointseven(i,1)<=k && pointseven(i,1)>(k-1))
            for j=1:bit4
                if( pointseven(i,2)<=j && pointseven(i,2)>(j-1))
                    sevenbit(k,j)=sevenbit(k,j)+1;
                end
            end
            
            
        end
    end
end
sevenbit(bit4/2,bit4/2)=sevenbit(bit4/2,bit4/2)-nostopcodon;
sevenbitp=sevenbit;% will be filled or notfilled ie 0 or 1.
sevenbitp(logical(sevenbitp))=1;

%eitbit
pointeit = zeros(niter,2) ;
pointeit(1,1)=bit5/2;pointeit(1,2)=bit5/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
     if vIdx~=5
    pointeit(i,:) = vrtxeit(vIdx,:) - (vrtxeit(vIdx,:) - pointeit(i-1,:))/2;
     else
      pointeit(i,1)=bit5/2;pointeit(i,2)=bit5/2;
     end
end
eitbit= zeros(bit5,bit5);
% k=1;
% j=1;

for i=1:niter
    for k=1:bit5
        if (pointeit(i,1)<=k && pointeit(i,1)>(k-1))
            for j=1:bit5
                if( pointeit(i,2)<=j && pointeit(i,2)>(j-1))
                    eitbit(k,j)=eitbit(k,j)+1;
                end
            end
            
            
        end
    end
end
eitbit(bit5/2,bit5/2)=eitbit(bit5/2,bit5/2)-nostopcodon;
eitbitp=eitbit; %same as frbitp
eitbitp(logical(eitbit))=1;

% nine bit
pointnine = zeros(niter,2) ;
pointnine(1,1)=bit6/2;pointnine(1,2)=bit6/2;
for i = 2:niter                                              % Generate the points
    vIdx = coarse(i-1,1);
     if vIdx~=5
    pointnine(i,:) = vrtxnine(vIdx,:) - (vrtxnine(vIdx,:) - pointnine(i-1,:))/2;
     else
     pointnine(i,1)=bit6/2;pointnine(i,2)=bit6/2;
     end
end

ninebit= zeros(bit6,bit6);
for i=1:niter
    for k=1:bit6
        if (pointnine(i,1)<=k && pointnine(i,1)>(k-1))
            for j=1:bit6
                if( pointnine(i,2)<=j && pointnine(i,2)>(j-1))
                    ninebit(k,j)=ninebit(k,j)+1;
                end
            end
            
            
        end
    end
end
ninebit(bit6/2,bit6/2)=ninebit(bit6/2,bit6/2)-nostopcodon;
ninebitp=ninebit;% will be filled or notfilled ie 0 or 1.
ninebitp(logical(ninebitp))=1;


A1=int32(niter1);
B1=int32(3);
codonniter1=idivide(A1,B1);
codoncoarsevalue1=zeros(codonniter1,1);

t=1;u=2;v=3;i=1;
while v <= niter1
codoncoarsevalue1(i)=(coarse1(t)*100) +(coarse1(u)*10)+ coarse1(v);
 i=i+1;t=t+3;u=u+3;v=v+3;
end


codonsize1=size(codoncoarsevalue1);
codonsz1=codonsize1(1,1);nostopcodon1=0;
codonabnpcount1=0;
for i=1:codonsz1
    for j=1:64
        if codoncoarsevalue1(i)==dnacodontotal(j)
            codonabnpcount1=codonabnpcount1+1;
        end
    end
end
    
codonniter1=codonabnpcount1;
codoncoarse1=zeros(codonsz1,1);

for i=1:codonsz1
    for j= 1:4
        if codoncoarsevalue1(i)==dnacodontotal(j)
            codoncoarse1(i)=codonbin(3);
        end
    end
    for j=5:14
         if codoncoarsevalue1(i)==dnacodontotal(j)
            codoncoarse1(i)=codonbin(2);
         end
    end
  for j=15:43
         if codoncoarsevalue1(i)==dnacodontotal(j)
            codoncoarse1(i)=codonbin(1);
         end
  end
  for j=44:61
         if codoncoarsevalue1(i)==dnacodontotal(j)
            codoncoarse1(i)=codonbin(4);
         end
         
  end
  
  for j=62:64
         if codoncoarsevalue1(i)==dnacodontotal(j)
            codoncoarse1(i)=codonbin(5);
            nostopcodon1=nostopcodon1+1;
         end
         
  end     
end


i=1;j=1;
codoncoarsenew1=zeros(codonniter1,1);
while i<=codonsz1
    if  codoncoarse1(i)~=0
        codoncoarsenew1(j)=codoncoarse1(i);
        i=i+1;
        j=j+1;
    else
        i=i+1;
    end
end
coarse1=codoncoarsenew1;
    niter1=codonniter1;
    
    %n01-no4 specify the number of CATand G respectivley
noo1=0;noo2=0;noo3=0;noo4=0;noo5=0;
for i=1:niter1
    if coarse1(i)==1
        noo1=noo1+1;
    elseif coarse1(i)==2
        noo2=noo2+1;
    elseif coarse1(i)==3
        noo3=noo3+1;
    elseif coarse1(i)==4
        noo4=noo4+1;
    else
        noo5=noo5+1;
    end
end

percent1noo1=noo1*100/(niter1-noo5);
percent1noo2=noo2*100/(niter1-noo5);
percent1noo3=noo3*100/(niter1-noo5);
percent1noo4=noo4*100/(niter1-noo5);




pointthree1 = zeros(niter1,2) ;
pointthree1(1,1)=bit03/2;pointthree1(1,2)=bit03/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    if vIdx1~=5
    pointthree1(i,:) = vrtxthree(vIdx1,:) - (vrtxthree(vIdx1,:) - pointthree1(i-1,:))/2;
    else
     pointthree1(i,1)=bit03/2;pointthree1(i,2)=bit03/2; 
    end
end


threebit1= zeros(bit03,bit03);
for i=1:niter1
    for k=1.00:bit03
        if (pointthree1(i,1)<=k && pointthree1(i,1)>(k-1))
            for j=1.00:bit03
                if( pointthree1(i,2)<=j && pointthree1(i,2)>(j-1))
                    threebit1(k,j)=threebit1(k,j)+1;
                end
            end
            
            
        end
    end
end
threebit1(bit03/2,bit03/2)=threebit1(bit03/2,bit03/2)-nostopcodon1;
threebitp1=threebit1;% will be filled or notfilled ie 0 or 1.
threebitp1(logical(threebitp1))=1;

%four bit

point1 = zeros(niter1,2);
point1(1,1)=bit1/2;point1(1,2)=bit1/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    if vIdx1~=5
    point1(i,:) = vrtx(vIdx1,:) - (vrtx(vIdx1,:) - point1(i-1,:))/2;
    else
    point1(i,1)=bit1/2;point1(i,2)=bit1/2;
    end
end

frbit1= zeros(bit1,bit1);
for i=1:niter1
    for k=1:bit1
        if (point1(i,1)<=k && point1(i,1)>(k-1))
            for j=1:bit1
                if( point1(i,2)<=j && point1(i,2)>(j-1))
                    frbit1(k,j)=frbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
frbit1(bit1/2,bit1/2)=frbit1(bit1/2,bit1/2)-nostopcodon1;
frbitp1=frbit1;%frbitp will be filled or notfilled ie 0 or 1.
frbitp1(logical(frbitp1))=1;

%five bit

pointfive1 = zeros(niter1,2) ;
pointfive1(1,1)=bit2/2;pointfive1(1,2)=bit2/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    if vIdx1~=5
    pointfive1(i,:) = vrtxfive(vIdx1,:) - (vrtxfive(vIdx1,:) - pointfive1(i-1,:))/2;
    else
     pointfive1(i,1)=bit2/2;pointfive1(i,2)=bit2/2;
    end
end

fivebit1= zeros(bit2,bit2);
for i=1:niter1
    for k=1:bit2
        if (pointfive1(i,1)<=k && pointfive1(i,1)>(k-1))
            for j=1:bit2
                if( pointfive1(i,2)<=j && pointfive1(i,2)>(j-1))
                    fivebit1(k,j)=fivebit1(k,j)+1;
                end
            end
            
            
        end
    end
end
fivebit1(bit2/2,bit2/2)=fivebit1(bit2/2,bit2/2)-nostopcodon1;
fivebitp1=fivebit1;% will be filled or notfilled ie 0 or 1.
fivebitp1(logical(fivebitp1))=1;

%sixbit
pointsix1 = zeros(niter1,2) ;
pointsix1(1,1)=bit3/2;pointsix1(1,2)=bit3/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    if vIdx1~=5
    pointsix1(i,:) = vrtxsix(vIdx1,:) - (vrtxsix(vIdx1,:) - pointsix1(i-1,:))/2;
    else
    pointsix1(i,1)=bit3/2;pointsix1(i,2)=bit3/2;
    end
end

sixbit1= zeros(bit3,bit3);
for i=1:niter1
    for k=1:bit3
        if (pointsix1(i,1)<=k && pointsix1(i,1)>(k-1))
            for j=1:bit3
                if( pointsix1(i,2)<=j && pointsix1(i,2)>(j-1))
                    sixbit1(k,j)=sixbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
sixbit1(bit3/2,bit3/2)=sixbit1(bit3/2,bit3/2)-nostopcodon1;
sixbitp1=sixbit1;% will be filled or notfilled ie 0 or 1.
sixbitp1(logical(sixbitp1))=1;

%sevenbit
pointseven1 = zeros(niter1,2) ;
pointseven1(1,1)=bit4/2;pointseven1(1,2)=bit4/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
     if vIdx1~=5
    pointseven1(i,:) = vrtxsev(vIdx1,:) - (vrtxsev(vIdx1,:) - pointseven1(i-1,:))/2;
     else
     pointseven1(i,1)=bit4/2;pointseven1(i,2)=bit4/2;
     end
end

sevenbit1= zeros(bit4,bit4);
for i=1:niter1
    for k=1:bit4
        if (pointseven1(i,1)<=k && pointseven1(i,1)>(k-1))
            for j=1:bit4
                if( pointseven1(i,2)<=j && pointseven1(i,2)>(j-1))
                    sevenbit1(k,j)=sevenbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
sevenbit1(bit4/2,bit4/2)=sevenbit1(bit4/2,bit4/2)-nostopcodon1;
sevenbitp1=sevenbit1;% will be filled or notfilled ie 0 or 1.
sevenbitp1(logical(sevenbitp1))=1;

%eitbit
pointeit1 = zeros(niter1,2) ;
pointeit1(1,1)=bit5/2;pointeit1(1,2)=bit5/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
    if vIdx1~=5
    pointeit1(i,:) = vrtxeit(vIdx1,:) - (vrtxeit(vIdx1,:) - pointeit1(i-1,:))/2;
    else
    pointeit1(i,1)=bit5/2;pointeit1(i,2)=bit5/2;
    end
end

eitbit1= zeros(bit5,bit5);
% k=1;
% j=1;

for i=1:niter1
    for k=1:bit5
        if (pointeit1(i,1)<=k && pointeit1(i,1)>(k-1))
            for j=1:bit5
                if( pointeit1(i,2)<=j && pointeit1(i,2)>(j-1))
                    eitbit1(k,j)=eitbit1(k,j)+1;
                end
            end
            
            
        end
    end
end
eitbit1(bit5/2,bit5/2)=eitbit1(bit5/2,bit5/2)-nostopcodon1;
eitbitp1=eitbit1; %same as frbitp
eitbitp1(logical(eitbit1))=1;


% nine bit
pointnine1 = zeros(niter1,2) ;
pointnine1(1,1)=bit6/2;pointnine1(1,2)=bit6/2;
for i = 2:niter1                                              % Generate the points
    vIdx1 = coarse1(i-1,1);
     if vIdx1~=5
    pointnine1(i,:) = vrtxnine(vIdx1,:) - (vrtxnine(vIdx1,:) - pointnine1(i-1,:))/2;
     else
      pointnine1(i,1)=bit6/2;pointnine1(i,2)=bit6/2; 
     end
end

ninebit1= zeros(bit6,bit6);
for i=1:niter1
    for k=1:bit6
        if (pointnine1(i,1)<=k && pointnine1(i,1)>(k-1))
            for j=1:bit6
                if( pointnine1(i,2)<=j && pointnine1(i,2)>(j-1))
                    ninebit1(k,j)=ninebit1(k,j)+1;
                end
            end
            
            
        end
    end
end
ninebit1(bit6/2,bit6/2)=ninebit1(bit6/2,bit6/2)-nostopcodon1;
ninebitp1=ninebit1;% will be filled or notfilled ie 0 or 1.
ninebitp1(logical(ninebitp1))=1;

threebitpercent=threebit*100/(niter-nostopcodon);threebitpercent1=threebit1*100/(niter1-nostopcodon1);
frbitpercent=frbit*100/(niter-nostopcodon); frbitpercent1=frbit1*100/(niter1-nostopcodon1);
fivebitpercent=fivebit*100/(niter-nostopcodon); fivebitpercent1=fivebit1*100/(niter1-nostopcodon1);
sixbitpercent=sixbit*100/(niter-nostopcodon); sixbitpercent1=sixbit1*100/(niter1-nostopcodon1);
sevenbitpercent=sevenbit*100/(niter-nostopcodon); sevenbitpercent1=sevenbit1*100/(niter1-nostopcodon1);
eitbitpercent=eitbit*100/(niter-nostopcodon); eitbitpercent1=eitbit1*100/(niter1-nostopcodon1);
ninebitpercent=ninebit*100/(niter-nostopcodon); ninebitpercent1=ninebit1*100/(niter1-nostopcodon1);

sthreebit=threebitpercent-threebitpercent1;
sfrbit=frbitpercent-frbitpercent1;
sfivebit=fivebitpercent-fivebitpercent1;
ssixbit=sixbitpercent-sixbitpercent1;
ssevenbit=sevenbitpercent-sevenbitpercent1;
seitbit=eitbitpercent-eitbitpercent1;
sninebit=ninebitpercent-ninebitpercent1;



posindexthreebit=0;
negindexthreebit=0;
for i= 1:bit03
    for j= 1:bit03
 if (sthreebit(i,j)<=0)
     negindexthreebit=negindexthreebit+sthreebit(i,j);
 else
posindexthreebit=posindexthreebit+sthreebit(i,j);
 end
    end
end

posindexfrbit=0;
negindexfrbit=0;


for i= 1:bit1
    for j= 1:bit1
 if (sfrbit(i,j)<=0)
     negindexfrbit=negindexfrbit+sfrbit(i,j);
 else
posindexfrbit=posindexfrbit+sfrbit(i,j);
 end
    end
end

posindexfivebit=0;
negindexfivebit=0;


for i= 1:bit2
    for j= 1:bit2
 if (sfivebit(i,j)<=0)
     negindexfivebit=negindexfivebit+sfivebit(i,j);
 else
posindexfivebit=posindexfivebit+sfivebit(i,j);
 end
    end
end

posindexsixbit=0;
negindexsixbit=0;


for i= 1:bit3
    for j= 1:bit3
 if (ssixbit(i,j)<=0)
     negindexsixbit=negindexsixbit+ssixbit(i,j);
 else
posindexsixbit=posindexsixbit+ssixbit(i,j);
 end
    end
end

posindexsevenbit=0;
negindexsevenbit=0;


for i= 1:bit4
    for j= 1:bit4
 if (ssevenbit(i,j)<=0)
     negindexsevenbit=negindexsevenbit+ssevenbit(i,j);
 else
posindexsevenbit=posindexsevenbit+ssevenbit(i,j);
 end
    end
end

posindexeitbit=0;
negindexeitbit=0;


for i= 1:bit5
    for j= 1:bit5
 if (seitbit(i,j)<=0)
     negindexeitbit=negindexeitbit+seitbit(i,j);
 else
posindexeitbit=posindexeitbit+seitbit(i,j);
 end
    end
end

posindexninebit=0;
negindexninebit=0;


for i= 1:bit6
    for j= 1:bit6
 if (sninebit(i,j)<=0)
     negindexninebit=negindexninebit+sninebit(i,j);
 else
posindexninebit=posindexninebit+sninebit(i,j);
 end
    end
end
%%%% matrx with similarity index %%%%%%%
similarityindexmatr(1,1)=posindexthreebit;similarityindexmatr(2,1)=negindexthreebit;
similarityindexmatr(1,2)=posindexfrbit;similarityindexmatr(2,2)=negindexfrbit;
similarityindexmatr(1,3)=posindexfivebit;similarityindexmatr(2,3)=negindexfivebit;
similarityindexmatr(1,4)=posindexsixbit;similarityindexmatr(2,4)=negindexsixbit;
similarityindexmatr(1,5)=posindexsevenbit;similarityindexmatr(2,5)=negindexsevenbit;
similarityindexmatr(1,6)=posindexeitbit;similarityindexmatr(2,6)=negindexeitbit;
similarityindexmatr(1,7)=posindexninebit;similarityindexmatr(2,7)=negindexninebit;

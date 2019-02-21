clc;
close all;
clear all;
syms 'x';syms 'y';
%----------------------------------------------------------------------------
%   plane stress analysis of a solid using Q4            
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   disp = system nodal displacement vector
%   eldisp = element nodal displacement vector
%   stress = matrix containing stresses
%   strain = matrix containing strains
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

nel=39;                  % number of elements 12
nnel=4;                  % number of nodes per element
ndof=2;                  % number of dofs per node
nnode=56;                % total number of nodes in system 21
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=2e11;            % elastic modulus
poisson=0.3;             % Poisson's ratio
Length=5;
D=1;
P=10000;
MoI=(D^3)/12;
q=-(P/(2*MoI))*(((D^2)/4)-(y^2));
%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------
gcoord=xlsread('Data.xlsx',1,'B2:C57');
nodes=xlsread('Data.xlsx',2,'B1:E39');

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------
%-------------------------------------
%  input data for boundary conditions
%-------------------------------------
syms 'y1';
U(y1)=(-2.300092004*10^-7)*((y1^3)-(y1/4));
V(y1)=(-4.50180072*10^-7)*(y1^2);
bcdof=[35 36 61 62 63 64 1 2];         % first two dofs are constrained
ydiscoor=xlsread('Data.xlsx',3,'B1:B4');           % whose described values are 0 
bcval=zeros(1,length(bcdof));
for t=1:length(ydiscoor)
    if(t==1)
    bcval(t)=subs(U(y1),y1,gcoord(ydiscoor(t),2));
    bcval(t+1)=subs(V(y1),y1,gcoord(ydiscoor(t),2));    
    else
bcval(t+1)=subs(U(y1),y1,gcoord(ydiscoor(t),2));
bcval(t+2)=subs(V(y1),y1,gcoord(ydiscoor(t),2));
    end
end
%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
displ=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nel,4);    % matrix containing stress components
strain=zeros(nel,4);    % matrix containing strain components
index=zeros(edof,1);    % index vector
Bmtx2=zeros(3,edof);    % kinematic matrix
matmtx=zeros(3,3);      % constitutive matrix
NodeStress=zeros(3,nnode);
NodeStrain=zeros(3,nnode);
NumEle=3;
% a=D/NumEle*2;
%----------------------------
%  force vector
%----------------------------

ForceNodes=xlsread('Data.xlsx',3,'A1:A4');
Force1=(1/18)*[2 1;1 2]*[subs(q,y,gcoord(ForceNodes(1),2));subs(q,y,gcoord(ForceNodes(2),2))];
Force2=(1/18)*[2 1;1 2]*[subs(q,y,gcoord(ForceNodes(2),2));subs(q,y,gcoord(ForceNodes(3),2))];
Force3=(1/18)*[2 1;1 2]*[subs(q,y,gcoord(ForceNodes(3),2));subs(q,y,gcoord(ForceNodes(4),2))];
ForceEqu=[Force1(1),Force1(2)+Force2(1),Force2(2)+Force3(1),Force3(2)];
ff(ForceNodes(1:length(ForceNodes))*2)=ForceEqu(1:4);
disp(double(ForceEqu));                       
% %-----------------------------------------------------------------
% %  computation of element matrices and vectors and their assembly
% %-----------------------------------------------------------------
% 
matmtx=fematiso(1,emodule,poisson);        % compute constitutive matrix
ShapFun(1)=(1/4)*(1-x)*(1-y);
ShapFun(2)=(1/4)*(1+x)*(1-y);
ShapFun(3)=(1/4)*(1+x)*(1+y);
ShapFun(4)=(1/4)*(1-x)*(1+y);
DiffMatr=[diff(ShapFun(1),x),0,diff(ShapFun(2),x),0,diff(ShapFun(3),x),0,diff(ShapFun(4),x),0;
          diff(ShapFun(1),y),0,diff(ShapFun(2),y),0,diff(ShapFun(3),y),0,diff(ShapFun(4),y),0;
          0,diff(ShapFun(1),x),0,diff(ShapFun(2),x),0,diff(ShapFun(3),x),0,diff(ShapFun(4),x);
          0,diff(ShapFun(1),y),0,diff(ShapFun(2),y),0,diff(ShapFun(3),y),0,diff(ShapFun(4),y)];
prompt='Enter the desired Gaussian Quadrature points 1 or 2 or 3 = ';
point=input(prompt);
for iel=1:nel % loop for the total number of elements
nd=nodes(iel,1:4);
for j=1:nnel
Xcoord(j)=gcoord(nd(j),1);
Ycoord(j)=gcoord(nd(j),2);% coord values of 1st node
end
index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
%-------------------------------------------------------
%  find the derivatives of shape functions
%-------------------------------------------------------
coeffMat=[1 0 0 0 ;0 0 0 1;0 1 1 0];
[Jacobian]=fekine2d(Xcoord,Ycoord,ShapFun);
TransfMatrx=[inv(Jacobian) zeros(2,2);zeros(2,2) inv(Jacobian)];
Bmtx2=coeffMat*TransfMatrx*DiffMatr;       % compute kinematic matrix
Bmtx3=Bmtx2;
[k]=EleStiffCalc(Bmtx2,matmtx,Jacobian,point);% element stiffnes matrix
kk=feasmbl1(kk,k,index);                   % assemble element matrices 
end

%-----------------------------
%   apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);
%----------------------------
%  solve the matrix equation
%----------------------------

displ=kk\ff;   
DisplX=displ(1:2:nnode*2-1);
DisplY=displ(2:2:nnode*2);
%---------------------------------------
%  element stress computation
%---------------------------------------

for ielp=1:nel                                       % loop for the total number of elements
for j=1:nnel
nd(j)=nodes(ielp,j);                                 % 1st connected node for (iel)-th element
Xcoord(j)=gcoord(nd(j),1); Ycoord(j)=gcoord(nd(j),2);% coord values of 1st node
end
index=feeldof(nd,nnel,ndof);                         % extract system dofs associated with element
%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------
t=sqrt(3);
MultMatr=[-1 -1;1 -1;1 1;-1 1];
for i=1:edof
eldisp(i)=displ(index(i));
end
estr=Bmtx2*eldisp;
% disp('The element '+ielp);
for j=1:4
   NodeStrain(1:3,nd(j))=subs(estr,[x,y],[MultMatr(j,1),MultMatr(j,2)]);
   NodeStress(1:3,nd(j))=matmtx*NodeStrain(1:3,nd(j));
end
t=sqrt(3)/2;
s=1/2;
% conver=[1+t -s 1-t -s;-s 1+t -s 1-t;1-t -s 1+t -s;-s 1-t -s 1+t];
estrain=eleStrain(MultMatr,DiffMatr,coeffMat,eldisp,Xcoord,Ycoord,ShapFun);
ElemeStrains(1:3,ielp)=estrain;
ElemeStresses(1:3,ielp)=matmtx*estrain;
end
disp('The Elemental Strains are');
disp(double(ElemeStrains))
disp('The Elemental Stresses are');
disp(double(ElemeStresses));
disp('The Nodal Strains are');
disp(double(NodeStrain));
disp('The Nodal Stresses are');
disp(double(NodeStress));
disp('The Nodal Displacements are');
disp('Ux     Uy');
disp([DisplX,DisplY]);
%print stresses
%---------------------------------------------------------------





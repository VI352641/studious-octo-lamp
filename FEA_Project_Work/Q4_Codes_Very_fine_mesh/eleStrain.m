function estrain=eleStrain(MultMatr,DiffMatr,coeffMat,eldisp,Xcoord,Ycoord,ShapFun)
syms 'x';syms 'y';
DiffMatr1=DiffMatr;
o=[1,3,5,7];e=[2,4,6,8];
for m=1:4
    if(m==1||m==2)
    for n=1:4
DiffMatr1(m,o(n))=subs(DiffMatr1(m,o(n)),[x,y],[MultMatr(n,1),MultMatr(n,2)]);
    end
    else
    for n=1:4
DiffMatr1(m,e(n))=subs(DiffMatr1(m,e(n)),[x,y],[MultMatr(n,1),MultMatr(n,2)]);
    end    
    end
end
CoordMatr=zeros(4,2);
for i=1:4
CoordMatr(i,1)=Xcoord(i);
CoordMatr(i,2)=Ycoord(i);
end
DiffShap=[diff(ShapFun(1),x),diff(ShapFun(2),x),diff(ShapFun(3),x),diff(ShapFun(4),x);
          diff(ShapFun(1),y),diff(ShapFun(2),y),diff(ShapFun(3),y),diff(ShapFun(4),y)];
      for j=1:2
      for i=1:4
          mum=DiffShap(j,i);
          Diffshap(j,i)=subs(mum,x,MultMatr(i,1));
          Diffshap(j,i)=subs(mum,y,MultMatr(i,2));
      end
      end
[Jacobian1]=DiffShap*CoordMatr;
TransfMatrx1=[inv(Jacobian1) zeros(2,2);zeros(2,2) inv(Jacobian1)];
[estrain]=coeffMat*TransfMatrx1*DiffMatr1*eldisp;





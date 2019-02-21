function [k]=EleStiffCalc(Bmtx2,matmtx,Jacobian,point)
syms 'x';syms 'y';
[k]=Bmtx2'*matmtx*Bmtx2*det(Jacobian);
[X,W]=lgwt(point,-1,1);
switch point
    case 1
    for i=1:length(k)
    for j=1:length(k)
        k(i,j)=W(1)*subs(W(1)*subs(k(i,j),x,X(1)),y,X(1));
    end
    end
    case 2   
for i=1:length(k)
for j=1:length(k)
  k(i,j)= W(1)*subs(W(1)*subs(k(i,j),x,X(1))+W(2)*subs(k(i,j),x,X(2)),y,X(1))+W(2)*subs(W(1)*subs(k(i,j),x,X(1))+W(2)*subs(k(i,j),x,X(2)),y,X(2));
end 
end
    case 3
    for i=1:length(k)
    for j=1:length(k)
     Val=W(1)*subs(k(i,j),x,X(1))+W(2)*subs(k(i,j),x,X(2))+W(3)*subs(k(i,j),x,X(3));
     k(i,j)=W(1)*subs(Val,y,X(1))+W(2)*subs(Val,y,X(2))+W(3)*subs(Val,y,X(3));
    end
    end
end

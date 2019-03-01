AA=load('ResultsLog.txt');
%%
Al=AA(2:end,1); C_Val=AA(2:end,2); Stab=AA(2:end,3);
i=1;
Al2=Al(1);
%%
while Al2<=max(Al)
Stab2=Stab(Al==Al2);
C2=C_Val(Al==Al2);
[m,idx]=min(Stab2);
[n,idx2]=min(C2);
K(i,1:3)=[Al2,Stab(idx2),n];
Z(i,1:3)=[Al2,m,C2(idx)];
i=i+1;
C_Val=C_Val(Al>Al2);
Stab=Stab(Al>Al2);
Al=Al(Al>Al2);
Al2=Al(1);
end
%%
subplot(2,1,1)
plot(K(:,1),K(:,3))
title('Min C')
subplot(2,1,2)
loglog(Z(:,1),Z(:,3))
title('Min Stab')

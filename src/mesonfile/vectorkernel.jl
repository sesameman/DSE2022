timetest=time()
# 计算
# i & s 代表外动量
# j & m 代表内动量

Threads.@threads for i=1:length(meshk)
    k2=meshk[i]
    for s=1:length(meshz)
        zk=meshz[s]
        for j=1:length(meshk)
            q2=meshk[j]
            for m=1:length(meshz)
                zq=meshz[m]
                weightzq=weightz[m]
                weightq2=weightk[j]
                A1in=A1[j,m]
                A2in=A2[j,m]
                B1in=B1[j,m]
                B2in=B2[j,m]
                allweight=-weightzq*weightq2*q2/(16*pi^3)/((P2/4+q2+sqrt(P2*q2)*zq)*A1in^2+B1in^2)/((P2/4+q2-sqrt(P2*q2)*zq)*A2in^2+B2in^2)*z2^2

                kernel11(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*A1in*A2in*P2*q2^2 + k2*P2*(8*B1in*B2in - 2*A1in*A2in*(P2 - 4*q2))*zk^2 - k2*(-5*A1in*A2in*P2^2 + 4*P2*(5*B1in*B2in + 3*A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2) + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(-4*B1in*B2in + A1in*A2in*(P2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 2*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(20*B1in*B2in + A1in*A2in*(-5*P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 2*P2*q2*zq^2*(4*B1in*B2in + A1in*A2in*(-P2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) - q2*(16*A1in*A2in*P2*q2*zq^2 + P2*(20*B1in*B2in + A1in*A2in*(-5*P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel12(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(4*sqrt(k2*P2)*(4*B1in*B2in - A1in*A2in*P2)*sqrt(P2*q2)*zk*zq*(-(P2*q2) + P2*q2*zq^2) - k2*P2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*zk^2*(P2*q2 + 2*P2*q2*zq^2) + k2*P2*q2*(3*A1in*A2in*P2^2 + 4*P2*(-3*B1in*B2in + A1in*A2in*q2) - 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-2*(-(P2*q2) + P2*q2*zq^2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in - 2*A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2) + sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(-3*A1in*A2in*P2^2 + 4*P2*(3*B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)) + (-(P2*q2) + P2*q2*zq^2)*(2*(4*B1in*B2in - A1in*A2in*(P2 + 8*q2))*(P2*q2 - P2*q2*zq^2) + k2*(A1in*A2in*P2^2 - 4*P2*(B1in*B2in + 3*A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2))))/(9. *P2^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel13(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*sqrt(P2*q2)*zq*(k2*P2*sqrt(P2*q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*zk^2*zq + sqrt(k2*P2)*zk*(P2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) - 2*P2*(4*B1in*B2in + A1in*A2in*P2)*q2*zq^2 - sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(P2*q2)*zq*((4*B1in*B2in + A1in*A2in*(-4*k2 + P2 - 8*q2))*(-(P2*q2) + P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in - 5*A1in*A2in*q2) + 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 8*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel14(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*A1in*A2in*(k2*P2*q2*zk^2 + (2*k2 + q2)*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(q2)*(-(P2*q2) + sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(9. *(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel15(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*((A2in*B1in - A1in*B2in)*k2*P2*sqrt(P2*q2)*zk^2*zq - ((A2in*B1in + A1in*B2in)*k2 + 3*(A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(4*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + P2*(-6*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 2*(A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + sqrt(k2*P2)*zk*(2*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in - A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 - 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel16(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-16*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(k2*P2*q2*zk^2 + (2*k2 + q2)*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(q2)*(-(P2*q2) + sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel17(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + (2*k2 + q2)*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(q2)*(-(P2*q2) + sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel18(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*(2*k2*P2*sqrt(P2*q2)*zk^2*zq*(2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) - 4*sqrt(k2*P2)*zk*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2) + (P2*q2 - P2*q2*zq^2)*(-2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + k2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + q2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + 2*k2*P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*((3*(A2in*B1in + A1in*B2in)*P2 + 4*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 - P2*q2*zq^2) + 2*sqrt(k2*P2)*zk*((-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)))))/(9. *P2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel21(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((k2^2*P2^2*(-8*B1in*B2in + 2*A1in*A2in*(P2 - 4*q2))*zk^4 + k2*P2*zk^2*((-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*(P2*q2 + 2*P2*q2*zq^2) + 4*k2*(-(A1in*A2in*P2^2) + P2*(4*B1in*B2in + 5*A1in*A2in*q2) + 2*A1in*A2in*P2*q2*zq^2) - 4*sqrt(k2)*sqrt(q2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 16*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + 4*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(4*B1in*B2in + A1in*A2in*(-P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*(2*k2^2*(A1in*A2in*P2^2 - 2*P2*(2*B1in*B2in + 3*A1in*A2in*q2) + 2*A1in*A2in*P2*q2*zq^2) + k2*((4*B1in*B2in - A1in*A2in*(P2 + 8*q2))*(P2*q2 - P2*q2*zq^2) + 4*sqrt(k2)*sqrt(q2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in + A1in*A2in*q2) - 4*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 28*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2 - 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + 2*P2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(k2*(-8*B1in*B2in + 2*A1in*A2in*(P2 - 10*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(4*B1in*B2in + A1in*A2in*(-P2 - 8*q2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel22(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*k2^2*P2^2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*zk^4*(P2*q2 + 2*P2*q2*zq^2) - 2*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(P2*(4*B1in*B2in - A1in*A2in*(P2 - 36*q2))*q2 + 8*P2*(4*B1in*B2in - A1in*A2in*P2)*q2*zq^2 + 4*sqrt(k2)*sqrt(q2)*(-3*A1in*A2in*P2^2 + 4*P2*(3*B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*P2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(P2*q2*(28*B1in*B2in - A1in*A2in*(7*P2 + 20*q2)) + 4*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2 + sqrt(k2)*sqrt(q2)*(15*A1in*A2in*P2^2 + P2*(-60*B1in*B2in + 4*A1in*A2in*q2) - 64*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*((4*B1in*B2in - A1in*A2in*(P2 - 36*q2))*(P2*q2 - P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(-21*A1in*A2in*P2^2 + 4*P2*(21*B1in*B2in + A1in*A2in*q2) + 80*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*(k2^2*(P2*q2 - P2*q2*zq^2)*(-7*A1in*A2in*P2^2 + 4*P2*(7*B1in*B2in + 3*A1in*A2in*q2) + 16*A1in*A2in*P2*q2*zq^2) - 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(P2*q2*(28*B1in*B2in - A1in*A2in*(7*P2 + 20*q2)) + 4*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2 + 4*sqrt(k2)*sqrt(q2)*(3*A1in*A2in*P2^2 + 4*P2*(-3*B1in*B2in + A1in*A2in*q2) - 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - k2*((-(P2*q2) + P2*q2*zq^2)*(P2*q2*(28*B1in*B2in - A1in*A2in*(7*P2 + 20*q2)) + 4*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2) - 8*sqrt(k2)*sqrt(q2)*(-(P2*q2) + P2*q2*zq^2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 7*k2*P2*q2*(-3*A1in*A2in*P2^2 + 4*P2*(3*B1in*B2in - A1in*A2in*q2) + 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) + k2*P2*zk^2*((P2*q2 + 2*P2*q2*zq^2)*(P2*q2*(-28*B1in*B2in + A1in*A2in*(7*P2 + 20*q2)) - 4*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2) + 4*k2*P2*q2*(-3*A1in*A2in*P2^2 + 4*P2*(3*B1in*B2in - A1in*A2in*q2) + 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-(A1in*A2in*P2^3*q2) + 16*P2^2*q2*(2*B1in*B2in + A1in*A2in*q2)*zq^2 + 16*A1in*A2in*P2^2*q2^2*zq^4 + 4*P2^2*(B1in*B2in*q2 + A1in*A2in*(q2^2 - 2*P2*q2*zq^2))) + k2*(11*A1in*A2in*P2^3*q2 - 8*P2^2*q2*(5*B1in*B2in + 3*A1in*A2in*q2)*zq^2 - 32*A1in*A2in*P2^2*q2^2*zq^4 + 2*P2^2*(-22*B1in*B2in*q2 + A1in*A2in*(-14*q2^2 + 5*P2*q2*zq^2)))))/(18. *P2*(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel23(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((sqrt(P2*q2)*zq*(-2*k2^2*P2^2*sqrt(P2*q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*zk^4*zq - k2*P2*sqrt(P2*q2)*zk^2*zq*((4*B1in*B2in + A1in*A2in*(P2 - 8*q2))*(P2*q2 + 2*P2*q2*zq^2) - 2*k2*(A1in*A2in*P2^2 + P2*(4*B1in*B2in + 6*A1in*A2in*q2) + 4*A1in*A2in*P2*q2*zq^2) + 8*sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) + 4*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 16*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + P2*sqrt(k2*P2)*zk*(-(k2*(P2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(4*B1in*B2in + A1in*A2in*(P2 - 12*q2))*q2*zq^2 + 2*sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) + 20*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*P2*(4*B1in*B2in + A1in*A2in*(P2 - 8*q2))*q2*zq^2 + sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) + 32*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + (k2*P2)^1.5*zk^3*(4*A1in*A2in*P2*q2^2 + 2*sqrt(k2)*P2*(4*B1in*B2in + A1in*A2in*P2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 4*P2*q2*zq^2*(4*B1in*B2in + A1in*A2in*(P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*sqrt(P2*q2)*zq*(4*A1in*A2in*k2^2*(-(P2*q2) + P2*q2*zq^2) + k2*((4*B1in*B2in + A1in*A2in*(P2 - 8*q2))*(P2*q2 - P2*q2*zq^2) + 2*sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 28*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) - 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(4*B1in*B2in + A1in*A2in*(P2 - 8*q2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel24(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*A1in*A2in*P2*(-4*k2^2*P2^2*q2*zk^4 + 2*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(3*q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*P2*zk^2*(P2*q2^2 - 4*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 - 2*P2*q2*zq^2*(k2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + q2*(2*P2*q2*zq^2 + P2*(9*k2 - 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + 2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(3*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-q2 + 5*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(3*P2*q2*zq^2 - P2*(3*q2 + 7*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*(5*k2^2*(-(P2*q2) + P2*q2*zq^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(q2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(-(P2*q2^2) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-8*P2*q2*zq^2 + 7*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + q2*(P2*q2*zq^2 + 8*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel25(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(2*(-(A2in*B1in) + A1in*B2in)*k2^2*P2^2*sqrt(P2*q2)*zk^4*zq + (k2*P2)^1.5*zk^3*(4*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in - A1in*B2in)*P2 - 4*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(k2)*sqrt(q2)*((A2in*B1in - A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*P2*zk^2*(-((-3*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 + 2*P2*q2*zq^2)) + k2*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + P2*((A2in*B1in + A1in*B2in)*q2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)) - 8*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*((A2in*B1in - A1in*B2in)*P2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 4*(A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + P2*sqrt(k2*P2)*zk*(k2*((-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + q2*((-(A2in*B1in) + A1in*B2in)*P2 + 4*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 - 5*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*sqrt(P2*q2)*zq*(-3*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((A2in*B1in - A1in*B2in)*P2 + 8*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*((A2in*B1in + A1in*B2in)*k2^2*(-(P2*q2) + P2*q2*zq^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(3*(A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq - 4*(A2in*B1in + A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(-((-3*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2)) + 2*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*((A2in*B1in - A1in*B2in)*P2 - 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 7*(A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel26(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(4*k2^2*P2^2*q2*zk^4 - 2*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(3*q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(3*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(q2 - 5*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(-3*P2*q2*zq^2 + P2*(3*q2 + 7*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + k2*P2*zk^2*(-(P2*q2^2) + 4*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + 2*P2*q2*zq^2*(k2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + q2*(-2*P2*q2*zq^2 + P2*(-9*k2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*(-5*k2^2*(-(P2*q2) + P2*q2*zq^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(-q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(P2*q2^2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(8*P2*q2*zq^2 - 7*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - q2*(P2*q2*zq^2 + 8*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel27(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(4*k2^2*P2^2*q2*zk^4 - 2*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(3*q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(3*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(q2 - 5*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(-3*P2*q2*zq^2 + P2*(3*q2 + 7*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + k2*P2*zk^2*(-(P2*q2^2) + 4*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + 2*P2*q2*zq^2*(k2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + q2*(-2*P2*q2*zq^2 + P2*(-9*k2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*(-5*k2^2*(-(P2*q2) + P2*q2*zq^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(-q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(P2*q2^2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(8*P2*q2*zq^2 - 7*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - q2*(P2*q2*zq^2 + 8*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel28(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*k2^2*P2^2*sqrt(P2*q2)*zk^4*zq*(2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*(k2*P2)^1.5*zk^3*(2*(A2in*B1in - A1in*B2in)*P2*q2^2 - 4*(A2in*B1in + A1in*B2in)*(P2*q2)^1.5*zq^3 + q2*sqrt(P2*q2)*zq*((A2in*B1in + A1in*B2in)*P2 + 4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + 4*sqrt(k2)*sqrt(q2)*((A2in*B1in - A1in*B2in)*P2*q2*zq^2 + P2*((A2in*B1in - A1in*B2in)*q2 - (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*P2*zk^2*((P2*q2 + 2*P2*q2*zq^2)*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 - q2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + k2*(-6*(A2in*B1in + A1in*B2in)*P2^2*q2*zq^2 + 4*(A2in*B1in - A1in*B2in)*(P2*q2)^1.5*zq^3 - P2*q2*((A2in*B1in + A1in*B2in)*P2 + 10*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + 8*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*(2*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(4*(-(A2in*B1in) + A1in*B2in)*q2 + 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 4*k2*P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + P2*(k2^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 - P2*q2*zq^2) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 - q2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + 4*sqrt(k2)*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*((-(P2*q2) + P2*q2*zq^2)*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 - q2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) - 8*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*(-(P2*q2) + P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 7*k2*P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) + 2*P2*sqrt(k2*P2)*zk*(2*(-(A2in*B1in) + A1in*B2in)*k2*P2*q2^2 + sqrt(P2*q2)*zq*((A2in*B1in + A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(7*k2 - 15*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(5*k2 - 12*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + (A2in*B1in + A1in*B2in)*P2*q2*zq^2*(k2 - 6*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) - q2*((A2in*B1in + A1in*B2in)*P2*sqrt(P2*q2)*zq*(k2 - 9*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*(A2in*B1in - A1in*B2in)*P2*q2*zq^2*(k2 - 3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*k2 - 3*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(9. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel31(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*((k2*P2)^1.5*(-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*zk^3 + sqrt(k2*P2)*zk*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2 + k2*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in + A1in*A2in*q2) + 4*A1in*A2in*P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 8*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + 2*k2*P2*sqrt(P2*q2)*zk^2*zq*(4*B1in*B2in + A1in*A2in*(-P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*sqrt(P2*q2)*zq*(k2*(-4*B1in*B2in + A1in*A2in*(P2 - 4*q2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(4*B1in*B2in + A1in*A2in*(-P2 - 8*q2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel32(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*((k2*P2)^1.5*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*zk^3*(P2*q2 + 2*P2*q2*zq^2) + P2*sqrt(P2*q2)*zq*(k2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*(P2*q2 - P2*q2*zq^2) + 2*sqrt(k2)*sqrt(q2)*(4*B1in*B2in - A1in*A2in*(4*k2 + P2 + 8*q2))*(P2*q2 - P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*q2*(3*A1in*A2in*P2^2 + 4*P2*(-3*B1in*B2in + 5*A1in*A2in*q2) - 32*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + sqrt(k2*P2)*zk*(k2*P2*q2*(-3*A1in*A2in*P2^2 + 4*P2*(3*B1in*B2in - A1in*A2in*q2) + 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 - (-(P2*q2) + P2*q2*zq^2)*(2*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*q2))*zq^2 + k2*(A1in*A2in*P2^2 - 4*P2*(B1in*B2in + A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)) + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(A1in*A2in*P2^3*q2 + 8*P2^2*q2*(2*B1in*B2in - A1in*A2in*q2)*zq^2 + 16*A1in*A2in*P2^2*q2^2*zq^4 - 4*P2^2*(B1in*B2in*q2 + A1in*A2in*(-q2^2 + P2*q2*zq^2)))) - k2*P2*sqrt(P2*q2)*zk^2*zq*(12*A1in*A2in*P2*q2^2 - 6*sqrt(k2)*P2*(-4*B1in*B2in + A1in*A2in*P2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 4*P2*q2*zq^2*(4*B1in*B2in + A1in*A2in*(-P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(3. *P2^2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel33(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-0.3333333333333333*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(4*A1in*A2in*P2*q2^2 + 2*k2*P2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*zk^2 + k2*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2) - 2*sqrt(k2)*P2*(4*B1in*B2in + A1in*A2in*P2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(4*B1in*B2in + A1in*A2in*(P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 2*P2*q2*zq^2*(4*B1in*B2in + A1in*A2in*(P2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + q2*(-16*A1in*A2in*P2*q2*zq^2 + P2*(4*B1in*B2in + A1in*A2in*(P2 - 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel34(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*A1in*A2in*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel35(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*A2in*B1in + 2*A1in*B2in)*k2*P2*zk^2 + 2*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + k2*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + q2*((-(A2in*B1in) + A1in*B2in)*P2 + 6*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 8*(A2in*B1in + A1in*B2in)*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 4*sqrt(k2*P2)*zk*(-((A2in*B1in + A1in*B2in)*q2) + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq + (A2in*B1in + A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel36(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((16*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel37(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel38(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(k2*P2*zk^2*(2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + ((A2in*B1in - A1in*B2in)*k2 + (A2in*B1in - A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(4*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(2*(A2in*B1in - A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2*P2)*zk*(-2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*(-((A2in*B1in + A1in*B2in)*P2) + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel41(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*A1in*A2in*(sqrt(k2*P2)*q2*zk*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - sqrt(k2)*sqrt(q2)*(P2*q2 + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - P2*(q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel42(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-4*A1in*A2in*(-(P2*q2) + P2*q2*zq^2)*(sqrt(k2*P2)*q2*zk*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - sqrt(k2)*sqrt(q2)*(P2*q2 + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - P2*(q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *P2*(k2*P2 - k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel43(y) = 0
                
                kernel44(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(sqrt(k2*P2)*q2*zk*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - sqrt(k2)*sqrt(q2)*(P2*q2 + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - P2*(q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel45(y) = 0
                
                kernel46(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-4*sqrt(P2*q2)*zq*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(sqrt(k2*P2)*q2*zk*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - sqrt(k2)*sqrt(q2)*(P2*q2 + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - P2*(q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *P2*(k2*P2 - k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel47(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*(A2in*B1in + A1in*B2in)*(-(P2*q2) + P2*q2*zq^2)*(sqrt(k2*P2)*q2*zk*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - sqrt(k2)*sqrt(q2)*(P2*q2 + 4*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk + sqrt(P2*q2)*zq) - P2*(q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *P2*(k2*P2 - k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel48(y) = 0
                
                kernel51(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(A2in*B1in + A1in*B2in)*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(-(k2*P2) + k2*P2*zk^2))
                
                kernel52(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*(A2in*B1in + A1in*B2in)*(P2*q2 - P2*q2*zq^2)*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(P2*(k2*P2 - k2*P2*zk^2)))
                
                kernel53(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*sqrt(P2*q2)*zq*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(-(k2*P2) + k2*P2*zk^2))
                
                kernel54(y) = 0
                
                kernel55(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(-(k2*P2) + k2*P2*zk^2))
                
                kernel56(y) = 0
                
                kernel57(y) = 0
                
                kernel58(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*A1in*A2in*(P2*q2 - P2*q2*zq^2)*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(-(k2*P2) + k2*P2*zk^2))
                
                kernel61(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*((A2in*B1in + A1in*B2in)*(k2*P2)^1.5*zk^3 - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + (A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - k2*P2*zk^2*((-(A2in*B1in) + A1in*B2in)*q2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*((-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + (A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2*P2)*zk*(k2*(-((A2in*B1in + A1in*B2in)*P2) + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(P2*q2)*zq*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel62(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-0.3333333333333333*((k2*P2)^1.5*zk^3*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 6*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + k2*(P2*q2 - P2*q2*zq^2)*(4*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(2*(A2in*B1in - A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) + 2*sqrt(k2)*P2*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*k2 + (-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 - P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(4*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 + P2*(2*(A2in*B1in - A1in*B2in)*q2 - 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*P2*zk^2*(-((2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + 4*P2*q2*zq^2)) + 2*sqrt(k2)*sqrt(q2)*(4*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 + P2*(2*(A2in*B1in - A1in*B2in)*q2 - 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2*P2)*zk*((-(P2*q2) + P2*q2*zq^2)*(k2*((A2in*B1in + A1in*B2in)*P2 + 4*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(P2*q2)*zq*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) + 2*sqrt(k2)*sqrt(q2)*(-(P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)) + 4*P2*q2*zq^2*((A2in*B1in + A1in*B2in)*P2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel63(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel64(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(k2*P2)^1.5*zk^3*((A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2*P2)*zk*((A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(k2*P2 - 2*P2*q2*zq^2) + q2*(-2*(A2in*B1in + A1in*B2in)*k2*P2 + (A2in*B1in - A1in*B2in)*P2*sqrt(P2*q2)*zq + 2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2) + 4*sqrt(k2)*P2*sqrt(q2)*((A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 2*(A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + k2*P2*zk^2*(-(q2*((A2in*B1in - A1in*B2in)*P2 + 4*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) + 2*(A2in*B1in - A1in*B2in)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*(k2*(-((A2in*B1in - A1in*B2in)*(-(P2*q2) + P2*q2*zq^2)) + sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-(q2*((A2in*B1in - A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) + (A2in*B1in - A1in*B2in)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel65(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*A1in*A2in*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel66(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((sqrt(P2*q2)*zq*(-2*(k2*P2)^1.5*(4*B1in*B2in + A1in*A2in*P2)*sqrt(P2*q2)*zk^3*zq + k2*P2*zk^2*(-4*A1in*A2in*P2*q2^2 + 2*(4*B1in*B2in + A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - P2*q2*(4*B1in*B2in + A1in*A2in*(P2 - 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(k2*P2*(4*B1in*B2in + A1in*A2in*(P2 - 4*q2)) + P2*(4*B1in*B2in + A1in*A2in*P2)*q2 + 4*A1in*A2in*P2*q2^2 - 2*P2*(4*B1in*B2in + A1in*A2in*P2)*q2*zq^2 - 4*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(4*B1in*B2in + A1in*A2in*(P2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*(k2*((4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(P2*q2 - P2*q2*zq^2) - sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-4*A1in*A2in*P2*q2^2 + (4*B1in*B2in + A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - P2*q2*(4*B1in*B2in + A1in*A2in*(P2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel67(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-((k2*P2)^1.5*zk^3*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2)) + P2*sqrt(P2*q2)*zq*(k2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*(P2*q2 - P2*q2*zq^2) + 2*sqrt(k2)*(4*B1in*B2in + A1in*A2in*(4*k2 - P2))*sqrt(q2)*(-(P2*q2) + P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + sqrt(k2*P2)*zk*(-2*sqrt(k2)*P2*sqrt(q2)*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + (-(P2*q2) + P2*q2*zq^2)*(-4*A1in*A2in*k2*P2*q2 + (4*B1in*B2in - A1in*A2in*P2)*(k2*P2 - 2*P2*q2*zq^2))) + k2*P2*sqrt(P2*q2)*zk^2*zq*(4*A1in*A2in*P2*q2^2 + 2*(4*B1in*B2in - A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*q2*(-12*B1in*B2in + A1in*A2in*(3*P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *P2*sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel68(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(6. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel71(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(A2in*B1in + A1in*B2in)*k2^2*P2^2*zk^4 - 2*(k2*P2)^1.5*zk^3*((-(A2in*B1in) + A1in*B2in)*q2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*sqrt(k2*P2)*zk*(-(sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*((A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + (A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + k2*((-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*((-(A2in*B1in) + A1in*B2in)*q2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + k2*P2*zk^2*(k2*(-3*(A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + 2*(sqrt(P2*q2)*zq*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*(k2^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + (A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(2*(A2in*B1in - A1in*B2in)*q2*sqrt(P2*q2)*zq - (A2in*B1in + A1in*B2in)*(P2*q2*zq^2 + 2*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel72(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*k2^2*P2^2*zk^4*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 6*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + 2*(k2*P2)^1.5*zk^3*(-((2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + 4*P2*q2*zq^2)) + 2*sqrt(k2)*sqrt(q2)*(4*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 + P2*(2*(A2in*B1in - A1in*B2in)*q2 - 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*(k2^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 - P2*q2*zq^2) + k2*P2*q2*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 6*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*(-((-(P2*q2) + P2*q2*zq^2)*(2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 - q2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq))) + 4*(A2in*B1in + A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(P2*q2 - P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 3*k2*P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) + 2*sqrt(k2*P2)*zk*(2*(A2in*B1in - A1in*B2in)*k2*P2^2*q2^2 + sqrt(P2*q2)*zq*(4*(A2in*B1in - A1in*B2in)*k2*(P2*q2)^1.5*zq^3 + (A2in*B1in + A1in*B2in)*P2^2*q2*zq^2*(k2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*sqrt(P2*q2)*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(k2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 3*(A2in*B1in + A1in*B2in)*sqrt(k2)*P2^2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(k2 - sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) - P2*q2*(6*(A2in*B1in - A1in*B2in)*P2*q2*zq^2*(k2 - sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*(A2in*B1in - A1in*B2in)*sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*k2 - sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + (A2in*B1in + A1in*B2in)*P2*sqrt(P2*q2)*zq*(k2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + k2*P2*zk^2*(4*(A2in*B1in + A1in*B2in)*P2^2*q2^2*zq^4 + P2*q2^2*(3*(A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + 16*(A2in*B1in + A1in*B2in)*sqrt(k2)*P2^2*q2^1.5*zq^2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 6*(A2in*B1in + A1in*B2in)*k2*P2^2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + 12*(-(A2in*B1in) + A1in*B2in)*k2*P2*q2*sqrt(P2*q2)*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 - 8*(A2in*B1in - A1in*B2in)*(P2*q2)^1.5*zq^3*(k2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - q2*(4*(A2in*B1in + A1in*B2in)*P2^2*q2*zq^2 + 4*(A2in*B1in - A1in*B2in)*(P2*q2)^1.5*zq^3 - 2*(A2in*B1in - A1in*B2in)*P2*sqrt(P2*q2)*zq*(7*k2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + (A2in*B1in + A1in*B2in)*P2^2*(3*k2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(6. *P2*(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel73(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(P2*(-k2 + q2) + 2*sqrt(k2*P2)*zk*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq))*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel74(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*k2^2*P2^2*zk^4*(-((A2in*B1in + A1in*B2in)*q2) + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + k2*P2*zk^2*(k2*P2*(3*(A2in*B1in + A1in*B2in)*q2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + (-((A2in*B1in + A1in*B2in)*q2) + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + 2*P2*q2*zq^2) - 4*sqrt(k2)*P2*sqrt(q2)*((A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 2*(A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + (k2*P2)^1.5*zk^3*(q2*((A2in*B1in - A1in*B2in)*P2 + 4*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) - 2*(A2in*B1in - A1in*B2in)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) - P2*(-((A2in*B1in + A1in*B2in)*k2^2*(-(P2*q2) + P2*q2*zq^2)) + k2*P2*q2*((A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*((-((A2in*B1in + A1in*B2in)*q2) + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2) + 2*sqrt(k2)*P2*sqrt(q2)*(-2*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + (A2in*B1in + A1in*B2in)*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) + P2*sqrt(k2*P2)*zk*(k2*(3*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 - q2*((A2in*B1in - A1in*B2in)*P2 + 4*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(k2)*sqrt(q2)*((A2in*B1in - A1in*B2in)*P2 - (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*(A2in*B1in + A1in*B2in)*q2*sqrt(P2*q2)*zq - (A2in*B1in - A1in*B2in)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel75(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-0.3333333333333333*(A1in*A2in*(P2*(-k2 + q2) + 2*sqrt(k2*P2)*zk*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq))*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/((-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel76(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((sqrt(P2*q2)*zq*(2*k2^2*P2^2*(4*B1in*B2in + A1in*A2in*P2)*sqrt(P2*q2)*zk^4*zq + k2*P2*sqrt(P2*q2)*zk^2*zq*(-2*k2*P2*(4*B1in*B2in + A1in*A2in*(P2 - 2*q2)) + (4*B1in*B2in + A1in*A2in*P2)*(-(P2*q2) + 2*P2*q2*zq^2) + 4*sqrt(k2)*P2*(4*B1in*B2in + A1in*A2in*P2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 8*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + P2*sqrt(P2*q2)*zq*(4*A1in*A2in*k2^2*(-(P2*q2) + P2*q2*zq^2) + k2*P2*(4*B1in*B2in + A1in*A2in*P2)*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*((4*B1in*B2in + A1in*A2in*P2)*(P2*q2 - P2*q2*zq^2) - 2*sqrt(k2)*P2*(4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 4*A1in*A2in*k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) + (k2*P2)^1.5*zk^3*(4*A1in*A2in*P2*q2^2 - 2*(4*B1in*B2in + A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*q2*(4*B1in*B2in + A1in*A2in*(P2 - 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + P2*sqrt(k2*P2)*zk*(k2*(-(P2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))) + P2*(12*B1in*B2in + A1in*A2in*(3*P2 - 4*q2))*q2*zq^2 + 2*sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 + 4*P2*(B1in*B2in + A1in*A2in*q2) - 4*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-4*A1in*A2in*sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - (4*B1in*B2in + A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(3. *P2*(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel77(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*k2^2*P2^2*zk^4*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2) + k2*P2*zk^2*((-(P2*q2) + 2*P2*q2*zq^2)*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2) + k2*P2*(-3*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + 4*P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 2*q2))*zq^2) + 4*sqrt(k2)*P2*sqrt(q2)*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - 2*k2*P2*q2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + P2*(k2^2*(-(P2*q2) + P2*q2*zq^2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in - A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2) + k2*P2*q2*(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + P2*(8*B1in*B2in - 2*A1in*A2in*P2)*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*((-(P2*q2) + P2*q2*zq^2)*(-(P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2))) + P2*(-8*B1in*B2in + 2*A1in*A2in*P2)*q2*zq^2) - 4*sqrt(k2)*P2*sqrt(q2)*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(P2*q2 - P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)) - 2*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq*(4*A1in*A2in*P2*q2^2 + 2*(4*B1in*B2in - A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*q2*(-12*B1in*B2in + A1in*A2in*(3*P2 + 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) - 2*P2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(k2*((-12*B1in*B2in + A1in*A2in*(3*P2 + 4*q2))*(-(P2*q2) + P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(A1in*A2in*P2^2 - 4*P2*(B1in*B2in + 3*A1in*A2in*q2) + 8*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(4*A1in*A2in*P2*q2^2 + (4*B1in*B2in - A1in*A2in*P2)*(2*P2*q2*zq^2 + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*q2*(-4*B1in*B2in + A1in*A2in*(P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(6. *P2*(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel78(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(P2*(-k2 + q2) + 2*sqrt(k2*P2)*zk*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq))*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(12. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel81(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(4*(A2in*B1in - A1in*B2in)*(k2*P2)^1.5*q2*zk^3 - k2^2*P2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) - 2*(A2in*B1in + A1in*B2in)*k2*P2^2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 - k2*P2*zk^2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 4*sqrt(k2)*P2*sqrt(k2*P2)*sqrt(q2)*zk*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*((A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + (A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(k2*P2*zk^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + 4*(A2in*B1in - A1in*B2in)*sqrt(k2*P2)*zk*(-(P2*q2) + P2*q2*zq^2) + P2*(-2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(k2)*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel82(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-2*(-(k2^2*P2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 - P2*q2*zq^2)) - 2*(k2*P2)^1.5*q2*zk^3*(4*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(2*(-(A2in*B1in) + A1in*B2in)*q2 + 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) + k2*P2^2*q2*((A2in*B1in + A1in*B2in)*P2*q2 + 6*(-(A2in*B1in) + A1in*B2in)*q2*sqrt(P2*q2)*zq + 2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + k2*P2*zk^2*(q2*(-(P2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)) + 4*P2*q2*zq^2*((A2in*B1in + A1in*B2in)*P2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)) + 4*sqrt(k2)*sqrt(q2)*(2*P2*q2*zq^2*((A2in*B1in + A1in*B2in)*P2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + P2*q2*((A2in*B1in + A1in*B2in)*P2 + 4*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*sqrt(k2)*P2*sqrt(k2*P2)*sqrt(q2)*zk*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*(A2in*B1in + A1in*B2in)*(P2*q2)^1.5*zq^3 + q2*sqrt(P2*q2)*zq*((A2in*B1in + A1in*B2in)*P2 + 6*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*(4*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(2*(-(A2in*B1in) + A1in*B2in)*q2 + 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(k2*P2*zk^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(P2*q2 + 2*P2*q2*zq^2) + P2*(P2*q2 - P2*q2*zq^2)*(-2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)) + 2*sqrt(k2*P2)*zk*(P2*q2 - P2*q2*zq^2)*(4*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + P2*(2*(-(A2in*B1in) + A1in*B2in)*q2 + 3*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)) - 2*sqrt(k2)*P2*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(2*P2*q2 + 3*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - 2*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 3*k2*P2^2*q2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2)))/(3. *P2*(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel83(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel84(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-4*P2*((A2in*B1in + A1in*B2in)*k2 - (A2in*B1in + A1in*B2in)*q2 - (A2in*B1in - A1in*B2in)*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq))*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel85(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*A1in*A2in*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel86(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*sqrt(P2*q2)*zq*(sqrt(k2*P2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*zk - (4*B1in*B2in + A1in*A2in*(4*k2 + P2))*sqrt(P2*q2)*zq)*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel87(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(-4*A1in*A2in*P2*q2^2 + 2*(4*B1in*B2in - A1in*A2in*P2)*sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk - sqrt(P2*q2)*zq) + q2*(4*B1in*B2in*P2 - A1in*A2in*P2^2 + 8*A1in*A2in*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + k2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2))*(k2*P2*q2*zk^2 + k2*(-(P2*q2) + P2*q2*zq^2) - 2*sqrt(k2)*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + k2*P2*q2*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel88(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)^2*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))

                kernel[i,s,j,m,1,1]=allweight*gausslegendreint64(kernel11)

                kernel[i,s,j,m,1,2]=allweight*gausslegendreint64(kernel12)
               
                kernel[i,s,j,m,1,3]=allweight*gausslegendreint64(kernel13)
               
                kernel[i,s,j,m,1,4]=allweight*gausslegendreint64(kernel14)
               
                kernel[i,s,j,m,1,5]=allweight*gausslegendreint64(kernel15)
               
                kernel[i,s,j,m,1,6]=allweight*gausslegendreint64(kernel16)
               
                kernel[i,s,j,m,1,7]=allweight*gausslegendreint64(kernel17)
               
                kernel[i,s,j,m,1,8]=allweight*gausslegendreint64(kernel18)
               
                kernel[i,s,j,m,2,1]=allweight*gausslegendreint64(kernel21)
               
                kernel[i,s,j,m,2,2]=allweight*gausslegendreint64(kernel22)
               
                kernel[i,s,j,m,2,3]=allweight*gausslegendreint64(kernel23)
               
                kernel[i,s,j,m,2,4]=allweight*gausslegendreint64(kernel24)
               
                kernel[i,s,j,m,2,5]=allweight*gausslegendreint64(kernel25)
               
                kernel[i,s,j,m,2,6]=allweight*gausslegendreint64(kernel26)
               
                kernel[i,s,j,m,2,7]=allweight*gausslegendreint64(kernel27)
               
                kernel[i,s,j,m,2,8]=allweight*gausslegendreint64(kernel28)
               
                kernel[i,s,j,m,3,1]=allweight*gausslegendreint64(kernel31)
               
                kernel[i,s,j,m,3,2]=allweight*gausslegendreint64(kernel32)
               
                kernel[i,s,j,m,3,3]=allweight*gausslegendreint64(kernel33)
               
                kernel[i,s,j,m,3,4]=allweight*gausslegendreint64(kernel34)
               
                kernel[i,s,j,m,3,5]=allweight*gausslegendreint64(kernel35)
               
                kernel[i,s,j,m,3,6]=allweight*gausslegendreint64(kernel36)
               
                kernel[i,s,j,m,3,7]=allweight*gausslegendreint64(kernel37)
               
                kernel[i,s,j,m,3,8]=allweight*gausslegendreint64(kernel38)
               
                kernel[i,s,j,m,4,1]=allweight*gausslegendreint64(kernel41)
               
                kernel[i,s,j,m,4,2]=allweight*gausslegendreint64(kernel42)
               
                kernel[i,s,j,m,4,3]=allweight*gausslegendreint64(kernel43)
               
                kernel[i,s,j,m,4,4]=allweight*gausslegendreint64(kernel44)
               
                kernel[i,s,j,m,4,5]=allweight*gausslegendreint64(kernel45)
               
                kernel[i,s,j,m,4,6]=allweight*gausslegendreint64(kernel46)
               
                kernel[i,s,j,m,4,7]=allweight*gausslegendreint64(kernel47)
               
                kernel[i,s,j,m,4,8]=allweight*gausslegendreint64(kernel48)
               
                kernel[i,s,j,m,5,1]=allweight*gausslegendreint64(kernel51)
               
                kernel[i,s,j,m,5,2]=allweight*gausslegendreint64(kernel52)
               
                kernel[i,s,j,m,5,3]=allweight*gausslegendreint64(kernel53)
               
                kernel[i,s,j,m,5,4]=allweight*gausslegendreint64(kernel54)
               
                kernel[i,s,j,m,5,5]=allweight*gausslegendreint64(kernel55)
               
                kernel[i,s,j,m,5,6]=allweight*gausslegendreint64(kernel56)
               
                kernel[i,s,j,m,5,7]=allweight*gausslegendreint64(kernel57)
               
                kernel[i,s,j,m,5,8]=allweight*gausslegendreint64(kernel58)
               
                kernel[i,s,j,m,6,1]=allweight*gausslegendreint64(kernel61)
               
                kernel[i,s,j,m,6,2]=allweight*gausslegendreint64(kernel62)
               
                kernel[i,s,j,m,6,3]=allweight*gausslegendreint64(kernel63)
               
                kernel[i,s,j,m,6,4]=allweight*gausslegendreint64(kernel64)
               
                kernel[i,s,j,m,6,5]=allweight*gausslegendreint64(kernel65)
               
                kernel[i,s,j,m,6,6]=allweight*gausslegendreint64(kernel66)
               
                kernel[i,s,j,m,6,7]=allweight*gausslegendreint64(kernel67)
               
                kernel[i,s,j,m,6,8]=allweight*gausslegendreint64(kernel68)
               
                kernel[i,s,j,m,7,1]=allweight*gausslegendreint64(kernel71)
               
                kernel[i,s,j,m,7,2]=allweight*gausslegendreint64(kernel72)
               
                kernel[i,s,j,m,7,3]=allweight*gausslegendreint64(kernel73)
               
                kernel[i,s,j,m,7,4]=allweight*gausslegendreint64(kernel74)
               
                kernel[i,s,j,m,7,5]=allweight*gausslegendreint64(kernel75)
               
                kernel[i,s,j,m,7,6]=allweight*gausslegendreint64(kernel76)
               
                kernel[i,s,j,m,7,7]=allweight*gausslegendreint64(kernel77)
               
                kernel[i,s,j,m,7,8]=allweight*gausslegendreint64(kernel78)
               
                kernel[i,s,j,m,8,1]=allweight*gausslegendreint64(kernel81)
               
                kernel[i,s,j,m,8,2]=allweight*gausslegendreint64(kernel82)
               
                kernel[i,s,j,m,8,3]=allweight*gausslegendreint64(kernel83)
               
                kernel[i,s,j,m,8,4]=allweight*gausslegendreint64(kernel84)
               
                kernel[i,s,j,m,8,5]=allweight*gausslegendreint64(kernel85)
               
                kernel[i,s,j,m,8,6]=allweight*gausslegendreint64(kernel86)
               
                kernel[i,s,j,m,8,7]=allweight*gausslegendreint64(kernel87)
               
                kernel[i,s,j,m,8,8]=allweight*gausslegendreint64(kernel88) 
            end
        end
    end
end

# 加一个尾巴
if dataset["mesonBSE"]["tailed"]==1
Threads.@threads for i = 1:length(meshk)
    k2 = meshk[i]
    for s = 1:length(meshz)
        zk = meshz[s]
        for j = 1:length(meshk)
            q2 = meshk[j]
            for m = 1:length(meshz)
                zq = meshz[m]
                weightzq = weightz[m]
                weightq2 = weightk[j]
                A1in = A1[j,m]
                A2in = A2[j,m]
                B1in = B1[j,m]
                B2in = B2[j,m]
                allweight = -weightzq*weightq2*q2/(16*pi^3)/((P2/4+q2+sqrt(P2*q2)*zq)*A1in^2+B1in^2)/((P2/4+q2-sqrt(P2*q2)*zq)*A2in^2+B2in^2)*z2^2*epsilon

                矢量介子的kernel尾巴如下所示

                kernel11(y) = 0
                
                kernel12(y) = 0
                
                kernel13(y) = 0
                
                kernel14(y) = 0
                
                kernel15(y) = 0
                
                kernel16(y) = 0
                
                kernel17(y) = 0
                
                kernel18(y) = 0
                
                kernel21(y) = 0
                
                kernel22(y) = 0
                
                kernel23(y) = 0
                
                kernel24(y) = 0
                
                kernel25(y) = 0
                
                kernel26(y) = 0
                
                kernel27(y) = 0
                
                kernel28(y) = 0
                
                kernel31(y) = 0
                
                kernel32(y) = 0
                
                kernel33(y) = 0
                
                kernel34(y) = 0
                
                kernel35(y) = 0
                
                kernel36(y) = 0
                
                kernel37(y) = 0
                
                kernel38(y) = 0
                
                kernel41(y) = 0
                
                kernel42(y) = 0
                
                kernel43(y) = 0
                
                kernel44(y) = 0
                
                kernel45(y) = 0
                
                kernel46(y) = 0
                
                kernel47(y) = 0
                
                kernel48(y) = 0
                
                kernel51(y) = 0
                
                kernel52(y) = 0
                
                kernel53(y) = 0
                
                kernel54(y) = 0
                
                kernel55(y) = 0
                
                kernel56(y) = 0
                
                kernel57(y) = 0
                
                kernel58(y) = 0
                
                kernel61(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-16*(A2in*B1in - A1in*B2in)*(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *(k2*P2)^1.5*zk*(-1 + zk^2)))
                
                kernel62(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((16*(A2in*B1in - A1in*B2in)*q2*(-1 + zq^2)*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *(k2*P2)^1.5*zk*(-1 + zk^2)))
                
                kernel63(y) = 0
                
                kernel64(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*P2*(A2in*B1in*(sqrt(k2*P2)*zk*zq*(sqrt(P2*q2) - 2*q2*zq) - sqrt(k2)*sqrt(q2)*(P2 - 2*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + A1in*B2in*(-(sqrt(k2*P2)*zk*zq*(sqrt(P2*q2) + 2*q2*zq)) + sqrt(k2)*sqrt(q2)*(P2 + 2*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))/(9. *(k2*P2)^1.5*zk*(-1 + zk^2)))
                
                kernel65(y) = 0
                
                kernel66(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*P2*sqrt(q2)*zq*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2 - 8*q2*zq^2))*(-(sqrt(k2*P2)*sqrt(q2)*zk*zq) + sqrt(k2)*sqrt(P2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *(k2*P2)^1.5*zk*(-1 + zk^2)))
                
                kernel67(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((64*A1in*A2in*P2*q2^1.5*zq*(-1 + zq^2)*(sqrt(k2*P2)*sqrt(q2)*zk*zq - sqrt(k2)*sqrt(P2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *(k2*P2)^1.5*zk*(-1 + zk^2)))
                
                kernel68(y) = 0
                
                kernel71(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*(A2in*B1in*(P2 - 2*sqrt(P2*q2)*zq) + A1in*B2in*(P2 + 2*sqrt(P2*q2)*zq)))/(9. *k2*P2*(-1 + zk^2)))
                
                kernel72(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-4*(A2in*B1in*(-6*sqrt(k2*P2)*sqrt(q2)*zk*zq*(sqrt(P2*q2) - 2*q2*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*q2*(P2 - 2*sqrt(P2*q2)*zq)*(-1 + zk^2 + zq^2 + 5*zk^2*zq^2 + 6*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + 3*y^2*(-1 + zk^2)*(-1 + zq^2))) + A1in*B2in*(-6*sqrt(k2*P2)*sqrt(q2)*zk*zq*(sqrt(P2*q2) + 2*q2*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*q2*(P2 + 2*sqrt(P2*q2)*zq)*(-1 + zk^2 + zq^2 + 5*zk^2*zq^2 + 6*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + 3*y^2*(-1 + zk^2)*(-1 + zq^2)))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel73(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*(A2in*B1in - A1in*B2in)*q2*zq*(-2*sqrt(k2*P2)*sqrt(q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*sqrt(P2*q2)*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*(-1 + zk^2)^2))
                
                kernel74(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-8*(A2in*B1in + A1in*B2in)*sqrt(q2)*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(1 - zk^2 + (-1 + 3*zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel75(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*A1in*A2in*(-2*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*q2*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel76(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((32*A1in*A2in*q2^1.5*zq^2*(2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(-1 + zk^2 + zq^2 - 3*zk^2*zq^2 - 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) - y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel77(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*sqrt(q2)*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2 - 8*q2*zq^2))*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(1 - zk^2 + (-1 + 3*zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel78(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2))*sqrt(q2)*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel81(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-16*(A2in*B1in*(P2 - 2*sqrt(P2*q2)*zq) + A1in*B2in*(P2 + 2*sqrt(P2*q2)*zq)))/(9. *k2*P2*(-1 + zk^2)))
                
                kernel82(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((16*(A2in*B1in*(-6*sqrt(k2*P2)*sqrt(q2)*zk*zq*(sqrt(P2*q2) - 2*q2*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*q2*(P2 - 2*sqrt(P2*q2)*zq)*(-1 + zk^2 + zq^2 + 5*zk^2*zq^2 + 6*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + 3*y^2*(-1 + zk^2)*(-1 + zq^2))) + A1in*B2in*(-6*sqrt(k2*P2)*sqrt(q2)*zk*zq*(sqrt(P2*q2) + 2*q2*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*q2*(P2 + 2*sqrt(P2*q2)*zq)*(-1 + zk^2 + zq^2 + 5*zk^2*zq^2 + 6*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + 3*y^2*(-1 + zk^2)*(-1 + zq^2)))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel83(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((32*(A2in*B1in - A1in*B2in)*q2*zq*(-2*sqrt(k2*P2)*sqrt(q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*sqrt(P2*q2)*(2*zk^2*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*(-1 + zk^2)^2))
                
                kernel84(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((32*(A2in*B1in + A1in*B2in)*sqrt(q2)*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel85(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-32*A1in*A2in*(-2*sqrt(k2*P2)*sqrt(q2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*q2*(2*zk^2*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel86(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((128*A1in*A2in*q2^1.5*zq^2*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel87(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((16*sqrt(q2)*(4*B1in*B2in - A1in*A2in*(P2 + 4*q2 - 8*q2*zq^2))*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(-1 + zk^2 + (1 + zk^2)*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                
                kernel88(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*(-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*sqrt(q2)*(-2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2)*P2*sqrt(q2)*(2*zk^2*zq^2 + 2*y*zk*sqrt(1 - zk^2)*zq*sqrt(1 - zq^2) + y^2*(-1 + zk^2)*(-1 + zq^2))))/(9. *k2^1.5*P2*(-1 + zk^2)^2))
                kernel[i,s,j,m,1,1]+=allweight*gausslegendreint64(kernel11)

                kernel[i,s,j,m,1,2]+=allweight*gausslegendreint64(kernel12)
               
                kernel[i,s,j,m,1,3]+=allweight*gausslegendreint64(kernel13)
               
                kernel[i,s,j,m,1,4]+=allweight*gausslegendreint64(kernel14)
               
                kernel[i,s,j,m,1,5]+=allweight*gausslegendreint64(kernel15)
               
                kernel[i,s,j,m,1,6]+=allweight*gausslegendreint64(kernel16)
               
                kernel[i,s,j,m,1,7]+=allweight*gausslegendreint64(kernel17)
               
                kernel[i,s,j,m,1,8]+=allweight*gausslegendreint64(kernel18)
               
                kernel[i,s,j,m,2,1]+=allweight*gausslegendreint64(kernel21)
               
                kernel[i,s,j,m,2,2]+=allweight*gausslegendreint64(kernel22)
               
                kernel[i,s,j,m,2,3]+=allweight*gausslegendreint64(kernel23)
               
                kernel[i,s,j,m,2,4]+=allweight*gausslegendreint64(kernel24)
               
                kernel[i,s,j,m,2,5]+=allweight*gausslegendreint64(kernel25)
               
                kernel[i,s,j,m,2,6]+=allweight*gausslegendreint64(kernel26)
               
                kernel[i,s,j,m,2,7]+=allweight*gausslegendreint64(kernel27)
               
                kernel[i,s,j,m,2,8]+=allweight*gausslegendreint64(kernel28)
               
                kernel[i,s,j,m,3,1]+=allweight*gausslegendreint64(kernel31)
               
                kernel[i,s,j,m,3,2]+=allweight*gausslegendreint64(kernel32)
               
                kernel[i,s,j,m,3,3]+=allweight*gausslegendreint64(kernel33)
               
                kernel[i,s,j,m,3,4]+=allweight*gausslegendreint64(kernel34)
               
                kernel[i,s,j,m,3,5]+=allweight*gausslegendreint64(kernel35)
               
                kernel[i,s,j,m,3,6]+=allweight*gausslegendreint64(kernel36)
               
                kernel[i,s,j,m,3,7]+=allweight*gausslegendreint64(kernel37)
               
                kernel[i,s,j,m,3,8]+=allweight*gausslegendreint64(kernel38)
               
                kernel[i,s,j,m,4,1]+=allweight*gausslegendreint64(kernel41)
               
                kernel[i,s,j,m,4,2]+=allweight*gausslegendreint64(kernel42)
               
                kernel[i,s,j,m,4,3]+=allweight*gausslegendreint64(kernel43)
               
                kernel[i,s,j,m,4,4]+=allweight*gausslegendreint64(kernel44)
               
                kernel[i,s,j,m,4,5]+=allweight*gausslegendreint64(kernel45)
               
                kernel[i,s,j,m,4,6]+=allweight*gausslegendreint64(kernel46)
               
                kernel[i,s,j,m,4,7]+=allweight*gausslegendreint64(kernel47)
               
                kernel[i,s,j,m,4,8]+=allweight*gausslegendreint64(kernel48)
               
                kernel[i,s,j,m,5,1]+=allweight*gausslegendreint64(kernel51)
               
                kernel[i,s,j,m,5,2]+=allweight*gausslegendreint64(kernel52)
               
                kernel[i,s,j,m,5,3]+=allweight*gausslegendreint64(kernel53)
               
                kernel[i,s,j,m,5,4]+=allweight*gausslegendreint64(kernel54)
               
                kernel[i,s,j,m,5,5]+=allweight*gausslegendreint64(kernel55)
               
                kernel[i,s,j,m,5,6]+=allweight*gausslegendreint64(kernel56)
               
                kernel[i,s,j,m,5,7]+=allweight*gausslegendreint64(kernel57)
               
                kernel[i,s,j,m,5,8]+=allweight*gausslegendreint64(kernel58)
               
                kernel[i,s,j,m,6,1]+=allweight*gausslegendreint64(kernel61)
               
                kernel[i,s,j,m,6,2]+=allweight*gausslegendreint64(kernel62)
               
                kernel[i,s,j,m,6,3]+=allweight*gausslegendreint64(kernel63)
               
                kernel[i,s,j,m,6,4]+=allweight*gausslegendreint64(kernel64)
               
                kernel[i,s,j,m,6,5]+=allweight*gausslegendreint64(kernel65)
               
                kernel[i,s,j,m,6,6]+=allweight*gausslegendreint64(kernel66)
               
                kernel[i,s,j,m,6,7]+=allweight*gausslegendreint64(kernel67)
               
                kernel[i,s,j,m,6,8]+=allweight*gausslegendreint64(kernel68)
               
                kernel[i,s,j,m,7,1]+=allweight*gausslegendreint64(kernel71)
               
                kernel[i,s,j,m,7,2]+=allweight*gausslegendreint64(kernel72)
               
                kernel[i,s,j,m,7,3]+=allweight*gausslegendreint64(kernel73)
               
                kernel[i,s,j,m,7,4]+=allweight*gausslegendreint64(kernel74)
               
                kernel[i,s,j,m,7,5]+=allweight*gausslegendreint64(kernel75)
               
                kernel[i,s,j,m,7,6]+=allweight*gausslegendreint64(kernel76)
               
                kernel[i,s,j,m,7,7]+=allweight*gausslegendreint64(kernel77)
               
                kernel[i,s,j,m,7,8]+=allweight*gausslegendreint64(kernel78)
               
                kernel[i,s,j,m,8,1]+=allweight*gausslegendreint64(kernel81)
               
                kernel[i,s,j,m,8,2]+=allweight*gausslegendreint64(kernel82)
               
                kernel[i,s,j,m,8,3]+=allweight*gausslegendreint64(kernel83)
               
                kernel[i,s,j,m,8,4]+=allweight*gausslegendreint64(kernel84)
               
                kernel[i,s,j,m,8,5]+=allweight*gausslegendreint64(kernel85)
               
                kernel[i,s,j,m,8,6]+=allweight*gausslegendreint64(kernel86)
               
                kernel[i,s,j,m,8,7]+=allweight*gausslegendreint64(kernel87)
               
                kernel[i,s,j,m,8,8]+=allweight*gausslegendreint64(kernel88)
            end
        end
    end
end

end # if tailed

# # print("完成kernel计算，用时",round((time()-timetest)*100)/100,"s \n")
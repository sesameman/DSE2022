timetest = time()
# 计算
# i & s 代表外动量 k2
# j & m 代表内动量 q2

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
                allweight = -weightzq*weightq2*q2/(16*pi^3)/((P2/4+q2+sqrt(P2*q2)*zq)*A1in^2+B1in^2)/((P2/4+q2-sqrt(P2*q2)*zq)*A2in^2+B2in^2)*z2^2
                kernel11(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))
                
                kernel12(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-2*(A2in*B1in + A1in*B2in)*P2 + 4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)
                
                kernel13(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-2*sqrt(P2*q2)*zq*(2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))
                
                kernel14(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(-4*A1in*A2in*(-(P2*q2) + P2*q2*zq^2))
                
                kernel21(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(k2^2*((A2in*B1in + A1in*B2in)*P2 + 2*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + k2*(-((A2in*B1in + A1in*B2in)*k2*P2*zk^2) + 2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in + A1in*B2in)*P2 + 6*(-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq) + 2*sqrt(k2)*sqrt(q2)*(-((A2in*B1in + A1in*B2in)*P2) + 4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 2*sqrt(k2*P2)*zk*(-((A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + (A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + sqrt(k2*P2)*zk*((A2in*B1in + A1in*B2in)*sqrt(k2*P2)*zk*(-q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(3*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + 4*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel22(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((k2^2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) - 8*A1in*A2in*P2*q2*zq^2) + sqrt(k2*P2)*zk*(sqrt(k2*P2)*(-4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*zk*(-q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 2*sqrt(k2)*sqrt(q2)*sqrt(P2*q2)*zq*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-4*B1in*B2in + A1in*A2in*(P2 - 8*q2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))) + k2*(4*A1in*A2in*P2*q2^2 + k2*P2*(4*B1in*B2in - A1in*A2in*(P2 + 4*q2))*zk^2 - 2*sqrt(k2)*P2*(-4*B1in*B2in + A1in*A2in*P2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 2*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(4*B1in*B2in + A1in*A2in*(-P2 - 4*q2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 2*P2*q2*zq^2*(-4*B1in*B2in + A1in*A2in*(P2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + q2*(-16*A1in*A2in*P2*q2*zq^2 + P2*(-4*B1in*B2in + A1in*A2in*(P2 - 8*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel23(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((sqrt(P2*q2)*zq*(k2^2*(-4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*sqrt(P2*q2)*zq + sqrt(k2*P2)*zk*(-2*A1in*A2in*sqrt(k2*P2)*q2*sqrt(P2*q2)*zk*zq - 4*k2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(12*B1in*B2in*q2 + A1in*A2in*(3*q2*(P2 + 4*q2) + 8*sqrt(k2*P2)*sqrt(P2*q2)*zk*zq - 4*P2*q2*zq^2))) + k2*(-2*A1in*A2in*k2*P2*sqrt(P2*q2)*zk^2*zq + sqrt(P2*q2)*zq*(-(q2*(12*B1in*B2in + A1in*A2in*(P2 + 12*q2))) + 4*A1in*A2in*P2*q2*zq^2 + 16*sqrt(k2)*sqrt(q2)*(B1in*B2in + A1in*A2in*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2*P2)*zk*(-4*A1in*A2in*P2*q2*zq^2 + sqrt(k2)*sqrt(q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel24(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(-2*(A2in*B1in + A1in*B2in)*k2^2*(-(P2*q2) + P2*q2*zq^2) + sqrt(k2*P2)*zk*(sqrt(k2*P2)*q2*zk*(-2*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*(-2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(2*sqrt(k2*P2)*zk - sqrt(P2*q2)*zq) + q2*(3*(-(A2in*B1in) + A1in*B2in)*P2 + 8*(A2in*B1in + A1in*B2in)*sqrt(k2*P2)*zk + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + 4*k2*q2*((A2in*B1in - A1in*B2in)*P2 - 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))^2) + k2*(k2*P2*zk^2*(-2*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) - 2*((A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(-(P2*q2) + P2*q2*zq^2) + 2*sqrt(k2)*sqrt(q2)*(4*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + P2*(-2*(A2in*B1in + A1in*B2in)*q2 + (-(A2in*B1in) + A1in*B2in)*sqrt(P2*q2)*zq))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2*P2)*zk*(2*sqrt(P2*q2)*zq*(-2*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel31(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-4*((A2in*B1in + A1in*B2in)*(k2*P2)^1.5*zk^3 - 2*k2*P2*zk^2*((-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + (A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2*P2)*zk*(-(k2*((A2in*B1in + A1in*B2in)*P2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)) + sqrt(P2*q2)*zq*(3*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((A2in*B1in + A1in*B2in)*P2 + 4*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*(k2*(2*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + 3*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(3*(-(A2in*B1in) + A1in*B2in)*q2 + (A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq + 4*(A2in*B1in - A1in*B2in)*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel32(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*((k2*P2)^1.5*(4*B1in*B2in - A1in*A2in*(P2 + 4*q2))*zk^3 + sqrt(k2*P2)*zk*(P2*(4*B1in*B2in - A1in*A2in*(P2 - 8*q2))*q2*zq^2 + k2*(A1in*A2in*P2^2 + P2*(-4*B1in*B2in + 4*A1in*A2in*q2) + 4*A1in*A2in*P2*q2*zq^2) + sqrt(k2)*sqrt(q2)*(-(A1in*A2in*P2^2) + 4*P2*(B1in*B2in - A1in*A2in*q2) - 16*A1in*A2in*P2*q2*zq^2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*k2*P2*sqrt(P2*q2)*zk^2*zq*(-4*B1in*B2in + A1in*A2in*(P2 + 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + P2*sqrt(P2*q2)*zq*(k2*(4*B1in*B2in - A1in*A2in*(P2 - 4*q2 + 12*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-4*B1in*B2in + A1in*A2in*(P2 - 8*q2 + 16*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel33(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((sqrt(P2*q2)*zq*(-4*A1in*A2in*(k2*P2)^1.5*sqrt(P2*q2)*zk^3*zq + sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(12*B1in*B2in*q2 + k2*(4*B1in*B2in + A1in*A2in*(5*P2 + 4*q2)) + A1in*A2in*(3*q2*(P2 + 4*q2) - 4*P2*q2*zq^2) - 8*sqrt(k2)*sqrt(q2)*(2*B1in*B2in + A1in*A2in*(P2 + 2*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*k2*P2*zk^2*(-(q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))) + 4*A1in*A2in*P2*q2*zq^2 + sqrt(k2)*sqrt(q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + P2*(k2*(2*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) - 4*A1in*A2in*P2*q2*zq^2 - 3*sqrt(k2)*sqrt(q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(-3*q2*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2)) + 4*A1in*A2in*P2*q2*zq^2 + 4*sqrt(k2)*sqrt(q2)*(4*B1in*B2in + A1in*A2in*(P2 + 4*q2))*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel34(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(k2*P2)^1.5*zk^3*(-2*(A2in*B1in + A1in*B2in)*q2 + (A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq) + 4*k2*P2*zk^2*(2*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + q2*((A2in*B1in - A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + 2*P2*(k2*(2*(A2in*B1in - A1in*B2in)*(-(P2*q2) + P2*q2*zq^2) + 3*sqrt(k2)*sqrt(q2)*((A2in*B1in - A1in*B2in)*P2 - 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))*(2*(-(A2in*B1in) + A1in*B2in)*P2*q2*zq^2 + q2*(3*(A2in*B1in - A1in*B2in)*P2 - 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq) + 4*sqrt(k2)*sqrt(q2)*((-(A2in*B1in) + A1in*B2in)*P2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + 2*sqrt(k2*P2)*zk*(sqrt(P2*q2)*zq*(2*(A2in*B1in - A1in*B2in)*P2*q2*zq^2 + 2*(A2in*B1in + A1in*B2in)*sqrt(P2*q2)*zq*(k2 - 4*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) - 3*(A2in*B1in - A1in*B2in)*P2*(k2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))) + q2*(3*(-(A2in*B1in) + A1in*B2in)*P2*sqrt(P2*q2)*zq + 2*(A2in*B1in + A1in*B2in)*P2*q2*zq^2 + 4*(A2in*B1in + A1in*B2in)*P2*(k2 - sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *sqrt(k2*P2)*zk*(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel41(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*A1in*A2in*(-2*k2*P2*q2*zk^2 - sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) + sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(sqrt(k2*P2)*zk - 2*sqrt(P2*q2)*zq) - P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel42(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((4*(A2in*B1in + A1in*B2in)*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel43(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((2*(A2in*B1in - A1in*B2in)*sqrt(P2*q2)*zq*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                
                kernel44(y) = D(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*(((4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(2*k2*P2*q2*zk^2 + sqrt(k2)*P2*q2^1.5*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)) - sqrt(k2*P2)*sqrt(P2*q2)*zk*zq*(q2 + 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))) + k2*(sqrt(P2*q2)*zq*(-(sqrt(k2*P2)*zk) + 2*sqrt(P2*q2)*zq) + P2*(-2*q2 + sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))))/(3. *(-(k2*P2) + k2*P2*zk^2)*(k2 + q2 - 2*sqrt(k2)*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))))
                kernel[i,s,j,m,1,1] = allweight*gausslegendreint64(kernel11)
                kernel[i,s,j,m,1,2] = allweight*gausslegendreint64(kernel12)
                kernel[i,s,j,m,1,3] = allweight*gausslegendreint64(kernel13)
                kernel[i,s,j,m,1,4] = allweight*gausslegendreint64(kernel14)
                kernel[i,s,j,m,2,1] = allweight*gausslegendreint64(kernel21)
                kernel[i,s,j,m,2,2] = allweight*gausslegendreint64(kernel22)
                kernel[i,s,j,m,2,3] = allweight*gausslegendreint64(kernel23)
                kernel[i,s,j,m,2,4] = allweight*gausslegendreint64(kernel24)
                kernel[i,s,j,m,3,1] = allweight*gausslegendreint64(kernel31)
                kernel[i,s,j,m,3,2] = allweight*gausslegendreint64(kernel32)
                kernel[i,s,j,m,3,3] = allweight*gausslegendreint64(kernel33)
                kernel[i,s,j,m,3,4] = allweight*gausslegendreint64(kernel34)
                kernel[i,s,j,m,4,1] = allweight*gausslegendreint64(kernel41)
                kernel[i,s,j,m,4,2] = allweight*gausslegendreint64(kernel42)
                kernel[i,s,j,m,4,3] = allweight*gausslegendreint64(kernel43)
                kernel[i,s,j,m,4,4] = allweight*gausslegendreint64(kernel44)
            end
        end
    end
end

# 加一个尾巴
if dataset["mesonBSE"]["epsilon"]!=0
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
                kernel11(y) = 0

                kernel12(y) = 0
                
                kernel13(y) = 0
                
                kernel14(y) = 0
                
                kernel21(y) = 0
                
                kernel22(y) = 0
                
                kernel23(y) = 0
                
                kernel24(y) = 0
                
                kernel31(y) = 0
                
                kernel32(y) = 0
                
                kernel33(y) = 0
                
                kernel34(y) = 0
                
                kernel41(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((-32*A1in*A2in*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *k2*P2*(-1 + zk^2)))

                kernel42(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((32*(A2in*B1in + A1in*B2in)*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *k2*P2*(-1 + zk^2)))

                kernel43(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((16*(A2in*B1in - A1in*B2in)*sqrt(q2)*zq*(-(sqrt(k2*P2)*sqrt(q2)*zk*zq) + sqrt(k2)*sqrt(P2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *k2*(-1 + zk^2)))

                kernel44(y) = D_infrared(k2+q2-2*sqrt(k2*q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2)))*((8*(4*B1in*B2in + A1in*A2in*(P2 - 4*q2))*(-(sqrt(k2*P2)*sqrt(P2*q2)*zk*zq) + sqrt(k2)*P2*sqrt(q2)*(zk*zq + y*sqrt(1 - zk^2)*sqrt(1 - zq^2))))/(9. *k2*P2*(-1 + zk^2)))
                kernel[i,s,j,m,1,1] += allweight*gausslegendreint64(kernel11)
                kernel[i,s,j,m,1,2] += allweight*gausslegendreint64(kernel12)
                kernel[i,s,j,m,1,3] += allweight*gausslegendreint64(kernel13)
                kernel[i,s,j,m,1,4] += allweight*gausslegendreint64(kernel14)
                kernel[i,s,j,m,2,1] += allweight*gausslegendreint64(kernel21)
                kernel[i,s,j,m,2,2] += allweight*gausslegendreint64(kernel22)
                kernel[i,s,j,m,2,3] += allweight*gausslegendreint64(kernel23)
                kernel[i,s,j,m,2,4] += allweight*gausslegendreint64(kernel24)
                kernel[i,s,j,m,3,1] += allweight*gausslegendreint64(kernel31)
                kernel[i,s,j,m,3,2] += allweight*gausslegendreint64(kernel32)
                kernel[i,s,j,m,3,3] += allweight*gausslegendreint64(kernel33)
                kernel[i,s,j,m,3,4] += allweight*gausslegendreint64(kernel34)
                kernel[i,s,j,m,4,1] += allweight*gausslegendreint64(kernel41)
                kernel[i,s,j,m,4,2] += allweight*gausslegendreint64(kernel42)
                kernel[i,s,j,m,4,3] += allweight*gausslegendreint64(kernel43)
                kernel[i,s,j,m,4,4] += allweight*gausslegendreint64(kernel44)
            end
        end
    end
end

end # tailed
# print("完成kernel计算，用时",round((time()-timetest)*100)/100,"s \n")
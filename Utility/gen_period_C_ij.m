function C_ij = gen_period_C_ij(U_turn_vel_x, U_turn_vel_y, T, tau, cyctime)
    C_ij = [];
    for k = 1:size(U_turn_vel_x, 2)
        for p = 1:size(U_turn_vel_x, 2)
            if k ~= p
                C_ij_tmp = [];
                tauTmp_list = floor([-tau:tau]);
                for tauTmp = tauTmp_list
                    tau_idx = find(tauTmp == tauTmp_list);
                    v_ij_tmp = 0;
                    cnt = 0;
                    if tauTmp <= 0
                        for t = abs(tauTmp) + T - tau:T
                            cnt = cnt + 1;
                            focal_vel = [U_turn_vel_x(t,k);U_turn_vel_y(t,k)];
                            neighbor_vel = [U_turn_vel_x(t + tauTmp,p); U_turn_vel_y(t + tauTmp,p)];
                            v_ij_tmp = v_ij_tmp + dot(focal_vel, neighbor_vel);
                        end
                    else
                        for t = T - tau:T - tauTmp
                            cnt = cnt + 1;
                            focal_vel = [U_turn_vel_x(t,k);U_turn_vel_y(t,k)];
                            neighbor_vel = [U_turn_vel_x(t + tauTmp,p); U_turn_vel_y(t + tauTmp,p)];
                            v_ij_tmp = v_ij_tmp + dot(focal_vel, neighbor_vel);
                        end
                    end
                    C_ij_tmp(tau_idx) = v_ij_tmp/cnt;
                end
                max_tau_ij = tauTmp_list(find(max(C_ij_tmp) == C_ij_tmp, 1)) * cyctime;
                C_ij(p,k) = max_tau_ij;
            else
                C_ij(p,k) = nan;
            end
        end
    end
end
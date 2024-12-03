classdef shift
    methods(Static)
        function [time_shifted, thrust_shifted, mass_shifted] = start(time, thrust, mass, tx_shift)
            % shift theoretical data to experimental start time
            time_shifted = time + tx_shift;
        
            if tx_shift < 0
                idx = find(time_shifted >= 0);
                time_shifted = time_shifted(idx);
                thrust_shifted = thrust(idx);
                mass_shifted = mass(idx);
            elseif tx_shift > 0
                del_t = mean(diff(time_shifted));
                time_pre = (0 : del_t : time_shifted(1)-del_t)';
                if ~isempty(time_pre)
                    time_shifted = [time_pre; time_shifted];
                    thrust_shifted = [zeros(size(time_pre)); thrust];
                    mass_shifted = [mass(1)*ones(size(time_pre)); mass];
                else
                    thrust_shifted = thrust;
                    mass_shifted = mass;
                end
            else
                thrust_shifted = thrust;
                mass_shifted = mass;
            end
        end
    end
end